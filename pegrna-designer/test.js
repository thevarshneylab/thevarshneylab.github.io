/*
 * Danio pegRNA Designer — test suite. Run with: node test.js
 * No dependencies. Exits non-zero on failure, so it works as a CI step.
 */
"use strict";

const PEG = require("./pegrna.js");

let passed = 0;
const failures = [];

function check(name, condition, detail) {
  if (condition) { passed++; return; }
  failures.push(name + (detail ? "\n    " + detail : ""));
}
function equal(name, actual, expected) {
  check(name, JSON.stringify(actual) === JSON.stringify(expected),
    "expected " + JSON.stringify(expected) + ", got " + JSON.stringify(actual));
}
function near(name, actual, expected, tol) {
  check(name, Math.abs(actual - expected) <= tol,
    "expected " + expected + " ± " + tol + ", got " + actual);
}

/* ------------------------------------------------------------------ *
 * Sequence basics
 * ------------------------------------------------------------------ */

equal("reverse complement", PEG.rc("ATGC"), "GCAT");
equal("rc is an involution", PEG.rc(PEG.rc("ACGTTGCA")), "ACGTTGCA");
equal("clean strips non-bases", PEG.clean("acg 12 t\nN"), "ACGTN");
check("NGG matches TGG", PEG.pamMatch("TGG", 0, "NGG"));
check("NGG rejects TGA", !PEG.pamMatch("TGA", 0, "NGG"));
check("NRN matches TAC", PEG.pamMatch("TAC", 0, "NRN"));

/* ------------------------------------------------------------------ *
 * HEK3 +1 CTT insertion — the published reference design.
 *
 * Protospacer GGCCCAGACTGAGCACGTGA / PAM TGG. With a 13 nt PBS and a
 * 13 nt RTT, Anzalone et al. 2019 give the 3' extension as
 * TCTGCCATCAAAG + CGTGCTCAGTCTG. If the nick position or the reverse
 * complement were off by one, neither half would match.
 * ------------------------------------------------------------------ */

const HEK3 =
  "TCTGGCTCTGACTACAGCTTGGAAGCTGGAGGCCCAGACTGAGCACG[-/CTT]" +
  "TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGGGCCCAGGCTTAGG";

const hek3 = PEG.design(HEK3, { pbsMin: 10, pbsMax: 15, rttMin: 10, rttMax: 20 });
check("HEK3 target parses", !hek3.error, hek3.error);
equal("HEK3 edit is an insertion", hek3.target.kind, "insertion");

const hek3site = hek3.sites.find(s => s.protospacer === "GGCCCAGACTGAGCACGTGA");
check("HEK3 protospacer found", !!hek3site);

if (hek3site) {
  equal("HEK3 PAM", hek3site.pamSeq, "TGG");
  equal("HEK3 is on the + strand", hek3site.strand, "+");
  equal("HEK3 nick sits at the edit", hek3site.nickToEdit, 0);
  check("HEK3 insertion disrupts the PAM", hek3site.pamChanged);

  const p = hek3site.pegs.find(x => x.pbsLen === 13 && x.rttLen === 13);
  check("HEK3 13/13 design exists", !!p);
  equal("HEK3 PBS (13 nt)", p.pbs, "CGTGCTCAGTCTG");
  equal("HEK3 RTT (13 nt)", p.rtt, "TCTGCCATCAAAG");
  equal("HEK3 3' extension is RTT then PBS", p.ext, "TCTGCCATCAAAG" + "CGTGCTCAGTCTG");
  equal("extension length is RTT + PBS", p.ext.length, p.rttLen + p.pbsLen);
  equal("3' homology excludes the insert", p.homology3, 13 - 0 - 3);
}

/* Every design, on either strand, must keep extension = RTT + PBS. */
let extOk = true;
hek3.sites.forEach(s => s.pegs.forEach(p => {
  if (p.ext !== p.rtt + p.pbs || p.ext.length !== p.rttLen + p.pbsLen) extOk = false;
}));
check("all extensions are RTT + PBS", extOk);

/* ------------------------------------------------------------------ *
 * Substitutions, deletions, and the reverse strand
 * ------------------------------------------------------------------ */

const SUB =
  "GCTGTGATCCTCACAGTCTGGAGAAACACTGTTTCTGCAGGACTGCAGTTTCAGAGGCTCTGGAAGCACTTTGGAAGCTGGAGG" +
  "CCCAGACTGAGCACGTGATGGCAGAGG[A/G]AAGGAAGCCCTGCTTCCTCCAGGGCCCAGGCTTAGGCTTCAGTCAGGCTCCAGG" +
  "AGAAGCCTGGTAGACAGTCAGTGGCCTCAAAGTTTGCAGACTGA";

const sub = PEG.design(SUB, {});
equal("substitution is typed", sub.target.kind, "substitution");
check("substitution finds candidates", sub.sites.length > 0);
check("both strands are searched", sub.sites.some(s => s.strand === "+") && sub.sites.some(s => s.strand === "−"));

/* A minus-strand design must be built from the minus-strand view: its
   spacer should not appear in the forward reference. */
const minus = sub.sites.find(s => s.strand === "−");
if (minus) {
  check("minus-strand spacer is on the minus strand",
    sub.target.ref.indexOf(minus.protospacer) === -1 &&
    PEG.rc(sub.target.ref).indexOf(minus.protospacer) !== -1);
}

const DEL = SUB.replace("[A/G]", "[AAGG/-]");
const del = PEG.design(DEL, {});
equal("deletion is typed", del.target.kind, "deletion");
check("deletion yields designs", del.sites.length > 0);
if (del.sites.length) {
  const s = del.sites[0], p = s.best;
  /* With nothing inserted, RTT homology past the edit is rttLen - nickToEdit. */
  equal("deletion 3' homology", p.homology3, p.rttLen - s.nickToEdit);
}

/* ------------------------------------------------------------------ *
 * Parser rejects what it cannot design
 * ------------------------------------------------------------------ */

check("missing bracket is rejected", !!PEG.parseTarget("ACGTACGTACGT").error);
check("two brackets are rejected", !!PEG.parseTarget("A".repeat(40) + "[C/T]" + "A".repeat(40) + "[G/A]" + "A".repeat(40)).error);
check("empty bracket is rejected", !!PEG.parseTarget("A".repeat(40) + "[-/-]" + "A".repeat(40)).error);
check("short flank is rejected", !!PEG.parseTarget("ACGT[C/T]" + "A".repeat(40)).error);

/* ------------------------------------------------------------------ *
 * Melting temperature — the basis of Tm-matched PBS design
 * (Ponnienselvan et al. 2023, NAR 51:6966)
 * ------------------------------------------------------------------ */

const PBS13 = "CGTGCTCAGTCTG";            /* HEK3, GC-rich */
let monotonic = true;
for (let n = 6; n < 13; n++) {
  const a = PEG.tmHybrid(PBS13.slice(-n), 2e-5);
  const b = PEG.tmHybrid(PBS13.slice(-(n + 1)), 2e-5);
  if (!(b > a)) monotonic = false;
}
check("hybrid Tm rises with PBS length", monotonic);
near("GC-rich 7 nt PBS melts near 28.5 C", PEG.tmHybrid(PBS13.slice(-7), 2e-5), 28.5, 2);
near("GC-rich 8 nt PBS melts near 37 C", PEG.tmHybrid(PBS13.slice(-8), 2e-5), 37, 2.5);

/* A/T-rich sequence must need a longer PBS to reach the same Tm than a
   G/C-rich one — the observation the Tm rule is built on. */
const AT = "TTATTAAGTTATA";
check("A/T-rich PBS is colder at equal length",
  PEG.tmHybrid(AT.slice(-10), 2e-5) < PEG.tmHybrid(PBS13.slice(-10), 2e-5));

check("hybrid Tm needs at least a dinucleotide", isNaN(PEG.tmHybrid("A", 2e-5)));
check("DNA model also runs", !isNaN(PEG.tmDna(PBS13, 0.25e-6, 0.05)));

if (hek3site) {
  const m = PEG.tmMatched(hek3site, Object.assign(PEG.defaults(), { targetTm: 28.5, minHomology: 7 }));
  check("a Tm-matched PBS is chosen", !!m);
  check("Tm-matched PBS is within range", m.pbsLen >= 10 && m.pbsLen <= 15);
  /* It must genuinely be the closest length available, not merely a short one. */
  let closest = true;
  hek3site.pegs.forEach(p => {
    if (p.homology3 < 7 || isNaN(p.pbsTm)) return;
    if (Math.abs(p.pbsTm - 28.5) < Math.abs(m.pbsTm - 28.5) - 1e-9) closest = false;
  });
  check("Tm-matched PBS is the closest available", closest);
}

/* ------------------------------------------------------------------ *
 * PBS-spacer auto-inhibition
 * ------------------------------------------------------------------ */

let autoOk = true;
sub.sites.forEach(s => s.pegs.forEach(p => {
  if (p.autoLen !== p.pbsLen + Math.min(3, s.nickToEdit)) autoOk = false;
}));
check("auto-inhibitory duplex = PBS + run-on into unedited RTT bases", autoOk);
if (hek3site) {
  const p = hek3site.pegs.find(x => x.pbsLen === 13 && x.rttLen === 13);
  /* The HEK3 insert sits at the nick, so nothing runs on past the PBS. */
  equal("no run-on when the edit starts at the nick", p.autoLen, 13);
}

/* ------------------------------------------------------------------ *
 * Nicking guides
 * ------------------------------------------------------------------ */

if (sub.sites.length) {
  const s = sub.sites[0];
  check("nicking guides are found", s.nicking.length > 0);
  const refRc = PEG.rc(s.view);
  let pe3bOk = true;
  s.nicking.forEach(g => {
    if (g.type === "PE3b" && refRc.indexOf(g.spacer + g.pam) !== -1) pe3bOk = false;
  });
  check("PE3b spacers are absent from the reference", pe3bOk);
  check("preferred window is 40-90 nt",
    s.nicking.every(g => !g.preferred || (g.dist >= 40 && g.dist <= 90)));
}

/* ------------------------------------------------------------------ *
 * Assembly
 * ------------------------------------------------------------------ */

if (hek3site) {
  const p = hek3site.pegs.find(x => x.pbsLen === 13 && x.rttLen === 13);
  const opt = Object.assign(PEG.defaults(), { g5: true, useMotif: true, term: false, scaffold: "std" });
  const a = PEG.assemble(hek3site, p, opt);
  equal("spacer already starts with G, so none is added", a.addG, false);
  equal("spacer is 20 nt", a.spacer.length, 20);
  equal("full pegRNA is the concatenation",
    a.full, a.spacer + PEG.SCAFFOLDS.std + p.ext + PEG.TEVOPREQ1);

  const noMotif = PEG.assemble(hek3site, p, Object.assign({}, opt, { useMotif: false, term: true }));
  equal("terminator is appended", noMotif.full.slice(-7), "TTTTTTT");
  check("motif is omitted when switched off", noMotif.full.indexOf(PEG.TEVOPREQ1) === -1);
}

/* A spacer that does not start with G must gain one for U6. */
const nonG = { protospacer: "ATCCCAGACTGAGCACGTGA" };
const gAdded = PEG.assemble(nonG, { ext: "AAAA" }, Object.assign(PEG.defaults(), { useMotif: false }));
equal("5' G is added for U6", gAdded.addG, true);
equal("spacer becomes 21 nt", gAdded.spacer.length, 21);
equal("original spacer is preserved", gAdded.spacer.slice(1), nonG.protospacer);

/* ------------------------------------------------------------------ *
 * GRCz11 coordinate lookup
 * ------------------------------------------------------------------ */

const fields = { chrom: "7", pos: 12345678, ref: "A", alt: "G", flank: 150 };
const w = PEG.windowBounds(fields);
equal("window is centred on the variant", [w.start, w.end], [12345528, 12345828]);
check("Ensembl URL is well formed",
  PEG.ensemblUrl(fields) ===
  "https://rest.ensembl.org/sequence/region/danio_rerio/7:12345528..12345828:1?content-type=text/x-fasta");

/* A synthetic window whose header declares its own coordinates. */
const WINSEQ = "A".repeat(150) + "A" + "CGTACGTACG".repeat(15);   /* pos 12345678 is an A */
const FASTA = ">chromosome:GRCz11:7:12345528:" + (12345528 + WINSEQ.length - 1) + ":1\n" + WINSEQ;

const built = PEG.buildBracketed(FASTA, fields);
check("bracketed target is built", !built.error, built.error);
if (!built.error) {
  equal("window start is taken from the header", built.start, 12345528);
  check("header is trusted", built.fromHeader);
  check("bracket is placed at the variant", built.target.indexOf("[A/G]") === 150);
  equal("bracketing preserves the sequence",
    built.target.replace("[A/G]", "A"), WINSEQ);
}

/* The reference check is the guard that matters most: a variant called on
   the wrong strand or assembly must fail loudly rather than design silently. */
const wrongRef = PEG.buildBracketed(FASTA, Object.assign({}, fields, { ref: "T" }));
check("reference mismatch is caught", !!wrongRef.error && wrongRef.mismatch);
equal("mismatch reports the real base", wrongRef.found, "A");

const wrongAsm = PEG.buildBracketed(FASTA.replace("GRCz11", "GRCz10"), fields);
check("wrong assembly is rejected", !!wrongAsm.error);

const truncated = PEG.buildBracketed(FASTA.slice(0, FASTA.length - 20), fields);
check("truncated paste is rejected", !!truncated.error);

const revStrand = PEG.buildBracketed(
  ">chromosome:GRCz11:7:12345528:" + (12345528 + WINSEQ.length - 1) + ":-1\n" + WINSEQ, fields);
check("reverse-strand window is rejected", !!revStrand.error);

/* An insertion is placed after the given position, not over it. */
const ins = PEG.buildBracketed(FASTA, Object.assign({}, fields, { ref: "", alt: "CTT" }));
check("insertion is bracketed after the position", !ins.error && ins.target.indexOf("[-/CTT]") === 151);

/* Designs built from a coordinate window must still work end to end. */
if (!built.error) {
  const fromCoord = PEG.design(built.target, {});
  check("coordinate-built target designs", !fromCoord.error);
}

/* ------------------------------------------------------------------ *
 * Report
 * ------------------------------------------------------------------ */

if (failures.length) {
  console.error("\n" + failures.length + " failed, " + passed + " passed\n");
  failures.forEach(f => console.error("  ✗ " + f));
  console.error("");
  process.exit(1);
}
console.log(passed + " checks passed");
