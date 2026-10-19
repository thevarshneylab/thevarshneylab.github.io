/*!
 * Danio pegRNA Designer — prime editing guide design for zebrafish.
 * https://github.com/thevarshneylab/thevarshneylab.github.io
 *
 * Runs in the browser (loaded by index.html) and in Node (required by test.js).
 * Design rules are documented inline against the papers they come from.
 */
(function (root, factory) {
  var api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.PEGRNA = api;
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /* ------------------------------------------------------------------ *
   * Sequence basics
   * ------------------------------------------------------------------ */

  var COMP = { A: "T", C: "G", G: "C", T: "A", N: "N" };

  function rc(s) {
    var o = "";
    for (var i = s.length - 1; i >= 0; i--) o += (COMP[s[i]] || "N");
    return o;
  }
  function clean(s) {
    return String(s || "").toUpperCase().replace(/[^ACGTN]/g, "");
  }
  function gc(s) {
    return s.length ? 100 * (s.match(/[GC]/g) || []).length / s.length : 0;
  }

  var IUPAC = {
    A: "A", C: "C", G: "G", T: "T",
    R: "AG", Y: "CT", S: "CG", W: "AT", K: "GT", M: "AC",
    B: "CGT", D: "AGT", H: "ACT", V: "ACG", N: "ACGT"
  };
  function pamMatch(seq, start, pam) {
    if (start + pam.length > seq.length) return false;
    for (var k = 0; k < pam.length; k++) {
      var set = IUPAC[pam[k]];
      if (!set || set.indexOf(seq[start + k]) === -1) return false;
    }
    return true;
  }

  /* ------------------------------------------------------------------ *
   * Melting temperatures
   *
   * The PBS is RNA annealing to the nicked DNA target strand, so the
   * hybrid model is the one that matters for choosing its length.
   * Parameters: Banerjee et al. 2020, Nucleic Acids Res 48:12042
   * (RNA/DNA nearest neighbours at 100 mM NaCl), keyed by the RNA
   * dinucleotide 5'->3'. The DNA/DNA model (SantaLucia 1998) is kept as
   * a second opinion only.
   * ------------------------------------------------------------------ */

  var NN_HYBRID = {
    AA: [-7.8, -22.9], AC: [-10.1, -27.3], AG: [-9.4, -26.2], AU: [-5.8, -17.5],
    CA: [-9.8, -27.4], CC: [-9.5, -24.8], CG: [-9.0, -24.3], CU: [-6.1, -17.9],
    GA: [-8.6, -22.7], GC: [-10.6, -27.7], GG: [-13.3, -35.7], GU: [-9.3, -25.5],
    UA: [-6.6, -19.7], UC: [-6.5, -16.3], UG: [-8.9, -23.3], UU: [-7.4, -24.3]
  };

  var NN_DNA = {
    AA: [-7.9, -22.2], TT: [-7.9, -22.2], AT: [-7.2, -20.4], TA: [-7.2, -21.3],
    CA: [-8.5, -22.7], TG: [-8.5, -22.7], GT: [-8.4, -22.4], AC: [-8.4, -22.4],
    CT: [-7.8, -21.0], AG: [-7.8, -21.0], GA: [-8.2, -22.2], TC: [-8.2, -22.2],
    CG: [-10.6, -27.2], GC: [-9.8, -24.4], GG: [-8.0, -19.9], CC: [-8.0, -19.9]
  };

  var R_GAS = 1.987;

  /* PBS (RNA, written 5'->3' as DNA letters) annealed to the target strand. */
  function tmHybrid(pbsDna, concM) {
    var r = String(pbsDna).replace(/T/g, "U");
    if (r.length < 2) return NaN;
    var dH = 0, dS = 0;
    for (var i = 0; i < r.length - 1; i++) {
      var p = NN_HYBRID[r.substr(i, 2)];
      if (!p) return NaN;
      dH += p[0]; dS += p[1];
    }
    /* Initiation: dH = 0; dS depends on whether either terminus is a G-C pair. */
    var ends = r[0] + r[r.length - 1];
    dS += /[GC]/.test(ends) ? -4.9 : -7.0;
    return (dH * 1000) / (dS + R_GAS * Math.log((concM || 2e-5) / 4)) - 273.15;
  }

  function tmDna(seq, concM, naM) {
    if (seq.length < 2) return NaN;
    var dH = 0, dS = 0;
    for (var i = 0; i < seq.length - 1; i++) {
      var p = NN_DNA[seq.substr(i, 2)];
      if (!p) return NaN;
      dH += p[0]; dS += p[1];
    }
    function ini(b) { return (b === "G" || b === "C") ? [0.1, -2.8] : [2.3, 4.1]; }
    var i5 = ini(seq[0]), i3 = ini(seq[seq.length - 1]);
    dH += i5[0] + i3[0];
    dS += i5[1] + i3[1];
    dS += 0.368 * (seq.length - 1) * Math.log(naM || 0.05);
    return (dH * 1000) / (dS + R_GAS * Math.log((concM || 0.25e-6) / 4)) - 273.15;
  }

  /* ------------------------------------------------------------------ *
   * Target parsing:  ACGT[C/T]ACGT  ·  [-/CTT] insert  ·  [AGG/-] delete
   * ------------------------------------------------------------------ */

  function parseTarget(raw) {
    var m = /\[([^\]\/]*)\/([^\]\/]*)\]/.exec(raw || "");
    if (!m) {
      return { error: "No edit found. Mark it in the sequence as [ref/alt] — for example [C/T], [-/CTT] to insert, or [AGG/-] to delete." };
    }
    if (/\[[^\]]*\]/.test(String(raw).replace(m[0], ""))) {
      return { error: "More than one bracket in the sequence. Design one edit at a time." };
    }
    var left = clean(String(raw).slice(0, m.index));
    var right = clean(String(raw).slice(m.index + m[0].length));
    var refA = clean(m[1] === "-" ? "" : m[1]);
    var altA = clean(m[2] === "-" ? "" : m[2]);
    if (!refA && !altA) {
      return { error: "The bracket is empty. Put the reference allele before the slash and the edited allele after it." };
    }
    if (left.length < 30) return { error: "Only " + left.length + " nt upstream of the edit. Paste at least 30, ideally 100." };
    if (right.length < 30) return { error: "Only " + right.length + " nt downstream of the edit. Paste at least 30, ideally 100." };
    var kind = !refA ? "insertion"
      : !altA ? "deletion"
      : (refA.length === altA.length ? (refA.length === 1 ? "substitution" : "multi-nt substitution") : "indel");
    return {
      ref: left + refA + right,
      alt: left + altA + right,
      editStart: left.length,
      refA: refA, altA: altA, kind: kind
    };
  }

  /* ------------------------------------------------------------------ *
   * Protospacer search
   *
   * SpCas9 nicks the PAM-containing strand between protospacer positions
   * 17 and 18, i.e. 3 bp 5' of the PAM (Anzalone et al. 2019). The edit
   * must lie 3' of that nick for the RT to template it.
   * ------------------------------------------------------------------ */

  var SPACER_LEN = 20;
  var SEED = 10;

  function scanStrand(ref, editStart, refA, altA, opt, strand) {
    var out = [], pam = opt.pam, sl = SPACER_LEN;
    for (var p = sl; p + pam.length <= ref.length; p++) {
      if (!pamMatch(ref, p, pam)) continue;
      var nick = p - 3;
      if (editStart < nick) continue;
      var protospacer = ref.slice(p - sl, p);
      if (/N/.test(protospacer)) continue;
      var nickToEdit = editStart - nick;
      if (nickToEdit > opt.maxNickToEdit) continue;

      var downstreamEdited = ref.slice(nick, editStart) + altA + ref.slice(editStart + refA.length);
      var editedSite = ref.slice(p - sl, nick) + downstreamEdited.slice(0, sl - (nick - (p - sl)) + pam.length);
      var refSite = ref.slice(p - sl, p + pam.length);

      out.push({
        strand: strand,
        protospacer: protospacer,
        pamSeq: ref.slice(p, p + pam.length),
        pamStart: p,
        nick: nick,
        nickToEdit: nickToEdit,
        downstreamEdited: downstreamEdited,
        pamChanged: editedSite.slice(sl, sl + pam.length) !== refSite.slice(sl, sl + pam.length),
        seedChanged: editedSite.slice(sl - SEED, sl) !== refSite.slice(sl - SEED, sl)
      });
    }
    return out;
  }

  /* ------------------------------------------------------------------ *
   * One pegRNA at a given PBS / RTT length
   *
   * The 3' extension reads 5'->3' as [RTT][PBS], and both come from one
   * reverse complement spanning the nick:
   *     extension = rc( PBS region + edited RTT region )
   * ------------------------------------------------------------------ */

  function buildPeg(site, ref, opt, pbsLen, rttLen) {
    if (site.nick - pbsLen < 0) return null;
    if (rttLen > site.downstreamEdited.length) return null;

    var pbsRegion = ref.slice(site.nick - pbsLen, site.nick);
    var rttRegion = site.downstreamEdited.slice(0, rttLen);
    var ext = rc(pbsRegion + rttRegion);
    var rtt = ext.slice(0, rttLen);
    var pbs = ext.slice(rttLen);
    var homology3 = rttLen - site.nickToEdit - site.altALen;

    /* PBS-spacer auto-inhibition (Ponnienselvan et al. 2023, NAR 51:6966):
       the PBS is complementary to the 3' end of its own spacer, and the
       pairing runs on into the first RTT-templated bases wherever those
       still match the target — up to three. */
    var runOn = Math.min(3, site.nickToEdit);
    var autoLen = pbsLen + runOn;
    var autoSeq = rc(ref.slice(site.nick - pbsLen, site.nick + runOn));

    var pbsTm = opt.tmModel === "dna"
      ? tmDna(pbsRegion, 0.25e-6, 0.05)
      : tmHybrid(pbs, opt.oligoConc || 2e-5);

    var flags = [];
    if (site.pamChanged) flags.push({ t: "good", m: "The edit changes the PAM — the editor cannot re-nick the corrected allele." });
    else if (site.seedChanged) flags.push({ t: "good", m: "The edit falls in the seed region, which disfavours re-nicking." });
    if (autoLen >= 17) flags.push({ t: "warn", m: "Long PBS–spacer duplex (" + autoLen + " bp) — this auto-inhibitory pairing competes with R-loop formation. A shorter PBS relieves it." });
    if (ext[0] === "C") flags.push({ t: "warn", m: "3′ extension starts with C, which lowers editing. Shift the RTT by a nucleotide." });
    if (/T{4,}/.test(site.protospacer) || /T{4,}/.test(ext)) flags.push({ t: "warn", m: "Poly-T run of 4 or more — Pol III can terminate early. A silent edit in the RTT can break it up." });
    if (site.protospacer[0] !== "G") flags.push({ t: "note", m: "Spacer starts with " + site.protospacer[0] + ". Add a 5′ G for U6 transcription, or GG for T7." });
    if (homology3 < opt.minHomology) flags.push({ t: "bad", m: "Only " + homology3 + " nt of homology 3′ of the edit." });

    var score = 0;
    score += site.nickToEdit <= 10 ? 30 : site.nickToEdit <= 20 ? 24 : site.nickToEdit <= 30 ? 16 : 8;
    score += Math.max(0, 20 - Math.abs((isNaN(pbsTm) ? 0 : pbsTm) - opt.targetTm) * 1.6);
    score += autoLen >= 18 ? -8 : autoLen >= 16 ? -4 : 0;
    score += homology3 >= 15 ? 15 : homology3 >= 10 ? 13 : homology3 >= 7 ? 10 : homology3 >= opt.minHomology ? 5 : -25;
    score += site.pamChanged ? 15 : site.seedChanged ? 9 : 0;
    score += (rttLen >= 10 && rttLen <= 30) ? 8 : rttLen <= 40 ? 4 : 0;
    score += (gc(pbsRegion) >= 30 && gc(pbsRegion) <= 70) ? 7 : 2;
    if (ext[0] === "C") score -= 10;
    if (/T{4,}/.test(site.protospacer) || /T{4,}/.test(ext)) score -= 15;
    if (homology3 < opt.minHomology) score -= 40;

    return {
      pbsLen: pbsLen, rttLen: rttLen, pbs: pbs, rtt: rtt, ext: ext,
      pbsRegion: pbsRegion, rttRegion: rttRegion,
      homology3: homology3, pbsGC: gc(pbsRegion), pbsTm: pbsTm,
      autoLen: autoLen, autoGC: gc(autoSeq),
      flags: flags,
      score: Math.max(0, Math.min(100, Math.round(score)))
    };
  }

  /* ------------------------------------------------------------------ *
   * Nicking sgRNAs for PE3 / PE3b
   * ------------------------------------------------------------------ */

  function nickingGuides(view, viewAlt, site, opt) {
    var list = [], sl = SPACER_LEN, pam = opt.pam, refRc = rc(view);

    function search(seq, tag) {
      for (var p = sl; p + pam.length <= seq.length; p++) {
        if (!pamMatch(seq, p, pam)) continue;
        var spacer = seq.slice(p - sl, p);
        if (/N/.test(spacer)) continue;
        var nickInView = seq.length - (p - 3);
        var d = nickInView - site.nick;
        if (Math.abs(d) < 10 || Math.abs(d) > 140) continue;
        list.push({
          spacer: spacer, pam: seq.slice(p, p + pam.length),
          dist: Math.abs(d), dir: d >= 0 ? "downstream" : "upstream",
          type: tag, preferred: Math.abs(d) >= 40 && Math.abs(d) <= 90
        });
      }
    }

    search(refRc, "PE3");
    var before = list.length;
    search(rc(viewAlt), "PE3b");
    /* A PE3b spacer must not exist in the reference — otherwise it nicks before the edit. */
    for (var i = before; i < list.length; i++) {
      if (refRc.indexOf(list[i].spacer + list[i].pam) !== -1) list[i].drop = true;
    }

    return list.filter(function (g) { return !g.drop; })
      .filter(function (g, i, a) { return a.findIndex(function (x) { return x.spacer === g.spacer; }) === i; })
      .sort(function (a, b) {
        return (b.type === "PE3b") - (a.type === "PE3b") || (b.preferred - a.preferred) || a.dist - b.dist;
      })
      .slice(0, 14);
  }

  /* ------------------------------------------------------------------ *
   * Tm-matched PBS (Ponnienselvan et al. 2023)
   * ------------------------------------------------------------------ */

  function tmMatched(site, opt) {
    var best = null;
    for (var i = 0; i < site.pegs.length; i++) {
      var p = site.pegs[i];
      if (isNaN(p.pbsTm)) continue;
      if (p.homology3 < opt.minHomology) continue;
      var d = Math.abs(p.pbsTm - opt.targetTm);
      if (!best || d < best.d) best = { d: d, pbsLen: p.pbsLen, pbsTm: p.pbsTm };
    }
    return best;
  }

  /* ------------------------------------------------------------------ *
   * Assembly
   * ------------------------------------------------------------------ */

  var SCAFFOLDS = {
    std: "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    fe: "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
  };
  var TEVOPREQ1 = "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA";

  function assemble(site, peg, opt) {
    var addG = opt.g5 && site.protospacer[0] !== "G";
    var spacer = (addG ? "G" : "") + site.protospacer;
    var motif = opt.useMotif ? clean(opt.motif) : "";
    var term = opt.term ? "TTTTTTT" : "";
    var scaffold = SCAFFOLDS[opt.scaffold] || SCAFFOLDS.std;
    return {
      spacer: spacer, addG: addG, scaffold: scaffold, motif: motif, term: term,
      full: spacer + scaffold + peg.ext + motif + term
    };
  }

  /* ------------------------------------------------------------------ *
   * Main entry point
   * ------------------------------------------------------------------ */

  function defaults() {
    return {
      pam: "NGG",
      pbsMin: 6, pbsMax: 16,
      rttMin: 10, rttMax: 34,
      minHomology: 7,
      maxNickToEdit: 30,
      targetTm: 28.5,          /* zebrafish embryos are raised at 28.5 C */
      tmModel: "hybrid",
      oligoConc: 2e-5,
      g5: true, term: false,
      useMotif: true, motif: TEVOPREQ1,
      scaffold: "std"
    };
  }

  function design(raw, options) {
    var opt = Object.assign(defaults(), options || {});
    var t = parseTarget(raw);
    if (t.error) return t;

    var L = t.ref.length, sites = [];

    scanStrand(t.ref, t.editStart, t.refA, t.altA, opt, "+").forEach(function (s) {
      s.view = t.ref; s.viewAlt = t.alt; s.altALen = t.altA.length;
      s.genomicNick = s.nick;
      sites.push(s);
    });

    var rRef = rc(t.ref), rAlt = rc(t.alt);
    scanStrand(rRef, L - (t.editStart + t.refA.length), rc(t.refA), rc(t.altA), opt, "−").forEach(function (s) {
      s.view = rRef; s.viewAlt = rAlt; s.altALen = t.altA.length;
      s.genomicNick = L - s.nick;
      sites.push(s);
    });

    var results = [];
    sites.forEach(function (s) {
      var pegs = [];
      for (var b = opt.pbsMin; b <= opt.pbsMax; b++) {
        for (var r = opt.rttMin; r <= opt.rttMax; r++) {
          var peg = buildPeg(s, s.view, opt, b, r);
          if (peg) pegs.push(peg);
        }
      }
      var ok = pegs.filter(function (p) { return p.homology3 >= opt.minHomology; });
      if (!ok.length) return;
      s.pegs = pegs;
      s.best = ok.slice().sort(function (a, b) { return b.score - a.score; })[0];
      s.nicking = nickingGuides(s.view, s.viewAlt, s, opt);
      results.push(s);
    });

    results.sort(function (a, b) { return b.best.score - a.best.score; });
    return { target: t, sites: results, opt: opt };
  }

  /* ------------------------------------------------------------------ *
   * GRCz11 coordinate lookup
   * ------------------------------------------------------------------ */

  function windowBounds(f) {
    var span = Math.max(f.ref.length, 1);
    return { start: Math.max(1, f.pos - f.flank), end: f.pos + span - 1 + f.flank };
  }

  function ensemblUrl(f) {
    var w = windowBounds(f);
    return "https://rest.ensembl.org/sequence/region/danio_rerio/" +
      encodeURIComponent(f.chrom) + ":" + w.start + ".." + w.end +
      ":1?content-type=text/x-fasta";
  }

  function parseFasta(raw) {
    var lines = String(raw || "").split(/\r?\n/), header = null, seq = "";
    for (var i = 0; i < lines.length; i++) {
      if (/^[>;]/.test(lines[i])) { if (!header) header = lines[i].trim(); continue; }
      seq += lines[i];
    }
    return { header: header, seq: clean(seq) };
  }

  /* Turn a fetched window plus a variant call into a bracketed target,
     verifying the reference allele rather than trusting the coordinate. */
  function buildBracketed(fastaText, f) {
    if (!f.chrom || isNaN(f.pos)) return { error: "Enter a chromosome and a position first." };
    if (!f.ref && !f.alt) return { error: "Give a reference allele, an alternate allele, or both. Use - for an empty one." };

    var p = parseFasta(fastaText);
    if (p.seq.length < 80) return { error: "Paste at least 80 nt of sequence." };

    var start = windowBounds(f).start, assembly = "GRCz11", strand = 1, hdrEnd = null, fromHeader = false;
    if (p.header) {
      var h = p.header.replace(/^[>;]\s*/, "").split(":");
      var nums = h.filter(function (x) { return /^-?\d+$/.test(x); }).map(Number);
      var asmIdx = -1;
      for (var i = 0; i < h.length; i++) if (/^GRCz|^danRer/i.test(h[i])) { asmIdx = i; break; }
      if (asmIdx >= 0) assembly = h[asmIdx];
      if (nums.length >= 3) {
        start = nums[nums.length - 3];
        hdrEnd = nums[nums.length - 2];
        strand = nums[nums.length - 1];
        fromHeader = true;
      } else {
        var m = /(\d+)[-–.]{1,2}(\d+)/.exec(p.header);
        if (m) start = parseInt(m[1], 10);
      }
    }

    if (strand === -1) return { error: "That window was fetched on the reverse strand. Re-fetch with :1 so the sequence is on the forward strand." };
    if (!/^GRCz11$/i.test(assembly) && !/danRer11/i.test(assembly)) {
      return { error: "The pasted sequence is labelled " + assembly + ", not GRCz11. Coordinates differ between assemblies." };
    }
    if (hdrEnd && p.seq.length !== hdrEnd - start + 1) {
      return { error: "The header describes " + (hdrEnd - start + 1) + " nt but " + p.seq.length + " nt were pasted, so every coordinate would be off. Copy the response again." };
    }

    var idx = f.pos - start;
    if (idx < 30 || idx > p.seq.length - 30) {
      return { error: "Position " + f.pos + " falls outside the usable part of the pasted window (starts at " + start + ", " + p.seq.length + " nt). Check the coordinate or use a larger flank." };
    }

    var left, right, bracket;
    if (!f.ref) {
      left = p.seq.slice(0, idx + 1);
      right = p.seq.slice(idx + 1);
      bracket = "[-/" + f.alt + "]";
    } else {
      var found = p.seq.slice(idx, idx + f.ref.length);
      if (found !== f.ref) {
        return {
          error: "Reference mismatch: GRCz11 has " + (found || "—") + " at chr" + f.chrom + ":" + f.pos +
            ", but you entered " + f.ref + ". Check the position and strand of your variant call.",
          mismatch: true, found: found
        };
      }
      left = p.seq.slice(0, idx);
      right = p.seq.slice(idx + f.ref.length);
      bracket = "[" + f.ref + "/" + (f.alt || "-") + "]";
    }

    return {
      target: left + bracket + right,
      start: start,
      end: start + p.seq.length - 1,
      assembly: assembly,
      fromHeader: fromHeader
    };
  }

  return {
    rc: rc, clean: clean, gc: gc, pamMatch: pamMatch,
    tmHybrid: tmHybrid, tmDna: tmDna,
    parseTarget: parseTarget,
    buildPeg: buildPeg,
    nickingGuides: nickingGuides,
    tmMatched: tmMatched,
    assemble: assemble,
    design: design,
    defaults: defaults,
    windowBounds: windowBounds,
    ensemblUrl: ensemblUrl,
    parseFasta: parseFasta,
    buildBracketed: buildBracketed,
    SCAFFOLDS: SCAFFOLDS,
    TEVOPREQ1: TEVOPREQ1
  };
});
