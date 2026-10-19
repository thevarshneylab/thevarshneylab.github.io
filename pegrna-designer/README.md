# Danio pegRNA Designer

Prime editing guide RNA design for zebrafish — pegRNA spacers, 3′ extensions, and PE3/PE3b nicking guides, with the PBS length chosen by melting temperature rather than by convention.

**→ [Use it here](https://thevarshneylab.github.io/pegrna-designer/)**

It runs entirely in your browser. Nothing is uploaded, apart from the Ensembl request you trigger yourself.

---

## Why this exists

Most pegRNA tools hand you a 13 nt PBS because that is the number in the original prime editing paper. For zebrafish that is usually the wrong number.

The PBS is complementary to the 3′ end of its own spacer, so a long PBS folds the pegRNA back on itself and competes with R-loop formation at the target — an auto-inhibitory duplex that the pegRNA carries by construction. The pairing also runs on into the first three RTT-templated bases wherever those still match the target. Ponnienselvan et al. showed that shortening the PBS until its hybrid with the nicked target melts near the growth temperature relieves this, and that in zebrafish embryos the right answer is **6 or 7 nt**, not 13.

A design tool whose PBS range starts at 10 cannot express that answer. This one starts at 6, computes the PBS:target melting temperature with RNA/DNA nearest neighbours rather than a DNA approximation, and marks the Tm-matched length for 28.5 °C.

## What it does

- **Finds every protospacer** that places the edit 3′ of the nick, on both strands, for SpCas9 (NGG) or the PAM-flexible SpG and SpRY variants.
- **Sweeps PBS × RTT** as a scored matrix — every cell is a real design, clickable, and exportable as a TSV with one row per combination.
- **Marks the Tm-matched PBS** for your incubation temperature, with one click to switch to it.
- **Reports the PBS–spacer duplex** in base pairs, including the run-on into unedited RTT bases.
- **Flags re-nicking**: designs whose edit lands in the PAM or seed cannot be re-cut once correct.
- **Finds nicking guides** for PE3, and marks PE3b spacers — the ones that match only after the edit is installed.
- **Looks up GRCz11 coordinates** through the Ensembl REST API, verifying that the reference allele you entered is actually the base at that position before designing anything.

## Worked example

Open the tool and press **Load example**. It loads a `+1 CTT` insertion at the human *HEK3* locus — the canonical prime editing test case — and returns:

```
GGCCCAGACTGAGCACGTGA TGG      + strand, nick→edit 0 nt, PAM disrupted
```

At PBS 13 / RTT 13 the 3′ extension is:

```
TCTGCCATCAAAG CGTGCTCAGTCTG
     RTT            PBS
```

which is the extension published for this site. At a 28.5 °C target the tool instead marks a shorter PBS, because the site is GC-rich: 7 nt melts at 28.5 °C and 8 nt at 37 °C.

The same arithmetic is asserted in `test.js`, so a change that breaks the nick position or the reverse complement fails CI rather than quietly shipping wrong oligos.

## Input format

Mark the edit in the sequence with brackets:

| Edit | Notation |
|---|---|
| Substitution | `…ACGT[C/T]ACGT…` |
| Insertion | `…ACGT[-/CTT]ACGT…` |
| Deletion | `…ACGT[AGG/-]ACGT…` |

Give at least 100 nt on each side so nicking guides can be found.

Or switch to the **GRCz11 coordinates** tab, enter `chr:position ref>alt`, and press **Fetch from Ensembl**. The tool reads the window coordinates out of the FASTA header rather than trusting the form, rejects a window returned on the reverse strand or under a different assembly, checks the pasted length against the header's span, and refuses to design if GRCz11 does not carry your reference allele at that position — the usual way a pegRNA order goes wrong.

## Files

| File | |
|---|---|
| `index.html` | markup and styles |
| `pegrna.js` | design logic; runs in the browser and in Node |
| `app.js` | browser UI |
| `test.js` | test suite — `node test.js` |

No build step and no dependencies.

## Scope

The score is a triage heuristic built from published design rules — nick-to-edit distance, PBS Tm, 3′ homology, PAM and seed disruption, poly-T runs, a 3′ extension starting with C. It is not a trained efficiency prediction, and the papers are consistent that locus behaviour varies: test several designs in parallel. For a trained model run a shortlist through PRIDICT or DeepPrime, neither of which was trained on zebrafish embryos.

Off-target search is out of scope. Copy a spacer into [CRISPOR](https://crispor.gi.ucsc.edu/) or [Cas-OFFinder](http://www.rgenome.net/cas-offinder/) against *Danio rerio* GRCz11/danRer11, and check the nicking sgRNA as well as the pegRNA spacer.

## Citing

See `CITATION.cff`. The design rules come from:

- Anzalone et al. 2019, *Nature* 576:149 — prime editing and PE3. [doi:10.1038/s41586-019-1711-4](https://doi.org/10.1038/s41586-019-1711-4)
- Ponnienselvan et al. 2023, *Nucleic Acids Res* 51:6966 — PBS–spacer auto-inhibition, Tm-matched PBS, zebrafish PE2/PE3 RNP. [doi:10.1093/nar/gkad456](https://doi.org/10.1093/nar/gkad456)
- Nelson et al. 2022, *Nat Biotechnol* 40:402 — epegRNAs and tevopreQ1. [doi:10.1038/s41587-021-01039-7](https://doi.org/10.1038/s41587-021-01039-7)
- Doman et al. 2022, *Nat Protoc* 17:2431 — pegRNA design and troubleshooting. [doi:10.1038/s41596-022-00724-4](https://doi.org/10.1038/s41596-022-00724-4)

PBS:target hybrid melting temperatures use RNA/DNA nearest-neighbour parameters from Banerjee et al. 2020, *Nucleic Acids Res* 48:12042 ([doi:10.1093/nar/gkaa572](https://doi.org/10.1093/nar/gkaa572)); the DNA/DNA option uses SantaLucia 1998. Ponnienselvan et al. computed theirs in MELTING 5 with the Sugimoto 1995 set, so absolute values here may differ by a degree or two — the ranking of lengths, which is what picks the pegRNA, does not.

## Licence

MIT — see `LICENSE`.

This tool is developed using Claude Opus 5.
