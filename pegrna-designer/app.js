/* Danio pegRNA Designer — browser UI. Design logic lives in pegrna.js. */
(function () {
  "use strict";

  var EXAMPLE =
    "GCTGTGATCCTCACAGTCTGGAGAAACACTGTTTCTGCAGGACTGCAGTTTCAGAGGCTCTGGAAGCACTTTGGAAGCTGGAGG" +
    "CCCAGACTGAGCACG[-/CTT]TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGGGCCCAGGCTTAGGCTTCAGTCAGGCT" +
    "CCAGGAGAAGCCTGGTAGACAGTCAGTGGCCTCAAAG";

  var CITATIONS = [
    { text: "Anzalone et al. 2019, Nature 576:149 — prime editing and PE3.", doi: "10.1038/s41586-019-1711-4" },
    { text: "Ponnienselvan et al. 2023, Nucleic Acids Res 51:6966 — PBS–spacer auto-inhibition and Tm-matched PBS.", doi: "10.1093/nar/gkad456" },
    { text: "Nelson et al. 2022, Nat Biotechnol 40:402 — epegRNAs and the tevopreQ1 motif.", doi: "10.1038/s41587-021-01039-7" },
    { text: "Doman et al. 2022, Nat Protoc 17:2431 — pegRNA design and troubleshooting.", doi: "10.1038/s41596-022-00724-4" },
    { text: "Banerjee et al. 2020, Nucleic Acids Res 48:12042 — RNA/DNA hybrid nearest-neighbour parameters used for the PBS Tm.", doi: "10.1093/nar/gkaa572" }
  ];

  var state = { res: null, si: 0, peg: null, tab: "design", coord: null };
  var CLIP = {};
  var $ = function (id) { return document.getElementById(id); };

  function esc(s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c];
    });
  }

  /* ---------- small helpers ---------- */

  function toast(msg) {
    var t = $("toast");
    t.textContent = msg;
    t.classList.add("on");
    clearTimeout(t._t);
    t._t = setTimeout(function () { t.classList.remove("on"); }, 1600);
  }

  function copy(text, what) {
    function fallback() {
      var ta = document.createElement("textarea");
      ta.value = text;
      document.body.appendChild(ta);
      ta.select();
      try { document.execCommand("copy"); toast(what + " copied"); }
      catch (e) { toast("Select the sequence and copy manually"); }
      ta.remove();
    }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(function () { toast(what + " copied"); }, fallback);
    } else fallback();
  }

  function download(text, filename, what) {
    try {
      var blob = new Blob([text], { type: "text/tab-separated-values" });
      var url = URL.createObjectURL(blob);
      var a = document.createElement("a");
      a.href = url; a.download = filename;
      document.body.appendChild(a); a.click();
      setTimeout(function () { URL.revokeObjectURL(url); a.remove(); }, 0);
      toast(what + " downloaded");
    } catch (e) { copy(text, what); }
  }

  function slug(s) {
    return String(s).replace(/[^A-Za-z0-9]+/g, "-").replace(/^-|-$/g, "").slice(0, 40) || "pegrna";
  }

  /* ---------- options ---------- */

  function opts() {
    function n(id, d) { var v = parseInt($(id).value, 10); return isNaN(v) ? d : v; }
    var pbsMin = n("pbsMin", 6), pbsMax = n("pbsMax", 16);
    var rttMin = n("rttMin", 10), rttMax = n("rttMax", 34);
    if (pbsMax < pbsMin) { var a = pbsMin; pbsMin = pbsMax; pbsMax = a; }
    if (rttMax < rttMin) { var b = rttMin; rttMin = rttMax; rttMax = b; }
    var tgt = parseFloat($("tgtTm").value);
    return Object.assign(PEGRNA.defaults(), {
      pam: $("pam").value,
      pbsMin: pbsMin, pbsMax: pbsMax, rttMin: rttMin, rttMax: rttMax,
      minHomology: n("minHom", 7), maxNickToEdit: n("maxN2E", 30),
      targetTm: isNaN(tgt) ? 28.5 : tgt,
      tmModel: $("tmModel").value,
      g5: $("g5").checked, term: $("term").checked,
      useMotif: $("useMotif").checked, motif: $("motif").value,
      scaffold: $("scaffold").value
    });
  }

  /* ---------- run ---------- */

  function run() {
    var o = opts();
    renderLocus();
    var res = PEGRNA.design($("seq").value, o);

    if (res.error) {
      $("candPanel").hidden = true;
      $("benchTag").textContent = "";
      $("benchBody").innerHTML = '<div class="flag bad">' + esc(res.error) + "</div>";
      return;
    }
    if (!res.sites.length) {
      $("candPanel").hidden = true;
      $("benchTag").textContent = "";
      $("benchBody").innerHTML = '<div class="empty">No protospacer places the edit within reach of the nick under these settings.' +
        "<br><br>Raise <b>max nick→edit</b>, lower <b>min 3′ homology</b>, or switch the PAM to NG or NRN for a PAM-flexible editor.</div>";
      return;
    }
    state.res = res; state.si = 0; state.peg = res.sites[0].best; state.tab = "design";
    renderCands();
    selectSite(0);
  }

  function renderCands() {
    var sites = state.res.sites;
    $("candPanel").hidden = false;
    $("candCount").textContent = "— " + sites.length + " found · " + state.res.target.kind;
    $("cands").innerHTML = sites.map(function (s, i) {
      return '<button class="cand" data-i="' + i + '" aria-current="' + (i === state.si) + '">' +
        '<div class="top"><span class="ps">' + esc(s.protospacer) +
        '<span class="pamx">' + esc(s.pamSeq) + "</span></span>" +
        '<span class="score">' + s.best.score + "</span></div>" +
        '<div class="meta"><span class="strandtag">' + s.strand + " strand</span>" +
        "<span>nick→edit " + s.nickToEdit + "</span>" +
        (s.pamChanged ? '<span style="color:var(--good)">PAM disrupted</span>'
          : s.seedChanged ? '<span style="color:var(--good)">seed hit</span>' : "") +
        "<span>" + s.nicking.length + " nicking guides</span></div></button>";
    }).join("");
    Array.prototype.forEach.call($("cands").querySelectorAll(".cand"), function (b) {
      b.onclick = function () { selectSite(+b.dataset.i); };
    });
  }

  function selectSite(i) {
    state.si = i;
    state.peg = state.res.sites[i].best;
    Array.prototype.forEach.call($("cands").querySelectorAll(".cand"), function (b) {
      b.setAttribute("aria-current", String(+b.dataset.i === i));
    });
    renderBench();
  }

  /* ---------- sequence track ---------- */

  function track(s, peg) {
    var view = s.view, sl = 20;
    var pamEnd = s.pamStart + s.pamSeq.length;
    var start = Math.max(0, s.pamStart - sl - 14);
    var end = Math.min(view.length, s.nick + Math.max(peg.rttLen, 24) + 10);
    var html = "";
    for (var i = start; i < end; i++) {
      var cls = [];
      if (i >= s.pamStart - sl && i < s.pamStart) cls.push("t-spacer");
      else if (i >= s.pamStart && i < pamEnd) cls.push("t-pam");
      else cls.push("t-dim");
      if (i >= s.nick - peg.pbsLen && i < s.nick) cls.push("t-pbs");
      html += '<span class="' + cls.join(" ") + '">' + view[i] + "</span>";
      if (i === s.nick - 1) html += '<span style="color:var(--target);font-weight:700">▼</span>';
    }
    var pad = s.nick - start;
    var line2 = '<span class="t-dim">' + new Array(Math.max(0, pad) + 1).join("·") + " </span>";
    var edited = s.downstreamEdited.slice(0, peg.rttLen);
    var e0 = s.nickToEdit, e1 = s.nickToEdit + s.altALen;
    for (var j = 0; j < edited.length; j++) {
      line2 += '<span class="t-rtt' + (j >= e0 && j < e1 ? " t-edit" : "") + '">' + edited[j] + "</span>";
    }
    return '<div class="track"><pre>' + html + "\n" + line2 + "</pre></div>" +
      '<div class="legend">' +
      '<span><b class="k-spacer">■</b> protospacer</span><span><b class="k-pam">■</b> PAM</span>' +
      '<span style="color:var(--target)"><b>▼</b> nick</span>' +
      '<span><b class="k-pbs">▁</b> PBS template</span><span><b class="k-rtt">▁</b> RTT template</span>' +
      '<span><b class="k-edit">■</b> edit</span></div>' +
      '<p class="hint" style="margin:9px 0 0">Shown 5′→3′ on the ' + (s.strand === "+" ? "top" : "bottom") +
      " strand — the strand the pegRNA nicks. The second line is the sequence the reverse transcriptase writes.</p>";
  }

  /* ---------- panels ---------- */

  function designPanel(s, peg, o) {
    var a = PEGRNA.assemble(s, peg, o);
    function seqBlock(label, seq, html) {
      return '<div class="seqrow"><div class="lbl"><h3>' + label + "</h3>" +
        '<span class="len">' + seq.length + ' nt</span><span class="sp"></span>' +
        '<button class="tiny" data-copy="' + esc(seq) + '" data-what="' + esc(label) + '">Copy</button></div>' +
        '<div class="seq">' + (html || esc(seq)) + "</div></div>";
    }
    function st(k, v) { return '<div class="stat"><div class="k">' + k + '</div><div class="v">' + v + "</div></div>"; }

    var match = PEGRNA.tmMatched(s, o);
    var isMatch = match && match.pbsLen === peg.pbsLen;

    var fullHtml = '<span class="s-spacer">' + esc(a.spacer) + "</span>" +
      '<span class="s-scaf">' + esc(a.scaffold) + "</span>" +
      '<span class="s-rtt">' + esc(peg.rtt) + "</span>" +
      '<span class="s-pbs">' + esc(peg.pbs) + "</span>" +
      (a.motif ? '<span class="s-motif">' + esc(a.motif) + "</span>" : "") +
      (a.term ? '<span class="s-term">' + esc(a.term) + "</span>" : "");

    return '<div class="stats">' +
      st("Score", peg.score) + st("PBS", peg.pbsLen + " nt") + st("RTT", peg.rttLen + " nt") +
      st("PBS Tm", isNaN(peg.pbsTm) ? "—" : peg.pbsTm.toFixed(1) + "°C") +
      st("PBS–spacer", peg.autoLen + " bp") +
      st("Nick→edit", s.nickToEdit + " nt") + st("3′ homology", peg.homology3 + " nt") +
      "</div>" +
      (match ? '<div class="flag ' + (isMatch ? "good" : "note") + '" style="margin-bottom:14px">' +
        (isMatch
          ? "This is the Tm-matched PBS for this protospacer"
          : "Tm-matched PBS for this protospacer is <b>" + match.pbsLen + " nt</b> (" +
            match.pbsTm.toFixed(1) + "°C) — you are viewing " + peg.pbsLen + " nt") +
        " at a target of " + o.targetTm + "°C." +
        (isMatch ? "" : ' <button class="tiny" data-usepbs="' + match.pbsLen + '">Switch to ' + match.pbsLen + " nt</button>") +
        "</div>" : "") +
      track(s, peg) +
      '<div style="margin-top:14px">' +
      seqBlock("Spacer" + (a.addG ? " (5′ G added)" : ""), a.spacer) +
      seqBlock("3′ extension (RTT + PBS)", peg.ext,
        '<span class="s-rtt">' + esc(peg.rtt) + '</span><span class="s-pbs">' + esc(peg.pbs) + "</span>") +
      seqBlock("RTT only", peg.rtt) +
      seqBlock("PBS only", peg.pbs) +
      seqBlock("Full pegRNA (DNA sense)", a.full, fullHtml) +
      "</div>" +
      (peg.flags.length ? '<div class="flags">' + peg.flags.map(function (f) {
        return '<div class="flag ' + f.t + '">' + esc(f.m) + "</div>";
      }).join("") + "</div>" : "") +
      '<p class="hint" style="margin:12px 0 0">Sequences are written as DNA. For a synthesised RNA oligo read T as U; the 5′ G matters only for U6 or T7 transcription.</p>';
  }

  function matrixPanel(s, o) {
    var match = PEGRNA.tmMatched(s, o);
    var head = "<tr><th></th>";
    for (var r = o.rttMin; r <= o.rttMax; r++) head += "<th>" + r + "</th>";
    head += "</tr>";

    var rows = "";
    for (var b = o.pbsMin; b <= o.pbsMax; b++) {
      var tmRow = s.pegs.filter(function (x) { return x.pbsLen === b; })[0];
      var hit = match && match.pbsLen === b;
      rows += "<tr><th" + (hit ? ' style="color:var(--accent);font-weight:700"' : "") +
        ' title="' + (tmRow && !isNaN(tmRow.pbsTm) ? "PBS Tm " + tmRow.pbsTm.toFixed(1) + "°C" : "") + '">' +
        b + (hit ? " ◀" : "") + "</th>";
      for (var c = o.rttMin; c <= o.rttMax; c++) {
        var p = null;
        for (var k = 0; k < s.pegs.length; k++) {
          if (s.pegs[k].pbsLen === b && s.pegs[k].rttLen === c) { p = s.pegs[k]; break; }
        }
        if (!p || p.homology3 < o.minHomology) {
          rows += '<td><button class="na" disabled title="not enough 3′ homology">·</button></td>';
          continue;
        }
        var cur = state.peg && state.peg.pbsLen === b && state.peg.rttLen === c;
        var hue = p.score >= 75 ? 158 : p.score >= 55 ? 45 : p.score >= 35 ? 28 : 5;
        rows += '<td><button data-b="' + b + '" data-r="' + c + '" aria-current="' + cur +
          '" style="background:hsl(' + hue + " 55% 45% / .17);color:hsl(" + hue + ' 60% 33%)" ' +
          'title="PBS ' + b + " / RTT " + c + " — score " + p.score + '">' + p.score + "</button></td>";
      }
      rows += "</tr>";
    }

    return '<div class="matrix-wrap"><table class="mx"><thead>' +
      '<tr><th colspan="99" style="text-align:left;padding-bottom:4px">RTT length →</th></tr>' +
      head + "</thead><tbody>" + rows + "</tbody></table></div>" +
      '<p class="hint" style="margin:10px 0 0">Rows are PBS length, columns RTT length. Click a cell to load that design; hover a row label for its PBS Tm. Dots fail the 3′ homology minimum.' +
      (match ? ' <b style="color:var(--accent)">◀</b> marks the Tm-matched PBS (' + match.pbsLen + " nt, " +
        match.pbsTm.toFixed(1) + "°C) for a target of " + o.targetTm + "°C." : "") + "</p>";
  }

  function nickPanel(s) {
    if (!s.nicking.length) {
      return '<div class="empty">No nicking guides found in the sequence provided. Paste another 100 nt on each side and design again.</div>';
    }
    return '<table class="data"><thead><tr><th>Spacer (opposite strand)</th><th>PAM</th><th>Nick offset</th><th>Mode</th><th></th></tr></thead><tbody>' +
      s.nicking.map(function (g) {
        return '<tr><td class="sp">' + esc(g.spacer) + "</td><td>" + esc(g.pam) + "</td>" +
          "<td>" + g.dist + " nt " + g.dir + (g.preferred ? ' <span class="tag pref">40–90</span>' : "") + "</td>" +
          "<td>" + (g.type === "PE3b" ? '<span class="tag pe3b">PE3b</span>' : '<span class="tag">PE3</span>') + "</td>" +
          '<td style="text-align:right"><button class="tiny" data-copy="' + esc(g.spacer) + '" data-what="Nicking spacer">Copy</button></td></tr>';
      }).join("") + "</tbody></table>" +
      '<p class="hint" style="margin:10px 0 0">Offsets are the distance between the two nicks, measured along the pegRNA’s protospacer strand; 40–90 nt is the window to try first. A <b>PE3b</b> spacer matches only after the edit is installed, so the second nick waits for the edited strand — fewer indels, and worth preferring where one exists. Ponnienselvan et al. added the nicking sgRNA to the pegRNA RNP at a 1:10 molar ratio for zebrafish injections.</p>';
  }

  function orderPanel(s, o) {
    var rows = s.pegs.filter(function (p) { return p.homology3 >= o.minHomology; })
      .slice().sort(function (a, b) { return b.score - a.score; });

    var tsv = ["pbs_len\trtt_len\tscore\tpbs_tm_C\tpbs_spacer_bp\tnick_to_edit\tspacer\trtt\tpbs\t3prime_extension\tfull_pegrna"]
      .concat(rows.map(function (p) {
        var a = PEGRNA.assemble(s, p, o);
        return [p.pbsLen, p.rttLen, p.score, isNaN(p.pbsTm) ? "" : p.pbsTm.toFixed(1),
          p.autoLen, s.nickToEdit, a.spacer, p.rtt, p.pbs, p.ext, a.full].join("\t");
      })).join("\n");
    CLIP.tsv = tsv;

    return '<div class="row" style="margin-bottom:14px">' +
      '<button class="tiny" data-dl="tsv" data-what="All ' + rows.length + ' designs">Download all ' + rows.length + " as TSV</button>" +
      '<button class="tiny" data-clip="tsv" data-what="All ' + rows.length + ' designs (TSV)">Copy to clipboard</button></div>' +
      '<table class="data"><thead><tr><th>PBS</th><th>RTT</th><th>Score</th><th>Tm</th><th>3′ extension</th><th></th></tr></thead><tbody>' +
      rows.slice(0, 12).map(function (p) {
        return "<tr><td>" + p.pbsLen + "</td><td>" + p.rttLen + "</td><td>" + p.score + "</td>" +
          "<td>" + (isNaN(p.pbsTm) ? "—" : p.pbsTm.toFixed(1)) + "</td>" +
          '<td class="sp">' + esc(p.ext) + "</td>" +
          '<td style="text-align:right"><button class="tiny" data-copy="' + esc(p.ext) + '" data-what="3′ extension">Copy</button></td></tr>';
      }).join("") + "</tbody></table>" +
      '<p class="hint" style="margin:10px 0 0">Top 12 shown; the TSV holds all ' + rows.length +
      ", one row per PBS × RTT combination, ready for an order sheet.</p>";
  }

  function notesPanel() {
    return '<ul class="note-list">' +
      "<li><b>Match the PBS Tm to the incubation temperature.</b> The PBS is complementary to the 3′ end of its own spacer, so a long PBS folds the pegRNA back on itself and blocks R-loop formation; the pairing runs on into the first three RTT-templated bases wherever they still match the target. Shortening the PBS until its hybrid with the nicked target melts near the growth temperature relieves this. In zebrafish embryos Ponnienselvan et al. used 6 and 7 nt — Tm near 28.5 °C and 37 °C — and both reached roughly 20–26% precise editing with PE2 RNP.</li>" +
      "<li><b>PE3 in embryos.</b> The nicking sgRNA was added to the pegRNA RNP at a 1:10 molar ratio, raising precise editing about 1.2-fold over PE2, to roughly 26–33%. Embryos were injected at the one-cell stage and raised at 28.5 °C.</li>" +
      "<li><b>Screen lengths, do not pick one.</b> There is no universal optimum. A first pass of two or three PBS lengths around the Tm-matched value, against RTT lengths that each extend at least 7 nt past the edit — 20 nt or more for insertions — covers most sites.</li>" +
      "<li><b>Editor choice.</b> PE2 and PEmax have both worked in zebrafish embryos; a systematic comparison of PE2, PE6b, PE6c, PEmax and PE7 for missense edits found PEmax the most efficient for pure prime edits. PEn gives higher total editing but far more indels.</li>" +
      "<li><b>RNP delivery.</b> Prime editor protein plus synthetic pegRNA injected at the one-cell stage works, as does mRNA. Refolding the pegRNA before complexing is worth doing — misfolding driven by spacer–extension complementarity has been reported to cost a large fraction of the activity.</li>" +
      "<li><b>epegRNAs are not automatic.</b> tevopreQ1 is the default in mammalian protocols, but at least one zebrafish study found no benefit from epegRNAs or a La-binding motif. Include a plain pegRNA arm in the first experiment.</li>" +
      "<li><b>Transcription start.</b> Zebrafish U6 promoters want a 5′ G and T7 wants GG. Appending the extra G to make a 21-nt spacer is standard and usually tolerated.</li>" +
      "<li><b>Block re-nicking where you can.</b> A design whose edit lands in the PAM or seed cannot be re-cut once correct. Where the edit does not reach, a silent change inside the RTT often can — check the reading frame before adding one.</li>" +
      "<li><b>Genotyping.</b> Keep primers outside the RTT footprint and sequence the amplicon rather than relying on a restriction site, unless the edit creates one.</li>" +
      "</ul>";
  }

  function renderBench() {
    var s = state.res.sites[state.si], peg = state.peg, o = opts();
    $("benchTag").textContent = "— " + s.strand + " strand · nick→edit " + s.nickToEdit + " nt" + nickLocus(s);

    var tabs = [["design", "Design"], ["matrix", "PBS × RTT"],
      ["nick", "Nicking guides (" + s.nicking.length + ")"], ["order", "Order sheet"], ["notes", "Zebrafish notes"]];

    var body = state.tab === "design" ? designPanel(s, peg, o)
      : state.tab === "matrix" ? matrixPanel(s, o)
      : state.tab === "nick" ? nickPanel(s)
      : state.tab === "order" ? orderPanel(s, o)
      : notesPanel();

    $("benchBody").innerHTML = '<div class="tabs" role="tablist">' +
      tabs.map(function (t) {
        return '<button role="tab" data-tab="' + t[0] + '" aria-selected="' + (state.tab === t[0]) + '">' + t[1] + "</button>";
      }).join("") + "</div>" + body;

    var root = $("benchBody");
    Array.prototype.forEach.call(root.querySelectorAll("[data-tab]"), function (b) {
      b.onclick = function () { state.tab = b.dataset.tab; renderBench(); };
    });
    Array.prototype.forEach.call(root.querySelectorAll("[data-copy]"), function (b) {
      b.onclick = function () { copy(b.dataset.copy, b.dataset.what || "Sequence"); };
    });
    Array.prototype.forEach.call(root.querySelectorAll("[data-clip]"), function (b) {
      b.onclick = function () { copy(CLIP[b.dataset.clip] || "", b.dataset.what || "Table"); };
    });
    Array.prototype.forEach.call(root.querySelectorAll("[data-dl]"), function (b) {
      b.onclick = function () {
        download(CLIP[b.dataset.dl] || "",
          "pegrnas_" + slug(state.coord ? state.coord.label : s.protospacer) + ".tsv",
          b.dataset.what || "Table");
      };
    });
    Array.prototype.forEach.call(root.querySelectorAll("[data-usepbs]"), function (b) {
      b.onclick = function () {
        var want = +b.dataset.usepbs, o2 = opts();
        var cands = s.pegs.filter(function (x) { return x.pbsLen === want && x.homology3 >= o2.minHomology; });
        if (!cands.length) return;
        var same = cands.filter(function (x) { return x.rttLen === peg.rttLen; })[0];
        state.peg = same || cands.slice().sort(function (a, c) { return c.score - a.score; })[0];
        renderBench();
      };
    });
    Array.prototype.forEach.call(root.querySelectorAll("table.mx button[data-b]"), function (b) {
      b.onclick = function () {
        var p = s.pegs.filter(function (x) {
          return x.pbsLen === +b.dataset.b && x.rttLen === +b.dataset.r;
        })[0];
        if (p) { state.peg = p; state.tab = "design"; renderBench(); }
      };
    });
  }

  /* ---------- GRCz11 coordinates ---------- */

  function coordFields() {
    var flank = parseInt($("flank").value, 10);
    if (isNaN(flank)) flank = 150;
    return {
      chrom: String($("chrom").value || "").trim().replace(/^chr/i, ""),
      pos: parseInt(String($("pos").value || "").replace(/[,\s]/g, ""), 10),
      ref: PEGRNA.clean($("refA").value === "-" ? "" : $("refA").value),
      alt: PEGRNA.clean($("altA").value === "-" ? "" : $("altA").value),
      flank: Math.max(60, Math.min(500, flank))
    };
  }

  function renderUrl() {
    var f = coordFields(), box = $("urlbox"), link = $("openUrl");
    if (!f.chrom || isNaN(f.pos) || f.pos < 1) {
      box.innerHTML = "<span>Enter a chromosome and position to build the query.</span>";
      link.removeAttribute("href");
      link.style.opacity = .45;
      return;
    }
    var w = PEGRNA.windowBounds(f), url = PEGRNA.ensemblUrl(f);
    box.innerHTML = "<b>chr" + esc(f.chrom) + ":" + w.start.toLocaleString() + "–" + w.end.toLocaleString() +
      "</b> · " + (w.end - w.start + 1) + " nt · Ensembl REST, GRCz11<br>" + esc(url);
    link.href = url;
    link.style.opacity = 1;
  }

  function buildTarget() {
    var f = coordFields();
    var out = PEGRNA.buildBracketed($("fasta").value, f);
    if (out.error) {
      $("coordMsg").innerHTML = '<div class="flag bad" style="margin-top:12px">' + esc(out.error) + "</div>";
      return false;
    }
    $("seq").value = out.target;
    state.coord = {
      chrom: f.chrom,
      start: out.start,
      label: "chr" + f.chrom + ":" + f.pos.toLocaleString() + " " + (f.ref || "-") + "→" + (f.alt || "-")
    };
    $("coordMsg").innerHTML = '<div class="flag good" style="margin-top:12px">Target built from ' +
      (out.fromHeader ? "the window in the FASTA header" : "the coordinates entered") + ": chr" + esc(f.chrom) + ":" +
      out.start.toLocaleString() + "–" + out.end.toLocaleString() + ". Reference allele verified.</div>";
    setMode("seq");
    run();
    return true;
  }

  function fetchSeq() {
    var f = coordFields(), btn = $("fetchSeq");
    if (!f.chrom || isNaN(f.pos) || f.pos < 1) {
      $("fetchMsg").innerHTML = '<div class="flag bad" style="margin-top:12px">Enter a chromosome and position first.</div>';
      return;
    }
    btn.disabled = true;
    btn.textContent = "Fetching…";
    $("fetchMsg").innerHTML = "";
    fetch(PEGRNA.ensemblUrl(f), { headers: { Accept: "text/x-fasta" } })
      .then(function (r) {
        if (!r.ok) throw new Error("Ensembl returned " + r.status);
        return r.text();
      })
      .then(function (txt) {
        if (!/^>/.test(txt.trim())) throw new Error("the response was not FASTA");
        $("fasta").value = txt;
        $("fetchMsg").innerHTML = '<div class="flag good" style="margin-top:12px">Sequence retrieved from Ensembl.</div>';
        buildTarget();
      })
      .catch(function (err) {
        $("fetchMsg").innerHTML = '<div class="flag bad" style="margin-top:12px">Could not reach Ensembl: ' +
          esc(err.message) + ". Open the link above and paste the response below.</div>";
      })
      .then(function () {
        btn.disabled = false;
        btn.textContent = "Fetch from Ensembl";
      });
  }

  function renderLocus() {
    $("locusTag").innerHTML = state.coord
      ? '<p class="hint" style="margin:8px 0 0">GRCz11 · <b style="font-family:var(--font-mono);color:var(--ink)">' +
        esc(state.coord.label) + "</b></p>"
      : "";
  }

  function nickLocus(s) {
    if (!state.coord) return "";
    var g = state.coord.start + s.genomicNick;
    return " · nick chr" + state.coord.chrom + ":" + (g - 1).toLocaleString() + "/" + g.toLocaleString();
  }

  function setMode(m) {
    $("pane-seq").hidden = m !== "seq";
    $("pane-coord").hidden = m !== "coord";
    Array.prototype.forEach.call($("inputTabs").querySelectorAll("[data-mode]"), function (b) {
      b.setAttribute("aria-selected", String(b.dataset.mode === m));
    });
    if (m === "coord") renderUrl();
  }

  /* ---------- wiring ---------- */

  $("citations").innerHTML = "<b>Design rules from:</b> " + CITATIONS.map(function (c) {
    return esc(c.text) + ' <a href="https://doi.org/' + c.doi + '" target="_blank" rel="noopener">doi:' + c.doi + "</a>";
  }).join(" ") + " Genome sequence from " +
    '<a href="https://rest.ensembl.org" target="_blank" rel="noopener">Ensembl REST</a> (GRCz11).';

  $("run").onclick = run;
  $("loadEx").onclick = function () { $("seq").value = EXAMPLE; state.coord = null; setMode("seq"); run(); };
  $("buildTarget").onclick = buildTarget;
  $("fetchSeq").onclick = fetchSeq;
  $("copyUrl").onclick = function () {
    var f = coordFields();
    if (f.chrom && !isNaN(f.pos)) copy(PEGRNA.ensemblUrl(f), "Ensembl URL");
  };

  Array.prototype.forEach.call($("inputTabs").querySelectorAll("[data-mode]"), function (b) {
    b.onclick = function () { setMode(b.dataset.mode); };
  });
  ["chrom", "pos", "refA", "altA", "flank"].forEach(function (id) {
    $(id).addEventListener("input", renderUrl);
  });
  ["pbsMin", "pbsMax", "rttMin", "rttMax", "minHom", "maxN2E", "tgtTm", "tmModel", "pam",
    "scaffold", "motif", "useMotif", "g5", "term"].forEach(function (id) {
      $(id).addEventListener("change", function () { if (state.res) run(); });
    });
  Array.prototype.forEach.call(document.querySelectorAll("[data-tm]"), function (b) {
    b.onclick = function () { $("tgtTm").value = b.dataset.tm; if (state.res) run(); };
  });
  $("seq").addEventListener("input", function () {
    if (state.coord) { state.coord = null; renderLocus(); }
  });
  $("seq").addEventListener("keydown", function (e) {
    if ((e.metaKey || e.ctrlKey) && e.key === "Enter") run();
  });

  (function theme() {
    var root = document.documentElement, saved = null;
    try { saved = localStorage.getItem("peg-theme"); } catch (e) {}
    if (saved) root.setAttribute("data-theme", saved);
    $("themeBtn").onclick = function () {
      var cur = root.getAttribute("data-theme");
      var next = cur === "dark" ? "light" : cur === "light" ? "dark"
        : (window.matchMedia("(prefers-color-scheme: dark)").matches ? "light" : "dark");
      root.setAttribute("data-theme", next);
      try { localStorage.setItem("peg-theme", next); } catch (e) {}
    };
  })();

  renderUrl();
  $("seq").value = EXAMPLE;
  run();
})();
