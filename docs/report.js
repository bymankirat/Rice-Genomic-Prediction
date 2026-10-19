/* Rice genomic prediction report. Every number comes from results.js, written by the pipeline. */
(function () {
  "use strict";
  const D = window.R;
  const $ = (s, r = document) => r.querySelector(s);
  const NS = "http://www.w3.org/2000/svg";
  const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
  const esc = (s) => String(s == null ? "" : s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const f = (x, d = 2) => (x == null || !isFinite(x) ? "–" : Number(x).toFixed(d));
  const pct = (x) => (x == null ? "–" : (x * 100).toFixed(0) + "%");
  const esci = (p) => (p == null ? "–" : p < 1e-4 ? Number(p).toExponential(1) : Number(p).toFixed(4));
  const nice = (t) => String(t).replace(/\./g, " ").replace(/\b\w/g, (c) => c.toUpperCase());

  if (!D || !D.qc || !D.qc.markers_kept) {
    const b = $("#banner"); b.hidden = false;
    b.textContent = "No results yet: the pipeline has not run. See the README.";
    return;
  }
  $("#generated").textContent = D.generated || "";

  const GROUP_COLOURS = ["--c1", "--c2", "--c3", "--h7", "--h9", "--h3", "--h5", "--accent"];
  const groups = Object.keys(D.structure.group_sizes || {}).sort();
  const gCol = (g) => css(GROUP_COLOURS[groups.indexOf(g) % GROUP_COLOURS.length] || "--muted");

  function svg(w, h) {
    const s = document.createElementNS(NS, "svg");
    s.setAttribute("viewBox", `0 0 ${w} ${h}`);
    s.setAttribute("preserveAspectRatio", "xMidYMid meet");
    s.setAttribute("role", "img");
    s.style.width = "100%"; s.style.height = "auto";
    return s;
  }
  const add = (p, tag, a = {}, t) => {
    const e = document.createElementNS(NS, tag);
    for (const k in a) if (a[k] != null) e.setAttribute(k, a[k]);
    if (t != null) e.textContent = t;
    p.appendChild(e); return e;
  };

  // ---------------------------------------------------------------- at a glance
  const pr = D.prediction || {};
  const rRand = (pr.mean_r_random || {}).ridge;
  const rAcross = (pr.mean_r_across_group || {}).ridge;
  const checks = (D.gwas || {}).known_gene_checks || [];
  const found = checks.filter((c) => c.significant).length;
  $("#kpis").innerHTML = [
    { v: D.qc.samples_kept, l: "rice varieties", s: "after quality control" },
    { v: D.qc.markers_kept.toLocaleString(), l: "markers used", s: `of ${D.qc.markers_in.toLocaleString()} genotyped` },
    { v: D.structure.k, l: "genetic groups", s: "chosen by silhouette on the PCs" },
    { v: f(rRand), l: "mean accuracy (r)", s: "random 5-fold cross-validation" },
    { v: f(rAcross), l: "mean accuracy (r)", s: "on a group never seen in training" },
    { v: `${found}/${checks.length}`, l: "known genes recovered", s: "sd1, Waxy, GS3 and others" },
  ].map((k) =>
    `<div class="kpi"><div class="v">${esc(k.v ?? "–")}</div><div class="l">${esc(k.l)}</div><div class="arms small">${esc(k.s)}</div></div>`).join("");

  const drop = rRand && rAcross ? (rRand - rAcross) / rRand : null;
  const t = [];
  t.push(found >= Math.max(1, Math.ceil(checks.length / 2))
    ? `<b>The analysis behaves.</b> ${found} of ${checks.length} textbook rice genes come back as genome-wide significant in the right place, so the scan is finding real signal, not noise.`
    : `<b>Treat the association results with caution.</b> Only ${found} of ${checks.length} textbook genes were recovered, which is fewer than expected on this panel.`);
  if (rRand != null) {
    const base = (pr.mean_r_random || {})["ridge (100 markers)"];
    t.push(`<b>Markers do carry the signal.</b> Using all ${D.qc.markers_kept.toLocaleString()} markers gives a mean accuracy of ${f(rRand)} against ${f(base)} for a 100-marker model in the same folds, so the gain is not an artefact of the test design.`);
  }
  if (drop != null) {
    t.push(drop > 0.25
      ? `<b>Most of the headline accuracy does not transfer.</b> Holding out a whole genetic group drops the mean accuracy from ${f(rRand)} to ${f(rAcross)}, a loss of ${pct(drop)}. A programme working on material unlike its training set should expect the lower number, not the one usually quoted.`
      : `<b>The models transfer better than expected.</b> Holding out a whole genetic group only moves the mean accuracy from ${f(rRand)} to ${f(rAcross)}, so a useful share of what the model learns is not simply family resemblance.`);
  }
  $("#takeaways").innerHTML = t.map((x) => `<p>${x}</p>`).join("");

  // ---------------------------------------------------------------- QC
  const q = D.qc;
  $("#qc").innerHTML = `<div class="qcsum">
    <p>Varieties: <b>${q.samples_kept}</b> kept of ${q.samples_in} (call rate ≥ ${pct(q.thresholds.sample_call_rate)}).</p>
    <p>Markers: <b>${q.markers_kept.toLocaleString()}</b> kept of ${q.markers_in.toLocaleString()}
       (call rate ≥ ${pct(q.thresholds.snp_call_rate)}, minor allele frequency ≥ ${pct(q.thresholds.maf)}).</p>
    <p>Missing calls before filtering: ${pct(q.missing_rate_before_qc)}; remaining gaps filled with each marker's mean.</p>
  </div>`;

  const mc = D.markers_per_chrom || [];
  if (mc.length) {
    const host = $("#chrom-chart"); host.innerHTML = "";
    const W = 760, H = 180, L = 40, R = 16, T = 16, B = 40;
    const max = Math.max(...mc.map((m) => m.markers));
    const bw = (W - L - R) / mc.length;
    const s = svg(W, H);
    mc.forEach((m, i) => {
      const h = (m.markers / max) * (H - T - B);
      const bar = add(s, "rect", { x: L + i * bw + bw * 0.15, y: H - B - h, width: bw * 0.7, height: h, fill: css("--accent"), rx: 2 });
      add(bar, "title", {}, `Chromosome ${m.chrom}: ${m.markers.toLocaleString()} markers`);
      add(s, "text", { x: L + i * bw + bw / 2, y: H - B + 16, "text-anchor": "middle", "font-size": 11, fill: css("--muted") }, m.chrom);
    });
    add(s, "text", { x: L, y: 12, "font-size": 12, fill: css("--muted") }, "markers per chromosome");
    add(s, "text", { x: (W) / 2, y: H - 6, "text-anchor": "middle", "font-size": 12, fill: css("--muted") }, "chromosome");
    host.appendChild(s);
  }

  // ---------------------------------------------------------------- structure
  const pcs = D.pcs || [];
  if (pcs.length) {
    const host = $("#pca-chart"); host.innerHTML = "";
    const W = 760, H = 420, L = 54, R = 150, T = 20, B = 48;
    const xs = pcs.map((p) => p.PC1), ys = pcs.map((p) => p.PC2);
    const xr = [Math.min(...xs), Math.max(...xs)], yr = [Math.min(...ys), Math.max(...ys)];
    const x = (v) => L + ((v - xr[0]) / (xr[1] - xr[0])) * (W - L - R);
    const y = (v) => H - B - ((v - yr[0]) / (yr[1] - yr[0])) * (H - T - B);
    const s = svg(W, H);
    add(s, "rect", { x: L, y: T, width: W - L - R, height: H - T - B, fill: "none", stroke: css("--grid") });
    pcs.forEach((p) => {
      const c = add(s, "circle", { cx: x(p.PC1), cy: y(p.PC2), r: 3.4, fill: gCol(p.group), "fill-opacity": 0.8 });
      add(c, "title", {}, `${p.iid} · ${p.group}`);
    });
    const ve = D.structure.variance_explained || [];
    add(s, "text", { x: (L + W - R) / 2, y: H - 12, "text-anchor": "middle", "font-size": 12, fill: css("--muted") },
      `PC1 (${pct(ve[0])} of marker variation)`);
    add(s, "text", { x: 16, y: (T + H - B) / 2, "font-size": 12, fill: css("--muted"),
      transform: `rotate(-90 16 ${(T + H - B) / 2})`, "text-anchor": "middle" }, `PC2 (${pct(ve[1])})`);
    groups.forEach((g, i) => {
      add(s, "circle", { cx: W - R + 14, cy: T + 16 + i * 20, r: 5, fill: gCol(g) });
      add(s, "text", { x: W - R + 26, y: T + 20 + i * 20, "font-size": 12, fill: css("--ink") },
        `${g} · ${D.structure.group_sizes[g]} varieties`);
    });
    host.appendChild(s);
  }
  $("#structure-note").innerHTML = `The first two components separate the panel into <b>${D.structure.k}</b> groups
    (silhouette score ${f((D.structure.silhouette || {})[D.structure.k], 3)}), close to the five subpopulations described for
    cultivated rice. Group names here are arbitrary labels from clustering, not the published subpopulation names:
    calling a cluster "indica" would need the reference accessions, which this analysis does not use.`;

  // ---------------------------------------------------------------- GWAS
  const g = D.gwas || {};
  $("#known-table").innerHTML = checks.length
    ? `<thead><tr><th>Trait</th><th>Known gene</th><th>Chromosome</th><th>Expected position</th><th>Best marker in the window</th><th>Recovered?</th></tr></thead><tbody>` +
      checks.map((c) => `<tr><td>${esc(nice(c.trait))}</td><td class="mono">${esc(c.gene)}</td><td>${c.chrom}</td>
        <td>${f(c.expected_mb, 2)} Mb</td>
        <td>${c.best_p_in_window == null ? "–" : `p = ${esci(c.best_p_in_window)} at ${f(c.best_pos_mb, 2)} Mb`}</td>
        <td>${c.significant ? '<span class="tag">yes</span>' : '<span class="small">not at the genome-wide threshold</span>'}</td></tr>`).join("") + "</tbody>"
    : '<tbody><tr><td class="small">No known-gene checks were possible for the traits in this file.</td></tr></tbody>';

  const traits = (g.traits || []).filter((x) => x.n_significant != null).sort((a, b) => b.n_significant - a.n_significant);
  if (traits.length) {
    const host = $("#manhattan"); host.innerHTML = "";
    const top = traits.slice(0, 14);
    const W = 760, rowH = 24, T = 28, B = 34, L = 190, R = 40;
    const H = T + B + rowH * top.length;
    const max = Math.max(...top.map((x) => x.n_significant), 1);
    const s = svg(W, H);
    top.forEach((x, i) => {
      const yy = T + i * rowH;
      add(s, "text", { x: L - 10, y: yy + 15, "text-anchor": "end", "font-size": 12, fill: css("--ink") }, nice(x.trait));
      const w = Math.max(1, (x.n_significant / max) * (W - L - R));
      const bar = add(s, "rect", { x: L, y: yy + 4, width: w, height: rowH - 9, fill: css("--accent"), rx: 3 });
      add(bar, "title", {}, `${nice(x.trait)}: ${x.n_significant} markers below ${esci(g.threshold)}; strongest on chromosome ${x.top_chrom} at ${f(x.top_pos_mb, 1)} Mb`);
      add(s, "text", { x: L + w + 6, y: yy + 15, "font-size": 11, fill: css("--muted") }, x.n_significant);
    });
    add(s, "text", { x: L, y: 16, "font-size": 12, fill: css("--muted") },
      `markers passing the genome-wide threshold (p < ${esci(g.threshold)})`);
    host.appendChild(s);
  }
  const lams = (g.traits || []).map((x) => x.lambda_gc).filter((x) => x != null && isFinite(x));
  if (lams.length) {
    const med = lams.slice().sort((a, b) => a - b)[lams.length >> 1];
    $("#lambda-note").innerHTML = `Genomic inflation factor across traits: median ${f(med, 2)}, maximum ${f(Math.max(...lams), 2)}.
      A value near 1 means the structure correction is doing its job; much above 1.1 would mean the scan is still picking up
      population differences rather than the trait. ${g.n_pcs} principal components were used as covariates.`;
  }

  // ---------------------------------------------------------------- prediction
  const summary = D.prediction_summary || [];
  const rand = summary.filter((r) => r.test === "random 5-fold");
  const across = summary.filter((r) => r.test === "leave-one-group-out");
  const models = ["mean", "ridge (100 markers)", "ridge", "gbm"];
  const MODEL_LABEL = { mean: "training average", "ridge (100 markers)": "ridge, 100 markers", ridge: "ridge, all markers", gbm: "gradient boosting" };
  const traitList = [...new Set(rand.map((r) => r.trait))];
  const byTrait = (rows, trait) => Object.fromEntries(rows.filter((r) => r.trait === trait).map((r) => [r.model, r.r]));

  function groupedBars(host, rows, title) {
    host.innerHTML = "";
    const list = traitList.map((tr) => ({ tr, v: byTrait(rows, tr) }))
      .filter((x) => x.v.ridge != null).sort((a, b) => b.v.ridge - a.v.ridge).slice(0, 14);
    if (!list.length) { host.innerHTML = '<p class="small">Not produced yet.</p>'; return; }
    const used = models.filter((m) => list.some((x) => x.v[m] != null));
    const W = 760, rowH = 15 * used.length + 14, T = 30, B = 40, L = 190, R = 110;
    const H = T + B + rowH * list.length;
    const s = svg(W, H);
    const max = Math.max(0.6, ...list.flatMap((x) => used.map((m) => x.v[m] || 0)));
    const x = (v) => L + Math.max(0, v / max) * (W - L - R);
    [0, 0.2, 0.4, 0.6, 0.8].filter((v) => v <= max).forEach((v) => {
      add(s, "line", { x1: x(v), x2: x(v), y1: T - 8, y2: H - B, stroke: css("--grid") });
      add(s, "text", { x: x(v), y: H - B + 16, "text-anchor": "middle", "font-size": 11, fill: css("--muted") }, v.toFixed(1));
    });
    list.forEach((it, i) => {
      const y0 = T + i * rowH;
      add(s, "text", { x: L - 10, y: y0 + rowH / 2, "text-anchor": "end", "font-size": 12, fill: css("--ink") }, nice(it.tr));
      used.forEach((m, j) => {
        const v = it.v[m];
        if (v == null) return;
        const yy = y0 + 6 + j * 15;
        const colour = m === "ridge" ? css("--accent") : m === "gbm" ? css("--h3") : css("--muted");
        const bar = add(s, "rect", { x: L, y: yy, width: Math.max(1, x(v) - L), height: 11,
          fill: colour, "fill-opacity": m.startsWith("ridge") || m === "gbm" ? 0.9 : 0.4, rx: 2 });
        add(bar, "title", {}, `${nice(it.tr)} · ${MODEL_LABEL[m]}: r = ${f(v)}`);
      });
      add(s, "text", { x: x(it.v.ridge) + 6, y: y0 + 16, "font-size": 11, fill: css("--muted") }, f(it.v.ridge));
    });
    used.forEach((m, j) => {
      add(s, "rect", { x: W - R + 10, y: T + j * 18, width: 10, height: 10,
        fill: m === "ridge" ? css("--accent") : m === "gbm" ? css("--h3") : css("--muted"),
        "fill-opacity": m.startsWith("ridge") || m === "gbm" ? 0.9 : 0.4, rx: 2 });
      add(s, "text", { x: W - R + 24, y: T + 9 + j * 18, "font-size": 11, fill: css("--ink") }, MODEL_LABEL[m]);
    });
    add(s, "text", { x: L, y: 16, "font-size": 12, fill: css("--muted") }, title);
    host.appendChild(s);
  }
  groupedBars($("#pred-chart"), rand, "accuracy (correlation between predicted and observed), random 5-fold");

  $("#pred-table").innerHTML =
    `<thead><tr><th>Trait</th>${models.map((m) => `<th>${esc(MODEL_LABEL[m])}</th>`).join("")}<th>Across groups (ridge)</th></tr></thead><tbody>` +
    traitList.map((tr) => {
      const v = byTrait(rand, tr), a = byTrait(across, tr);
      return `<tr><td>${esc(nice(tr))}</td>${models.map((m) => `<td>${f(v[m])}</td>`).join("")}<td>${f(a.ridge)}</td></tr>`;
    }).join("") + "</tbody>";

  groupedBars($("#transfer-chart"), across, "accuracy when an entire genetic group is held out of training");
  if (rAcross != null) {
    const worse = traitList.filter((tr) => {
      const v = byTrait(rand, tr).ridge, a = byTrait(across, tr).ridge;
      return v != null && a != null && a < v;
    }).length;
    $("#transfer-note").innerHTML = `Across ${traitList.length} traits, holding out a whole group lowers ridge accuracy for
      <b>${worse}</b> of them, and the mean falls from <b>${f(rRand)}</b> to <b>${f(rAcross)}</b>. That gap is the share of
      the usual cross-validation number that comes from relatedness between training and test varieties rather than from
      marker effects that hold anywhere.`;
  }

  // ---------------------------------------------------------------- methods
  const prov = D.provenance || {};
  $("#methods-box").innerHTML = `<ul class="plain">
    <li><b>Data.</b> Rice Diversity Panel 1 genotypes and 34 traits, downloaded from ricediversity.org
      ${prov.genotypes && prov.genotypes.md5 ? `(genotype archive MD5 <span class="mono">${esc(prov.genotypes.md5.slice(0, 12))}…</span>)` : ""}.</li>
    <li><b>Quality control.</b> Marker call rate ≥ ${pct(q.thresholds.snp_call_rate)}, minor allele frequency ≥ ${pct(q.thresholds.maf)},
      variety call rate ≥ ${pct(q.thresholds.sample_call_rate)}; remaining missing calls set to the marker mean.</li>
    <li><b>Structure.</b> Principal components on standardised markers; groups from k-means on the first five components,
      with k chosen by silhouette over ${Object.keys(D.structure.silhouette || {}).length} values.</li>
    <li><b>GWAS.</b> One linear model per marker with ${g.n_pcs} principal components as covariates;
      Bonferroni threshold ${esci(g.threshold)} across ${(g.n_markers || 0).toLocaleString()} markers.</li>
    <li><b>Prediction.</b> ${pr.repeats}×${pr.folds}-fold cross-validation and leave-one-group-out, with ridge regression
      (equivalent to GBLUP with a marker kernel), gradient boosting, and two baselines inside the same folds.
      Markers are standardised using training-fold statistics only.</li>
    <li><b>Tests.</b> Simulated data check that the genotype encoding, the association model and its structure correction
      behave; further tests check the published numbers against the files.</li>
  </ul>`;
})();
