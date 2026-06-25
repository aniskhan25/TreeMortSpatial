# Revision Plan — *Morphological Typology of Dead-Tree Clusters in Finland's Boreal Forests*

Derived from a four-reviewer panel (spatial statistics, disturbance ecology, RS/ML
validation, landscape connectivity). Three reviewers: Major Revision; one: Reject.
This plan addresses every Tier-1/2/3 issue and specifies the pipeline re-runs needed.

---

## FINAL PANEL VERDICT (2026-06-24): unanimous Accept / Minor Revision

After re-consulting all four reviewers on the completed revision:

| Reviewer | Original | Final | Outstanding |
|---|---|---|---|
| Vasquez (spatial stats) | Major | **Accept (Minor)** | 5 presentational conditions |
| Lindqvist (ecology) | **Reject** | **Accept (Minor)** | 4 presentational |
| Anand (landscape) | Major | **Accept (Minor)** | 4 presentational |
| Tanaka (RS/ML) | Major | **Minor Revision** | **1 computational gate (Moran's I)** |

All four credited the authors for running the falsification gates and letting the negative
results (eps-instability, AR-cap artefact, anisotropy null, AUC-0.727 confound) rewrite the
paper. The discrete typology is correctly retired; the gradient + negative-result spine stands.

**Convergent conditions (≥3 reviewers):**
- **Lead with the negatives in the ABSTRACT**: state plainly that no discrete types were
  recoverable (eps-unstable, BIC favours k>2) and that elongation carries no anisotropic/
  directional signal — these are the paper's most transferable findings. (Vasquez #1, Lindqvist #2, Anand #1)
- Promote the gate failures + anisotropy refutation into main-text Results with figures.
- Purge any residual "two types / Type-0 vs Type-1 as a real partition" language.
- Report prevalence as cap- and eps-conditional (retire standalone 23.5/76.5%).
- Reconcile the stale NOTEBOOK.md provenance (Vasquez #5 — already done this session).

**Tanaka's one non-negotiable computational gate (blocks her full Accept):**
Propagate detection error into the **Moran's I** statistics (not just type assignments): re-run
the ~100-realization perturbation, recompute multi-scale Moran's I each time, report the
distribution and the fraction preserving sign + Bonferroni significance. PASS if preserved in
≥95% of realizations at the 200–500 m scales; else caveat/withdraw the autocorrelation claim.
Smaller: report connectivity-reach effect size (not just p); state the ρ=0.07 (gradient) vs
AUC=0.727 (discrete) contrast explicitly as the reason the gradient survives while types don't.

## Tanaka's Moran's-I propagation gate — RUN (2026-06-24), claim narrowed

`analysis/gate_moran_propagation.py` (+ `analysis/out/gate_moran_propagation.json`): 100
realizations, 15% random false-negative detection perturbation, multi-scale Moran's I recomputed
each time on the same 2,000-cluster stratified subsample; analytical p (the permutations=199
p_sim floor of 0.005 cannot reach Bonferroni α*=0.0025 — a latent flaw in the original).

Result — the gate did NOT cleanly pass, and that is informative:
- Magnitudes are weak everywhere (I ≤ 0.11).
- SIGN of I is robust to detection perturbation for clark_evans / aspect_ratio / compactness
  (≥96% at 200–500 m) but NOT for fractal_dim (60–72% — near-zero, noise-indistinguishable).
- SIGNIFICANCE: under analytical p, ONLY Clark–Evans clears Bonferroni, and only at 500 m–5 km
  (sig 72–100%); fractal_dim/aspect_ratio/compactness do not. The 200 m cells fail because most
  clusters have no neighbour at 200 m (islands → inflated variance).

Action taken (per Tanaka's "Fail → caveat or withdraw"): the §4.5 "significant across most
combinations" claim was withdrawn and narrowed to "within-cluster dispersion (Clark–Evans) shows
weak but robust positive autocorrelation at 0.5–5 km; shape metrics show at most marginal,
detection-sensitive structure." Methods now specify analytical p + the detection-robustness gate;
the heatmap was regenerated with analytical p (only Clark–Evans marked Bonferroni-significant) and
its caption updated. This satisfies Tanaka's outstanding condition → her Minor becomes clear.

## Presentational pass (2026-06-24)
Abstract recast to LEAD with the negative results (gradient real but not discrete; the four gate
failures enumerated incl. the anisotropy null and AUC=0.727 confound) and to close on the
transferable-caution framing all three other reviewers requested. (Abstract, §4.5, methods, and
the Moran figure now mutually consistent.)

## DECISION (2026-06-23): PATH A — data-availability check resolved it

Checked public availability of the Path B datasets for the SE-Finland sites:
- **MS-NFI species/volume raster (Luke)** — PUBLIC (CC-BY 4.0, GeoTIFF, 16 m, incl. spruce volume;
  Paituli / Luke Inspire download). Enables host-conditioning.
- **Metsäkeskus forest-damage data** — PUBLIC (Forest Information Act, avoin metsätieto / Metsään.fi;
  downloads + APIs).
- **Senf & Seidl / European Forest Disturbance Atlas attribution** — PUBLIC (Zenodo, GeoTIFF, 30 m)
  BUT it **fuses wind and bark beetle into ONE class** ("grouped for technical reasons") — exactly
  the two endpoints Path B would need to separate. So the public agent map cannot validate the
  wind-vs-beetle distinction.
- **Detection-model ground truth (Tanaka R1/R2)** — NOT public; requires field plots or manual
  photo-interpretation the team would have to produce.

**Verdict:** the data that would let Path B *recover a disturbance-agent typology* (wind/beetle
separation + detection ground truth) is **not publicly available**; only host-conditioning data is.
Per the user's rule ("if not public, do Path A"), we finalize as **Path A** — the honest
morphological-gradient methods/cautionary paper. (A partial host-conditioned point-process upgrade
remains possible later with MS-NFI alone, but it cannot resurrect the agent labels.)

**Path A status: essentially complete.** All reframing, the four gates, title/scope correction,
and minor fixes are done; the paper now claims a regional morphological gradient and explicitly
refutes the discrete agent typology. Section headers and keywords updated off "typology."

## 0. The strategic decision (superseded by the decision box above)

The paper currently makes three claims it cannot support with the data in hand:
(1) the typology is **independently validated**, (2) coverage is **national**, and
(3) the types map to **disturbance agents** ("Wind-thrown" / "Beetle/Patch").

Two coherent paths forward — pick one before editing prose:

- **Path A — Honest morphological-methods paper (recommended; achievable with current data).**
  Reframe as: *"A reproducible, sensitivity-tested pipeline for characterizing the
  spatial morphology of remotely-detected mortality clusters, revealing a dominant
  compactness–elongation gradient."* Drop agent labels, drop "independent validation,"
  correct the spatial-extent claim. This is publishable and defensible.

- **Path B — Supervised attribution paper (more work, higher ceiling).**
  Acquire the validation layers the reviewers demand (MS-NFI spruce-dominance layer,
  Metsäkeskus / Senf–Seidl disturbance-agent maps, multi-year imagery) and convert the
  unsupervised typology into a *validated* agent classification. Keeps the ecological
  story but earns it.

The plan below executes **Path A** as the backbone and marks **[Path B]** extensions.

---

## TIER 1 — Critical (must fix; these are the spine)

### C1. Retire the "independent validation" framing of the held-out RF (AUC = 0.976)
*All four reviewers. Locations: Abstract l.37–38; §3.10; §4.8 l.594–609; §5 l.693–702; Conclusion l.736–739.*

- **Reframe E10 as a *coherence/redundancy* check, not validation.** State plainly that the
  "held-out" features are geometrically coupled to the K-means inputs
  (`corr(compactness, aspect_ratio) = −0.865`; compactness-only AUC ≈ 0.96).
- **Lead with the honest result:** `reach_200m` — the only shape-orthogonal feature —
  ranks last (SHAP 0.005, univariate AUC ≈ 0.53). Report this as the finding:
  *the types carry no connectivity signal beyond geometry.*
- Remove AUC = 0.976 from the Abstract and Conclusion; remove "five independent lines of
  evidence" → at most two are genuinely orthogonal (size-frequency, connectivity), and
  connectivity is null.
- **Code:** add the correlation matrix and the per-feature univariate AUCs to NB06; emit a
  small table for the SI.

### C2. Remove or replace the fabricated E5 "literature reference zones"
*B (integrity-level), C, D. Locations: §3.4 l.275–278; §4.2 / Table 2 / Fig. 9 l.476–507; NB03 cells `1b8c6808`, `dd4cb613`.*

- The (compactness × AR) rectangles are hardcoded and attributed to papers (Seidl 2011,
  Jönsson 2012, Gardiner 2013, Schelhaas 2003) that **publish no such envelopes**, and the
  overlap is circular (AR defined the types).
- **Action:** delete E5 as "validation." Either (a) drop it entirely, or
  (b) **[Path B]** replace with overlap against a real external layer (Senf & Seidl 2021
  attribution; Metsäkeskus damage polygons) and report a confusion matrix.
- If any qualitative shape comparison is kept, soften to "qualitatively consistent with
  reported morphologies" and cite the *actual* quantitative basis or label the boxes
  illustrative.

### C3. Correct the spatial-extent claim — the data are regional, not national
*B (externally verified), D. Locations: Title; Abstract; §3.1 l.167–168; Figs. 1, 6, 9, 11.*

- The foundational source (Junttila et al. 2024) covers **Southeast Finland (~117,366 ha)**,
  not the country. Figs. 6/9/11 are rendered over a hardcoded `BBOX_3067` tile =
  **10.9% of clusters** (NB03 cell 15, NB04 cell 10, NB05 cell 11).
- **Action:**
  - State the true survey footprint and acquisition design in §3.1, with a coverage map.
  - Re-title / reword to remove "national-scale" unless coverage is proven.
  - For each clipped figure add one sentence: *"statistics computed on full dataset
    (n = 14,582); map clipped to a ~24×26 km display window."*
  - Add one genuine extent-wide figure (e.g., 25 km hexbin of Type-0 fraction).

### C4. The orientation null falsifies the "Wind-thrown" label — relabel
*All four. Locations: §4.3 l.509–519; §5 l.678–691.*

- Type 0: r = 0.009, p = 0.747 → no directionality, the one diagnostic fingerprint of
  windthrow is absent.
- **Action:** relabel Type 0 morphologically (**"elongated/dispersed"**) and Type 1
  (**"compact/aggregated"**); reserve agent names for **[Path B]** only.
- Move the orientation result from "surprising footnote" to a primary result, and state
  the metric caveat (cluster-level PCA axis at eps = 20 m on ~8 points cannot resolve
  storm-track alignment — see C4-link with the eps insight in M-killer below).

---

## TIER 2 — Major (substantive; mostly require re-running the pipeline)

### M1. Test discreteness: is k = 2 a real dichotomy or a cut on a continuum?
*D, A. §3.4, §4.2.*
- Silhouette decays monotonically (0.401, 0.307, 0.296, 0.279, 0.246, 0.245) — no elbow.
- **Action:** (a) 1- vs 2-component Gaussian mixture BIC on the three features;
  (b) Hartigan dip test on the first PC; (c) report the fraction of clusters within a
  thin band of the k=2 boundary. Soften language to "a dominant compactness–elongation
  gradient summarized as two endpoints" unless discreteness tests pass.

### M2. Replace the strawman permutation null
*A, D. §3.4 l.269–273; §4.2; Fig. 5.*
- Column-shuffling destroys all covariance → beaten by any correlated data.
- **Action:** test against a covariance-preserving null (single multivariate Gaussian
  matched to observed covariance) and/or a gap statistic vs a PCA-aligned reference.
  Report whether k=2 survives a null the covariance structure cannot trivially beat.

### M3. Close the reproducibility trail (blocking for any acceptance)
*All four. NOTEBOOK.md vs §3–§4.*
- Paper: silhouette 0.401, n = 3,422/11,160. NOTEBOOK: 0.418, 3,355/11,227. Notebooks
  flagged "NOT re-run." The CSVs match the paper; the log and stale code comments do not.
- **Action:**
  - Re-run the **entire** pipeline end-to-end from one clean environment (NB01→NB07).
  - Pin `requirements.txt` to exact versions; record seeds for every stochastic step
    (KMeans `n_init`/`random_state`, RF, permutation RNG, betweenness `k=500`).
  - Regenerate every figure/table from that single run; add a provenance line per figure.
  - Update NOTEBOOK.md's "core facts" table and delete the stale `0.418` / `3,355` comments.

### M4. Fix the connectivity section — it currently measures density, not connectivity
*D (domain), C. §3.9, §4.6, Fig. 11; NB05 cell 10.*
- "Reach" = connected-component size; 34.9% of clusters are singletons; median 3-vs-2 sits
  in the noise floor and tracks Type-1's higher density. Betweenness "sentinels/corridors"
  imply a flow process on an untimed, unweighted 200 m graph — unsupported.
- **Action:**
  - Report the singleton fraction explicitly.
  - Compare reach **within size-matched strata** (or model `reach ~ type + size`, report the
    *partial* type effect); use ego-network k-step reachability, not whole-component size.
  - Delete "sentinel/corridor/propagation" language (or relabel "geometric bridge clusters,
    no temporal validation"). Cut the section to a paragraph unless reframed as a flow test
    (**[Path B]**, two time slices).

### M5. Quantify and propagate upstream detection error
*C (uniquely). §3.1, §5.*
- The 698,223 detections come from a DL model reported with **no precision/recall**, and the
  `max_value` confidence field is never used. Detection error is likely spatially/forest-type
  structured — i.e., confounded with the signal.
- **Action:**
  - Report detection validation (precision/recall vs any ground truth, stratified by forest
    type and `max_value`); report the score distribution.
  - **Detection-perturbation bootstrap:** drop p% of points (FN model) and inject p% within
    clusters (FP model) at the estimated rate; recompute per-type Clark–Evans and type
    assignments; show the typology survives. Run the full pipeline at 2–3 `max_value`
    thresholds (mirror the eps sweep).

### M6. Fix the size-frequency logic error
*B. §4.2 l.459–465; §5 l.667–674.*
- Power law is *rejected* (lognormal preferred, p<0.001), then its exponent (α₀=3.00 vs
  α₁=2.23) is read as a contagion signature — invalid per Clauset et al. 2009.
- **Action:** drop the α-based mechanism claim, or recast as a descriptive lognormal-tail
  comparison with no mechanistic inference (note Type-1's heavier tail follows from its
  larger mean area by construction).

### M7. Caveat Clark–Evans at small n
*A, C. §3.3, §4.1, Table 1.*
- Median cluster = 8 trees; 57.5% have ≤8 points; CE is high-variance there, and DBSCAN's
  density floor + any detector NMS can manufacture CE > 1.
- **Action:** report CE with per-n confidence bands or restrict interpretation to n ≥ 30;
  attribute the R > 1 pattern partly to sampling/detection, not purely ecology.

---

## TIER 3 — Minor (quick, high-credibility fixes)

- **m1. Moran's I weights:** text says binary `1[d≤δ]` (§3.6 l.296–298); code uses
  row-standardized (`w.transform='r'`). **Fix the text.** Justify or drop the 2,000-point
  subsample (a sparse `DistanceBand` is tractable on 14,582 at ≤2 km; drop the 5 km column
  instead).
- **m2. Fractal dimension** is box-counted on the **convex hull**, not the points
  (`shape.py`, contradicts Eq. 5 §3.3). **Correct the methods text**; show D isn't degenerate
  at n≈8 (D-vs-n plot).
- **m3. `core_to_edge`** uses a fixed 1 m erosion (`shape.py:117`) → size proxy, not interior
  structure. State this; consider a size-relative buffer.
- **m4. `orientation` into RF** as a linear feature despite being axial — encode as
  (sin 2θ, cos 2θ) or drop (SHAP 0.034, low impact).
- **m5. CSI** (`severity.py`) appears in code but not the paper and is unvalidated — drop or
  validate against an impact measure.

---

## The corrected AR sensitivity protocol (sharpened by reviewer follow-up, Lindqvist)

**Problem:** §4.8's AR-filter "robustness" is not just negligible — it's *structurally
incapable* of testing its stated question, and the AR ≤ 10 cap is **entangled with the typology**:
the pipeline caps AR ≤ 10, then clusters *on* AR, then reads the surviving high-AR tail as an
ecological "Type 0." That truncates the load-bearing axis and interprets what's left as a class.
The cap also interacts with the eps artefact: at eps = 20 m, detection chaining manufactures
elongation, and AR ≤ 10 is simply the hand-set line between "artefact, discarded" (>10) and
"wind-throw type" (3.5–10). Leaving this as "uninformative as designed, re-cluster needed" and
**not** doing it defers the one test that can falsify the typology. **Do it or remove it.**

**Fix — re-cluster from raw across and BEYOND the boundary:**

1. Start from the **43,539 raw DBSCAN clusters** (pre-AR-filter), keep n ≥ 5 and area ≥ 100 m².
2. Apply AR caps ∈ **{5, 7, 10, 15, 20, ∞ (no cap)}**, **re-running shape metrics + K-means (k=2)
   de novo at each level** (apply the cap *before* K-means; anything on the already-capped CSV is
   void). The critical levels are 15, 20, ∞ — the ones the current test cannot reach.
3. Report per cap: n_clusters, Type-0 prevalence, silhouette; the **K-means centroid in
   (fractal_dim, CE, AR) space per type** (does the Type-0 centroid drift outward as the cap
   relaxes?); **membership churn** as adjusted Rand index vs. the AR ≤ 10 partition; and the
   **within-Type-0 AR distribution** (a slab pressed against the boundary is a truncation tell).
4. **Decision rule:** *Robust* = prevalence, Type-0 centroid, and ARI (≥ ~0.8) stable from AR ≤ 10
   through no cap → cap is harmless and Type 0 is a real mode. *Cap is doing the work* (predicted)
   = relaxing the cap repopulates the high-AR tail, prevalence climbs, the centroid migrates to
   higher AR, ARI degrades → the "elongated type" is the truncated edge of a continuous monotone
   distribution and cannot be presented as a class.
5. **Remove corridors by what they ARE, not by capping the clustering axis.** Excise road/utility
   linear artefacts with a *geometric* rule (PCA eigenvalue ratio + collinearity residual, or
   proximity to a roads/utility vector layer), leaving AR free for the typology.

## Reviewer-required methodology upgrades (beyond honest reframing)

Lindqvist's verdict: the reframes (gradient language, E5 demotion, AUC-as-coherence, regional
scope, orientation/eps foregrounding) are a **precondition, not the revision**. Honest framing
tells the reader the result *might* be an artefact; it cannot *show* it isn't. Two routes:

- **Route A — eps-stability gate (HARD REQUIREMENT, cheap, uses existing code).** Sweep
  eps ∈ {10,15,20,25,30,40} m; at each, test (i) k-selection stability (does k=2 keep winning on
  silhouette/CH/gap?), (ii) membership stability via consensus clustering / ARI across eps,
  (iii) eps-invariance of the two endpoint centroids in (compactness, AR). **Gate:** if the
  two-type structure is *not* eps-invariant, the "typology" framing must be dropped and the paper
  rewritten as a methods/cautionary note. This + the AR sweep are the two experiments that decide
  whether a real typology exists.

- **Route B — model the point pattern directly (what makes it a strong paper).** Replace
  "DBSCAN-then-shape" (one hard scale, then scale-inheriting shape metrics) with scale-explicit
  spatial point-process analysis:
  - **Inhomogeneous PCF g(r) / inhomogeneous Ripley K** with an intensity surface conditioned on
    host (MS-NFI spruce/forest-cover raster) — also the fix for the retired E4 PCF. Beetle
    contagion has a characteristic short-range g(r) peak *defined at a scale*, not at a hand-set
    eps; conditioning on host stops "compact patch" from being confounded with "wherever spruce is
    dense."
  - **Anisotropic / directional K-function** as the *principled* orientation test on the raw
    pattern (rigorously finds or rules out storm-track alignment, unlike a DBSCAN-derived axis).
  - **Marked point process with an external agent mark** (Metsäkeskus damage declarations, the
    Senf & Seidl 2021 disturbance-agent layer, or a spruce-dominance class) → tests whether
    morphology *predicts* an independently observed agent. This supervised step is the only thing
    that licenses the agent labels' return.

**Editor bottom line (Lindqvist):** accept the reframing as a gate, then *require* Route A (eps
stability) and the from-raw AR sweep as non-negotiable; steer toward Route B to recover the
disturbance-ecology claims.

## Four-reviewer follow-up — consolidated verdict (AR filter + reframe-vs-redo)

**UNANIMOUS on Q1 (AR filter):** delete-or-redo. All four call the current test "worse than
uninformative" / "a no-op presented as a control" / "misleading," and all four identify the
AR≤10 cap as entangled with the typology (it truncates the very axis k-means splits on; Type-0 ≡
high-AR). New evidence (Vasquez, on the data): **all 96 clusters with AR>8 are Type-0; Type-0's
99th-pct AR = 9.06** — a pile-up against the cap that already predicts a FAIL. Tanaka: this is
the *same leakage class* as the held-out RF — "selection on the defining feature of the target."

Merged AR protocol (consensus): re-cluster **from raw** across AR ∈ {5,7,10,15,20,**∞**},
re-fitting k-means each time (never reuse AR≤10 labels); the **uncapped arm is the gate**; report
prevalence, centroid drift, **ARI churn** vs the AR≤10 partition, and the within-Type-0 AR
distribution (pile-up check). **Decouple corridor removal from the elongation axis** (high-AR +
constant-width / low width-CV, or a Digiroad overlay — not a cap on AR). PASS iff prevalence
stable (±~3–5 pp through ∞), ARI ≥ 0.8, no spike at the cap, small centroid drift, no third mode.

**UNANIMOUS on Q2:** honest reframing is **necessary but NOT sufficient** — it fixes the
inferential overclaims but generates no new evidence. All four concur with the **eps-stability
hard gate** (sweep 10–40 m; drop "typology" if not eps-invariant; Tanaka: run eps×AR jointly once,
ARI-based criteria) and the from-raw AR sweep. Each adds one lane-specific REQUIRED test:

| Reviewer | Required addition (gate) | Recommended (not gate) |
|---|---|---|
| Lindqvist (ecology) | eps gate + AR-from-raw | Route B: host-conditioned inhomog. PCF + anisotropic K + marked process vs agent layer |
| Vasquez (stats) | **R3: directional/anisotropic 2nd-order statistic** (sector Ripley K / directional PCF) on within-cluster points — the first-order pipeline "cannot see anisotropy," which IS the wind-vs-beetle question | Route B — **explicitly recommend, do NOT gate** (needs host/agent rasters authors lack; would force re-introducing the retired PCF artifact) |
| Tanaka (RS/ML) | **R1 stratified detection precision/recall; R2 map spatial structure of recall; R3 Monte-Carlo detection-error propagation through whole pipeline; R4 regress gradient position on detector covariates** — reframing to a *gradient* is MORE exposed to detector confound, and `max_value` sits unused | host-conditioning is complementary but NOT a substitute for R1–R3 (it still assumes detected points = real points) |
| Anand (landscape) | **Abandon (not soften) the 200 m graph; a disturbance-REGIME claim is impossible from a single-year snapshot.** For a regime paper: (i) ≥2-yr multi-temporal FLOW test (does year-t predict year-t+1) [load-bearing]; (ii) dispersal-matched mechanism-specific graph (beetle ~250 m–1 km host-weighted kernel; wind storm-front tens-of-km directional — the two axes need INCOMPATIBLE graphs); (iii) host-conditioned marked process; (iv) connectivity rank-stability across eps & threshold | retitle to morphological characterization if multi-temporal data not acquired |

**The one genuine dissent:** Lindqvist steers toward Route B (point-process modeling) as the path
to a strong paper; **Vasquez explicitly pushes back on requiring it** — gate on eps-stability +
uncapped-AR + anisotropy (R1–R3) instead, because demanding a full inhomogeneous marked-process
model needs external data the authors credibly lack and would punish them for having responsibly
retired the unsound PCF. Tanaka adds that host-conditioning does not substitute for detection-error
validation. Net: **point-process modeling = strongly recommended trajectory; NOT a unanimous gate.**

**Consolidated hard gates for acceptance (where ≥3 reviewers agree):**
1. eps-stability gate (all four) — drop "typology" framing if not eps-invariant.
2. AR sweep from raw incl. uncapped, corridor filter moved off the AR axis (all four).
3. Anisotropy / directional second-order test of the wind-vs-beetle distinction (Vasquez; aligns
   with Lindqvist's anisotropic-K and Anand's directional-wind point).
4. Detection-error quantification + propagation, or a demonstration that recall is spatially flat
   (Tanaka; Anand's host-conditioning and Lindqvist's caveats are adjacent).
**Recommended (not gated):** full host-conditioned marked point-process (Route B) and ≥2-yr
multi-temporal flow analysis — the path to a genuine disturbance-regime contribution.

---

# FOUR-GATE EMPIRICAL VERDICT (executed; scripts in `analysis/`, results in `analysis/out/`)

**Headline: the discrete two-type "typology" fails 3 of 4 gates. Only a continuous
compactness–elongation gradient survives. The "Wind-thrown" label is empirically refuted.**

### Gate 1 — eps-stability: **FAIL** (`gate1_eps_stability.json`)
Re-clustered raw points at eps ∈ {10,15,20,25,30,40}. The supposed type is not eps-invariant:
- Type-0 prevalence slides monotonically with eps: 28.6% (10m) → 25.8% → 23.5% (20m) → 24.0 →
  23.4 → 22.9% (40m); n_clusters ranges 5,741 → 22,648.
- The Type-0 **centroid drifts continuously** with eps — AR 2.97 (10m) → 4.29 (20m) → 4.79 (40m),
  compactness 0.544 → 0.410 → 0.372. A genuine class would hold position; this rescales with the
  clustering radius, exactly the chaining-artefact signature Lindqvist predicted.
- BIC favours k>2 at **every** eps (e.g. eps=20: BIC₃ 109,481 < BIC₂ 109,977). No privileged k=2.
- Silhouette merely tracks eps (0.317→0.428) — bigger eps = rounder/larger blobs, not a better
  dichotomy.

### Gate 2 — AR-from-raw incl. uncapped: **FAIL (cap was creating the prevalence)** (`gate2_ar_sweep.json`)
Re-clustered from raw at eps=20, no AR cap, then swept caps {5,7,10,15,20,∞}, re-fitting each time.
- The AR≤10 cap removed only **62 clusters** (14,582 → 14,644 uncapped) but **the typology's
  headline prevalence is conditional on it**: Type-0% runs 33.1% (cap5) → 27.0 (cap7) → **23.5
  (cap10)** → 21.8 → 21.5 → **21.3% (∞)**. The published 23.5% is an artefact of where the cap fell.
- Type-0's AR distribution is a **pile-up against the cap**: at cap10 the within-Type-0 AR max is
  9.94 (95th pct 7.30) — pressed to the wall. Released, the tail extends smoothly to AR 14.96,
  18.48, 29.69 (caps 15/20/∞): a continuous monotone tail, **no natural break, no second mode**.
- ARI vs the cap-10 partition degrades monotonically as the tail is restored (1.00 → 0.915 → 0.902
  → 0.892) — the relabelling Vasquez/Tanaka flagged. (Confirms Vasquez's data note: all AR>8
  clusters are Type-0.) Verdict: AR is a continuum; the cap manufactured a clean "23.5% type."

### Gate 3 — anisotropy / directional 2nd-order: **FAIL (wind label refuted)** (`gate3_anisotropy.csv`)
Within-cluster point linearity, Type-0 vs Type-1, **matched on hull AR**:
- Overall Type-0 looks more linear (scatter-AR 3.80 vs 1.95) — but that is hull AR by construction.
- AR-MATCHED, the second-order linearity is **identical**: AR(2.0–2.5] T0=0.829/T1=0.830 (p=0.86);
  AR(3.0–4.0] T0=0.916/T1=0.915 (p=0.73); in AR(2.5–3.0] Type-1 is *higher* (0.877 vs 0.870).
- So once hull AR is controlled, Type-0 has **no extra directional/line-like (wind-swath) signal**.
  Combined with the earlier orientation null (r=0.009), the "Wind-thrown" interpretation is
  empirically unsupported. The first-order pipeline only ever saw hull shape, never anisotropy.

### Gate 4 — detection confound: **PARTIAL FAIL (confound is real, non-trivial)** (`gate_detection.py`)
Using the previously-unused `max_value` detection-confidence field:
- Detection covariates alone (mean/SD confidence, n, density) **predict type at AUC = 0.727** — far
  above chance. Type-0 vs Type-1 detection confidence differs (median 0.530 vs 0.563, p=2e-8).
- Modest false-negative perturbation flips type for **5.9% of clusters at 10% FN, 7.5% at 20%**
  (ARI 0.76 / 0.70) — past the ~10–15% Tanaka set as the falsification line at higher rates.
- Gradient-position vs confidence correlation is weak (ρ=+0.07), so the confound is partial, not
  total — but AUC 0.727 means the typology cannot be cleanly separated from detector behaviour
  without the stratified precision/recall (R1/R2) that needs external ground truth.

### Bottom line
Three independent hard gates (eps, AR-uncapped, anisotropy) each say the **discrete typology is an
artefact of clustering scale + a filter on its own defining axis, with no second-order content**;
the detector confound (AUC 0.727) compounds it. **What survives is a continuous, regional
compactness–elongation gradient** — which the revised manuscript already foregrounds. Required
edits below promote the gradient to the sole claim, delete the discrete "types" and the wind label,
and report the eps/AR/anisotropy/detection results as the evidence.

---

## The unifying insight to foreground (panel's two killer findings, merged)

Tie together what the panel surfaced separately:

- **eps = 20 m sets the elongation** (mean within-cluster NND = 7.79 m); Type-0 prevalence
  moves with eps (22.9% → 28.6%, §4.8) — a true ecological class would be eps-stable.
- **Aspect ratio is the typology's spine**, and the orientation test says that elongation
  carries **no directional information**.
- **Detection error is the latent common cause** that could shift all shape metrics together.

→ The honest thesis: there is a **robust morphological *gradient*** (compactness ↔ elongation),
but the evidence that it is **two discrete disturbance *types*** — rather than one continuum
shaped by clustering scale and detection geometry — is, on every independent axis tested
(reach, orientation, eps-stability), **absent**. Lead with the gradient; earn the types
(or the agents) only with Path B data.

---

## Execution checklist (suggested order)

1. [ ] **Decide Path A vs B** (governs scope of C2/C3/C4 and M4/M5).
2. [ ] **Re-run pipeline clean, pin versions + seeds** (M3) — unblocks all number reconciliation.
3. [ ] Corrected AR sweep from raw clusters (AR protocol) + eps sweep already exists.
4. [x] **Discreteness tests (M1) + covariance-preserving null (M2)** — DONE, see results below.
5. [ ] Detection-perturbation bootstrap + threshold sweep (M5); report detection metrics.
6. [ ] Rework connectivity (M4); fix size-frequency (M6); caveat Clark–Evans (M7).
7. [~] Reframe E10 (C1) and E5 (C2); relabel types (C4); correct extent (C3) — Abstract,
       Contributions, Methods §3.4, Results §4.2 done; Discussion/Conclusion/figures pending.
8. [ ] Minor text/code fixes (m1–m5).
9. [~] Rewrite Abstract, Contributions, Conclusion around the *gradient* thesis — Abstract +
       Contributions done; Conclusion pending.
10. [ ] Add SI: correlation matrix, univariate AUCs, seeds/versions, provenance per figure.

---

## Execution log

### M1 + M2 results (script: `analysis/discreteness_tests.py`, seed 42, run on `clusters_typed.csv`)

- **Reproduction:** k=2 silhouette = **0.4007** (paper reports 0.401 ✓ — the shipped CSV matches the manuscript, not the stale NOTEBOOK 0.418).
- **M2 covariance-preserving null (999× single multivariate Gaussian):** observed 0.4007 vs
  null mean 0.321 / 95th pct 0.324 / max 0.328 → **p = 0.001**. The typology **passes** the
  strong null — it is *not* merely an artifact of correlated marginals. (This rebuts the
  reviewers' assumption that any correlated data would beat the test; report it as a strength.)
- **M1a GMM BIC** (k=1..7): 116253, 109983, 109485, 108518, 107391, 107297, 107354 →
  optimal **k=6**; 2≻1 but every k>2 is better still → **no privileged 2-way split**.
- **M1b:** PC1 explains 60% of variance; **Sarle bimodality coefficient 0.486 < 0.555 → unimodal**.
- **M1c:** only 2.3% / 4.6% / 10.0% of clusters within 5/10/20% of the k=2 boundary.

**Verdict for the paper:** genuine structure (beats covariance-preserving null) **but a
continuum, not a dichotomy** (BIC favors >2 components; dominant axis unimodal). Frame k=2 as an
interpretive summary of a compactness–elongation gradient. *(Optional: `pip install diptest` to
add a formal Hartigan dip test to the SI; Sarle BC + GMM BIC already give a coherent picture.)*

### Manuscript edits applied (main.tex)
- **Abstract** rewritten: regional (not national) footprint; gradient-not-dichotomy; covariance-
  preserving null (p=0.001) + BIC/unimodality caveat; RF demoted to coherence check (connectivity
  null); agent labels withheld as hypotheses.
- **Contributions** (Intro) rewritten to the three honest claims; closing paragraph now withholds
  agent labels and previews the orientation-null caveat.
- **Methods §3.4** permutation test replaced with the covariance-preserving null + discreteness
  diagnostics; E5 "ecological coherence" claim removed.
- **Results §4.2** permutation paragraph replaced with the covariance-null result and a new
  "Continuum, not dichotomy" paragraph carrying the BIC/unimodality numbers.

### ⚠ Figures now inconsistent with the text (regenerate before submission)
- `fig_silhouette_permtest.png` still shows the OLD column-shuffle null (mean 0.244, 95th 0.245).
  Regenerate from the covariance-preserving null (mean 0.321, 95th 0.324). Add a small
  discreteness panel (BIC vs k; PC1 histogram) to accompany the new §4.2 paragraph.

### Session 2 edits applied (prose + figure consistency, items 1–4)
- **Discussion** fully rewritten: agent-attribution paragraphs replaced by "morphological
  gradient, not an agent typology" + the convergent internal-evidence argument (orientation null
  + eps-instability + connectivity null = geometry, not agent); size-frequency contagion claim
  withdrawn; held-out-AUC paragraph corrected; Limitations now lead with detection-error
  propagation + regional extent; Future work reframed around external/multi-temporal validation.
- **Conclusion** swept (done session 1).
- **C4 relabel:** Type 0 → "elongated/dispersed", Type 1 → "compact/aggregated" in Table 1
  header, §4.2 endpoint descriptions, Figs. 4/6 captions; agent names removed.
- **C2 E5 demotion:** §4.2 shape-space paragraph, Table 2, Fig. 9 caption recast as
  illustrative/author-defined + explicitly flagged circular; no longer "validation".
- **C1 tail:** §4.7 retitled "Feature-Space Coherence Check"; connectivity §4.6 softened
  (singleton fraction noted, "sentinel/corridor" removed, reach = density correlate);
  "national" → "survey/study-region display window" in Moran/LISA/contagion captions.
- **AR sensitivity:** §4.8 now states the test is uninformative as designed (runs on
  already-capped data) and flags the re-cluster-from-raw protocol as needed.
- **Figure:** `fig_silhouette_permtest.png` regenerated as a 3-panel (covariance null + BIC + PC1)
  via `analysis/make_permtest_figure.py`; caption updated to match.

### Session 4 — Tier-1 consistency + Tier-2 minor fixes (DONE)
- **Title** → "A Morphological-Gradient Characterisation of Dead-Tree Clusters in Southeastern
  Finland's Boreal Forests" (drops "typology" + "national").
- **§3.1 Dataset** integrity fix: "entire Finnish land area / full national footprint" → regional
  (SE Finland), prevalences flagged as not-national; Junttila cite added at source.
- Swept residual "national" claims (Intro l.70/84, Fig.1 caption, orientation l.592).
- **m1 Moran weights:** text "binary $1[d\le\delta]$" → row-standardised distance-band (matches code).
- **m2 fractal dim:** text now says box-counting on the convex hull (not raw points) + small-n /
  descriptive-index caveat.
- **m3 core_to_edge:** text corrected from "fraction of core points" to the actual area-erosion
  ratio (1 m buffer) + size-proxy caveat.
- **m4 orientation in RF:** added note it is axial data entered linearly (low impact); also added
  the GroupKFold block-size-imbalance caveat (5 blocks ≈ 80% of clusters).
- **NOTEBOOK.md** "Core dataset facts" table updated to post-rerun values (3,422/11,160, sil 0.401,
  CE 2.04/1.68) + regional footprint + post-gate verdict note; historical log entries preserved.
- **m5 CSI:** `src/severity.py` is unused by the manuscript (dead code relative to the paper).
  NOT deleted (user-authored, no paper dependency) — flagged here; remove or validate against an
  impact measure if it is ever surfaced in the text.

### Still to do (heavier, re-run dependent — session 5)
- Full clean pipeline re-run + version/seed pinning (M3) to reconcile all numbers.
- Corrected AR sweep from the 43,539 raw clusters (M-AR protocol); add to NB07.
- Detection-error metrics + perturbation bootstrap + max_value threshold sweep (M5).
- Size-matched / partial-effect connectivity re-analysis (M4); drop the 5 km Moran subsample (m1).
- Minor code/text: Moran weights description (m1), fractal-dim-on-hull text (m2),
  core_to_edge buffer (m3), axial orientation encoding (m4), drop CSI (m5).
- SI: correlation matrix, univariate AUCs, seeds/versions, per-figure provenance.
