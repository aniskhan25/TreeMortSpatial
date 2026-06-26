# Paper Plan — *Is it the ecology or the detector? A falsification harness for spatial-pattern inference on deep-learning-derived maps*

*(working title; "Tanaka's paper" — the RS/ML methods spin-out of the TreeMort Contribution 2)*

---

## 1. The gap and the one-sentence thesis

Ecologists increasingly compute **spatial statistics** (clustering, autocorrelation, point-process
structure, morphometrics, connectivity) on **maps produced by black-box detectors** — U-Net / Mask
R-CNN segmentations of trees, crowns, dead wood, animals, burn scars, kelp, nests. The detector is
treated as a transparent window onto the landscape. It is not: recall varies in space, object
size/shape couples to detection confidence, a confidence threshold manufactures structure, and
image-tile seams inject anisotropy. **There is no standard protocol that tells you whether a
published spatial-ecological pattern is a property of the ecosystem or a footprint of the
instrument.**

> **Thesis.** A second-order spatial-ecological signal on detector-derived maps is only credible
> after it survives a defined battery of falsification tests; we formalise that battery, prove on
> simulation that it catches instrument artifacts with known ground truth, ship it as an open package,
> and demonstrate it on real datasets — including one where it overturns a plausible published-style
> claim.

This is a **methods + software** paper. It generalises the one-dataset harness already executed in the
TreeMort study into a reusable, dataset-agnostic standard.

---

## 2. Why this is novel and fundable

- The validation literature stops at the **detector** (precision/recall, mAP, IoU). It almost never
  asks what detector error does to the **downstream spatial statistic** an ecologist actually reports.
- Spatial-ecology methods papers assume the point pattern is observed cleanly; they model edge effects
  and intensity, not **observation error with spatial and morphological structure**.
- The contribution is a *bridge*: a taxonomy of confound mechanisms + a falsification protocol +
  a simulation testbed + software, validated end-to-end. Nobody has packaged this.
- It is **data-light**: the methodological core runs on simulation (ground truth known by
  construction), so it does not depend on acquiring new field data — the chronic blocker.

---

## 3. The confound taxonomy (the conceptual backbone)

Five mechanisms by which a detector fabricates or destroys spatial-ecological signal. Each maps to a
test and to a simulation knob.

| # | Confound mechanism | Fake signal it creates | Test (battery item) |
|---|---|---|---|
| C1 | **Spatially heterogeneous recall** (detection prob varies with terrain, illumination, density, latitude) | spurious intensity gradients & positive autocorrelation | T2 spatial-recall map; T3 propagation |
| C2 | **Confidence–morphology coupling** (object size/shape predicts detection confidence) | spurious morphological "classes"/gradients | T4 confound regression (AUC/ρ) |
| C3 | **Threshold dependence** (the `score > τ` cutoff) | structure that appears/dissolves with τ | T5 threshold sweep |
| C4 | **Tiling / acquisition seams** (image-strip boundaries, flight lines) | spurious anisotropy & edge clustering | T6 seam-anisotropy test |
| C5 | **Detector prior / class imbalance** (the model's own footprint) | "ecology" that is the training distribution | T1 stratified P/R; null models |

---

## 4. The falsification battery (T1–T7)

Each test yields a **verdict** (survives / caveat / withdraw) with a quantitative threshold. Items in
**bold** already exist as TreeMort scripts and become reference implementations; the rest are new and
generalised.

- **T1 — Stratified precision/recall.** P/R conditioned on covariates (density, terrain, class).
  *Needs validation labels;* where absent, supplied by the simulation engine (§5) or a held-out
  hand-labelled tile. Output: is recall flat in the covariates that also drive the ecological claim?
- **T2 — Spatial recall surface.** Estimate recall as a smooth field \(r(u)\); test for spatial
  trend/autocorrelation. A non-flat \(r(u)\) is a *generator* of fake first- and second-order signal.
- **T3 — Monte-Carlo error propagation** *(generalises `gate_moran_propagation.py`).* Perturb the
  detected set by a fitted error model (false negatives at rate \(1-r(u)\), false positives, jitter)
  and re-compute the target statistic across realisations. Verdict on whether **sign and magnitude
  survive**. TreeMort precedent: Moran sign survived for Clark–Evans/AR/compactness, **failed for
  fractal dimension**.
- **T4 — Detector-confound regression** *(generalises `gate_detection.py`, `scope_experiments.py`).*
  Can detector covariates (mean confidence, size, local density) predict the ecological label/gradient
  position? Report discriminative AUC and the coupling ρ between the *continuous* statistic and
  confidence; partial out spatial covariates (the latitude-conditioned R4). TreeMort precedent:
  type predictable at **AUC 0.727**, continuous structure ρ=0.07, confound survives latitude
  (0.727→0.720).
- **T5 — Threshold-stability sweep.** Re-run the whole pipeline across detection thresholds τ; a real
  signal is τ-stable, an artifact rescales. (Parallels the **eps-stability** logic of
  `gate_eps_ar.py`, applied to the *detector* knob rather than the clustering knob.)
- **T6 — Seam/anisotropy test** *(generalises `gate_anisotropy.py`).* Directional second-order
  statistics referenced to the tiling/flight-line geometry; distinguishes ecological anisotropy from
  acquisition-grid anisotropy.
- **T7 — Substrate-conditioning** *(generalises the host-/species-conditioning in
  `host_pointprocess.py`, `species_validation.py`).* Condition the statistic on an **external substrate
  the detector cannot know** (here forest host/species); residual structure beyond substrate is the
  part not attributable to a substrate-correlated detection bias. The general principle: *find a
  covariate that drives the ecology but not the detector, and condition on it.*

Plus the **non-detector** nulls already built — covariance-preserving null + discreteness diagnostics
(`discreteness_tests.py`) and leave-one-block-out CV (`gate_block_cv.py`) — folded in as the
"is the structure even there before we worry about the detector?" pre-filter.

---

## 5. The simulation engine (methodological heart — new code)

The piece that makes this a *methods* paper rather than a checklist: a generator with **known ground
truth** so every test can be shown to have the right operating characteristics.

1. **Truth process.** Simulate a point/object field with *prescribed* second-order structure
   (CSR; inhomogeneous Poisson with a chosen intensity; Thomas/Matérn cluster; LGCP; anisotropic).
2. **Detector model.** Apply a parametric observation layer:
   - recall surface \(r(u)\) (flat → strongly trended);
   - confidence model \(p(\text{detect}\mid \text{size, shape, density})\) (C2 coupling strength as a knob);
   - threshold τ; false-positive rate; localisation jitter; tile-seam dropout.
3. **Two readouts per scenario:** (a) the **naive** spatial statistic (what a paper would report) and
   (b) the **harness verdict.**
4. **Operating characteristics.** Sweep confound strength from 0 → high and plot:
   - false-positive rate of the *naive* statistic (fake signal vs confound strength);
   - the harness's **sensitivity/specificity** in flagging it (ROC of "survives" vs truth);
   - power loss when the harness is applied to a *genuine* signal (does it over-reject?).

This yields the paper's killer figure: *naive inference fabricates significant structure once
confound strength crosses X; the harness flags it; and it does not destroy genuine signal.*

---

## 6. Demonstration datasets (generality — pick ≥2)

- **D1 — TreeMort (in hand).** The full worked example already exists: AUC-0.727 confound, Moran
  propagation, host/species conditioning, eps/threshold/anisotropy gates. Becomes the anchor case
  study and a reference implementation. *No new work to produce the results — only to re-present.*
- **D2 — A public detector-derived dataset *with* validation labels** (so T1/T2 run for real). Candidates:
  NEON crown maps, a published global/again-public tree-crown or canopy-mortality detection set, or a
  segmentation benchmark with held-out masks. Selection criterion: raw detections **+ confidence** +
  some ground truth + a plausible spatial-ecological question.
- **D3 (optional, for breadth) — a non-forest detector map** (camera-trap/aerial animal counts, burn-scar
  or kelp segmentation) to show the protocol is taxon/sensor-agnostic.

Recommendation: **D1 + D2** for the first submission; D3 as reviewer-proofing breadth if time allows.

---

## 7. Software deliverable

`detectorcheck` (or similar) — a small Python package:
`run_battery(points, confidence=None, covariates=None, window=None, truth=None, tests=[...])`
→ a per-test verdict table + diagnostic plots. Reuses `figstyle.py`. Ships the simulation engine as
`detectorcheck.simulate`. Open-source (the venue below expects it). This is what makes the paper get
*used* and *cited*.

---

## 8. Figures (target ~5 main)

1. **Confound taxonomy schematic** (C1–C5: how each fabricates a pattern) — conceptual.
2. **Simulation operating characteristics** — naive false-positive rate & harness ROC vs confound
   strength (the killer figure, §5).
3. **TreeMort anchor** — the four-gate + confound-regression panel (largely the existing
   `fig_dbscan_sensitivity` / `fig_anisotropy` / confound results, re-laid-out).
4. **D2 case study** — battery verdict table + spatial-recall surface for the second dataset.
5. **Decision flowchart** — the protocol as a referee/author checklist (survives / caveat / withdraw).

---

## 9. Target venue

- **First choice: *Methods in Ecology and Evolution*** — exactly its remit (a reusable method + open
  software + simulation validation + ecological case studies).
- Alternatives: *Remote Sensing of Environment* (if framed RS-first, heavier on detector modelling),
  *Ecological Informatics*, or *ISPRS J. Photogrammetry & RS*.

---

## 10. Relationship to the current (TreeMort) paper — avoiding self-overlap

- TreeMort *uses* the harness on one dataset and reports an **ecological** result (host-conditioned
  mortality structure). This paper *is* the harness: general protocol + simulation proof + software +
  multi-dataset.
- Cite TreeMort as the motivating application; lift only the **method**, not the ecological findings.
- Net new for this paper: the **confound taxonomy**, the **simulation engine + operating
  characteristics**, the **second dataset**, and the **package** — none of which exist yet.
- Authorship/credit: the harness framing is "Tanaka's"; confirm real-world author roles before drafting.

---

## 11. Risks and mitigations

| Risk | Mitigation |
|---|---|
| No public dataset with detections **+ confidence + ground truth** | Lead with the simulation engine (self-contained truth); D2 can use a hand-labelled held-out tile for T1/T2 only |
| "Just a checklist" reviewer reaction | The simulation operating-characteristics result + the case where it overturns a claim make it a *result*, not a checklist |
| Overlap-with-TreeMort / self-plagiarism | Strict method-only lift; different figures; TreeMort cited as application |
| Scope creep (too many datasets/tests) | Freeze at T1–T7 + D1–D2 + simulation for v1; D3 and extra detectors are explicitly "future" |
| Detector-error model too bespoke | Keep the observation model parametric and documented; report sensitivity to its assumptions |

---

## 12. Milestones (suggested)

1. **M1 — Scoping & lit map** (the gap, the 10–15 papers that compute spatial stats on detector maps
   without validating). Decide D2.
2. **M2 — Simulation engine + operating-characteristics experiment** (the heart; self-contained).
3. **M3 — Refactor the seven TreeMort gate scripts into the `detectorcheck` package** with a clean API.
4. **M4 — D2 acquisition + battery run** (and D3 if pursued).
5. **M5 — Draft** (taxonomy → battery → simulation → case studies → decision protocol), figures, package release.
6. **M6 — Internal review, then submit to MEE.**

*(No calendar dates attached — sequence only; slot into your schedule.)*

---

## 13. Compute & data footprint

- Simulation + battery: laptop-scale (the TreeMort gates already run in seconds–minutes; spatstat
  envelopes are the heaviest piece).
- D2/D3: one or two public rasters/detection tables; no new field campaign.
- Reuses the existing `.venv` (Python 3.9, scikit-learn, rasterio, pointpats, geopandas) and
  R/spatstat already installed and version-pinned in the TreeMort SI.

---

## 14. Open decisions to confirm before drafting

1. **Scope/ambition:** full framework + package + two datasets (recommended), **or** a lighter
   "methods note" = simulation engine + TreeMort case only (faster, lower-impact).
2. **Second dataset (D2):** which public detector-derived map (NEON crowns? a canopy-mortality set?).
3. **Venue:** MEE (recommended) vs RSE.
4. **Authorship/credit** for the "Tanaka" framing.
