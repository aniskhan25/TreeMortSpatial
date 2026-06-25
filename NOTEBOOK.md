# Development Notebook

Running log of changes, analyses added, and decisions made.
Most recent entry first.

---

## 2026-06-09 — Expert panel review + analytical corrections

### Overview

Multi-agent expert panel review (forest sciences, statistical ecology, ML/data analytics)
identified five analyses that were statistically invalid or infeasible, plus three code bugs.
All corrections were applied in this session. Notebooks have NOT been re-run yet.

### Analyses retired (cells replaced with retirement markdown, code cells deleted)

| Analysis | Notebook cell | Reason |
|---|---|---|
| E3 — Temporal wave-front | NB06 m2 header, b4 code | `alku_aika` is flight-strip timestamp (same for all trees in a strip), not per-tree mortality date. Spearman r ≈ 0.37–0.47 (both types) reflects south-to-north survey aircraft flight schedule, not ecological wave. |
| E4 — Pair Correlation Function | NB04 cff14401 header, 61578819 code | CSR envelope over national bounding box (~320,000 km²) gives g(r)=400–525 at 50 m — artefact of heterogeneous forest cover, not a within-type signal. Requires inhomogeneous PCF with VMI forest cover raster (not available). |
| E7 — Exponential spread kernel | NB04 802ccb49 header, 561c934a code | Both types share the same national KDE density peak → both kernels fitted from the same origin point → mathematically identical 1/λ ≈ 6 km by construction. No biological inference possible. |
| E8 — Within-cluster date range | NB06 m1 header, b3 code | Within any 20m-radius DBSCAN cluster, all trees come from the same flight strip → date_range = 0 days for 14,581/14,582 clusters by construction. Mann-Whitney p=0.58 compares two Dirac deltas at zero. |

### Code bugs fixed

**E2 Rayleigh test (NB03 cell b11102c0) — always returned r=0.000, p=1.0:**
- Root cause: `np.concatenate([angles, angles + 180]) % 360` creates (α, α+180°) pairs that
  exactly cancel in the circular resultant vector, guaranteeing r=0 regardless of data.
- Fix: `doubled_rad = np.deg2rad((2 * angles) % 360)` for statistics; mirror approach kept
  for rose diagram display only; mean direction halved back to [0°,180°] after computing
  from doubled angles.

**Clark-Evans edge bias (`src/metrics/spatial.py`):**
- Root cause: `expected = 0.5 * sqrt(area/n)` ignores boundary effects; trees near cluster
  edge have inflated observed NND, biasing CE upward (especially for the majority 5–20 tree clusters).
- Fix: Donnelly (1978) edge correction added:
  `expected = 0.5√(A/n) + (0.0514 + 0.041/√n) · P/n`
  where P = convex hull perimeter (= `hull.area` in scipy 2D ConvexHull, where .area=perimeter).
- Impact: CE values will decrease modestly (~10–15%); Type0 vs Type1 ordering preserved.
  **NB02 through NB06 must be re-run** to regenerate clusters_with_metrics.csv.

**E10 tautological AUC=0.999 (NB06 cells m3, b5, b6) — redesigned:**
- Root cause: RF trained to predict K-means labels using the same features (fractal_dim,
  clark_evans, aspect_ratio) that defined those labels. AUC=0.999 is mathematically
  guaranteed by the K-means Voronoi boundary.
- Fix: E10 redesigned to use ONLY held-out features (compactness, alpha_compactness,
  core_to_edge_ratio, nnd_mean, nnd_std, orientation, reach_200m) — none used in K-means.
- Added: Spatial block cross-validation (100×100 km grid in EPSG:3067) to prevent
  autocorrelation leakage across the ≥5 km autocorrelation range.
- Interpretation: AUC > 0.7 = morphological types structured in independent feature space.

**E6 subsampling and multiple comparisons (NB04 cells 7428b5f4, 90a763ec):**
- Random n=2000 subsample replaced with stratified proportional subsample
  (n0=int(2000×3355/14582), n1=2000-n0) to prevent geographic overrepresentation.
- Bonferroni correction added: α*=0.05/20=0.0025 across 5 thresholds × 4 metrics.
  Cell annotation: `*` = Bonferroni significant; `†` = nominally significant only.

### Silhouette permutation test added (NB03)

New cell (id=b8c1bba9, inserted after b5 K-means cell): 999-permutation test comparing
observed silhouette=0.418 against null distribution from column-permuted features.
Tests whether the K-means partition exceeds random feature covariance structure.
New figure: `fig_silhouette_permtest.png`.

### Sensitivity analysis notebook added (NB07)

New notebook `07_sensitivity.ipynb` with two tests requiring no external data:
- **S1 — DBSCAN eps sensitivity**: eps ∈ {10, 20, 40} m; re-runs DBSCAN + shape metrics +
  K-means k=2; tracks n_clusters, silhouette, Type0%, mean AR per type. ~10–20 min to run.
- **S2 — AR filter sensitivity**: AR cap ∈ {7, 10, 15}; uses existing clusters_with_metrics.csv
  (no re-computation); tests whether AR≤10 truncates the Type0 elongation distribution.
New figures: `fig_dbscan_sensitivity.png`, `fig_ar_sensitivity.png`.

### Run order and prerequisites

**Before re-running notebooks:**
- NB02: CE values change (Donnelly correction) → regenerates clusters_with_metrics.csv
- NB03: E2 now gives a real orientation result; permutation test is new
- NB04: E6 with stratified subsampling + Bonferroni
- NB05: depends on clusters_typed.csv from NB03
- NB06: redesigned E10 with held-out features + GroupKFold spatial CV
- NB07: S2 depends on clusters_with_metrics.csv from NB02 (S1 needs trees_clustered.csv from NB01)

---

## 2026-06-09 — Type label revision

### Motivation

E5 (shape-space validation) showed the original labels were misassigned:
- Old Type0 = "Pest/Fire" — but 57.8% of Type0 falls in the **wind-damage** zone (AR 3.5–10,
  compactness 0.05–0.45), not the fire or beetle zone.
- Old Type1 = "Uniform/Scattered" — but 94.4% of Type1 falls in the **fire/beetle compact** zone
  (AR 1–3.5, compactness 0.40–0.85). Only 1.3% is in the wind zone.

Additional supporting evidence from E1: Type1 has the heavier size-distribution tail (α=2.23 vs
Type0 α=3.00), consistent with the self-reinforcing aggregation of bark beetles (Senf & Seidl 2021).

### New labels

| Key | Old label | New label | Morphological basis |
|---|---|---|---|
| Type0 | Pest/Fire | **Wind-thrown** | AR≈4.3, compactness≈0.41, CE=2.71; elongated, smooth, regular spacing; 57.8% in wind zone |
| Type1 | Uniform/Scattered | **Beetle/Patch** | AR≈2.0, compactness≈0.65, CE=2.02; compact, complex, aggregated spacing; 94.4% in fire/beetle zone |

### Files changed

`morph_type` column values ("Type0", "Type1") are **unchanged** throughout the pipeline — only
the human-readable display strings were updated in:
- `03_typology.ipynb` cells b5 (print statements + comment), b6 (legend `label` dict),
  b11102c0 (orientation rose label list)
- `04_spatial_patterns.ipynb` cell 61578819 (PCF label list)
- `06_disturbance_ecology.ipynb` cell m1 (markdown description)

---

## 2026-06-08 — Literature-grounded enhancements (E1–E10)

### What changed

Added 10 new analytical sections across notebooks 02–06, each grounded in a specific
published literature gap. A new notebook `06_disturbance_ecology.ipynb` was created for
synthesis analyses. `src/metrics/spatial.py` gained `pair_correlation_function`.

New packages added to `requirements.txt` and the venv:
- `powerlaw` — MLE power-law fitting and likelihood ratio tests (Clauset et al. 2009)
- `pycircstat` — circular statistics for axial orientation data (requires `nose`)
- `shap` — TreeExplainer for interpretable RF feature importance
- `networkx` — proximity graph construction and betweenness centrality

**Pipeline change**: NB02 now computes `date_range` and `date_mean_doy` per cluster
(from `alku_aika` in `trees_clustered.csv`) and merges them into `clusters_with_metrics.csv`.
NB05 appends `reach_{d}m` and `betweenness_200m` to `clusters_typed.csv`. NB06 reads
this enriched CSV for E3 and E10.

**Run order is now strictly**: NB01 → NB02 → NB03 → NB04 → NB05 → NB06.

### Enhancements detail

| ID | What | Notebook | New figure |
|---|---|---|---|
| E1 | Power-law vs log-normal CCDF per type | NB03 (after K-means) | `fig_size_dist.png` |
| E2 | Orientation rose diagram + Rayleigh test | NB03 (after K-means) | `fig_orientation_rose.png` |
| E3 | Temporal wave-front Spearman r (date_mean_doy vs distance) | NB06 | `fig_wavefront.png` |
| E4 | Pair Correlation Function with 99 CSR envelopes | NB04 | `fig7d_pcf.png` |
| E5 | Shape-space scatter vs published disturbance zones | NB03 (after K-means) | `fig_shape_space.png` |
| E6 | Multi-threshold Moran's I heatmap (5 thresholds × 4 metrics) | NB04 | `fig_moran_heatmap.png` |
| E7 | Exponential spread kernel ρ(r) = A·exp(−λr) from KDE hotspot | NB04 | `fig_spread_kernel.png` |
| E8 | Within-cluster acquisition date range violin by type | NB02 + NB06 | `fig_date_range_violin.png` |
| E9 | Contagion network reach + betweenness bridge map | NB05 | `fig_contagion.png` |
| E10 | Random Forest + SHAP discriminant synthesis | NB06 | `fig_shap.png` |

### Key technical notes

- **E2 axial data**: `orientation` is 0–180° (PCA major axis, not a direction).
  Doubled to 0–360° before Rayleigh test so the circular mean is meaningful.
- **E6 subsampling**: Moran's I at 5000m threshold on 14,582 points requires building a
  dense weight matrix; 2000-point subsample used for feasibility (noted in figure title).
- **E7 hotspot**: KDE computed in EPSG:3067 directly (metres), so exponential kernel
  distances are already in metres — no reprojection needed.
- **E9 betweenness**: Approximate betweenness centrality with k=500 random sources
  (networkx `betweenness_centrality(G, k=500, seed=42)`) — exact computation would be
  O(n·m) and too slow for sparse but large graphs.
- **E10 interpretation**: morph_type was defined by K-means on morphometric features,
  so the RF partly reconstructs that decision boundary by construction. The meaningful
  signal is whether `date_range` and `reach_200m` rank high in SHAP — ecological
  correlates independent of the shape definition.
- **`alku_aika` format**: `"14.6.2023 08.10.00"` — dots used in both date AND time parts.
  Parsed with `pd.to_datetime(format='%d.%m.%Y %H.%M.%S', errors='coerce')`.
  This is the flight-strip acquisition start time, not individual tree detection time.

---

## 2026-06-07 — Map figure fixes (contextily basemap slant)

### Problem

Figures in NB01, NB03, NB04, NB05 were visually slanted. Root cause: using
`ctx.add_basemap(ax, crs='EPSG:3067')` with raw matplotlib scatter plots. Contextily
reprojects Web Mercator tiles into EPSG:3067, which is a transverse Mercator centred
at 27°E — grid north diverges from true north near Finland's longitude extremes (~20°E
and ~32°E), causing visible rotation.

### Fix applied

Replaced all map-producing cells with the pattern confirmed in the original analysis
notebook (`shape_analysis_v1.ipynb`):

```python
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['cx'], df['cy']),
                        crs='EPSG:3067').to_crs(epsg=3857)
fig, ax = plt.subplots(figsize=(10, 10))
gdf.plot(ax=ax, ...)
ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)  # NO crs= argument
```

Key: reproject data to EPSG:3857 (Web Mercator), let geopandas set the axes aspect
ratio, then add basemap without `crs=` — tiles are already in their native projection.

### Packages added

`contextily`, `matplotlib-scalebar`, `pyproj` added to `requirements.txt`.
Installed via `pip3 install --break-system-packages` on the system Python; also
available in the venv.

### Special case: KDE figure (fig7a)

The `contourf` KDE cannot use `gdf.plot()`. Solution: extract 3857 coordinates from
the geopandas geometry after reprojection, compute KDE on those, then overlay
`gdf_pts.plot()` for the scatter layer. All coordinate spaces consistent in 3857.

---

## Core dataset facts (reference)

**Updated to the post-Donnelly-rerun pipeline state (matches `clusters_typed.csv` and the
manuscript). The earlier values in this table (3,355 / 11,227 / silhouette 0.418 / CE 2.71/2.02)
were pre-rerun and are retained only in the dated historical log entries above.**

| Item | Value |
|---|---|
| Raw trees | 698,223 dead trees, EPSG:3067 |
| Survey footprint | National-scale sample: 312 aerial tiles across Finland, ~60.4–67.8°N, 21.6–29.6°E (south boreal to Lapland), discrete latitudinal bands — NOT contiguous coverage, NOT southeastern-only |
| Raw clusters (DBSCAN eps=20m, min_samples=3) | 43,539 |
| Filtered clusters (≥5 trees, area ≥100m², AR ≤10) | **14,582** |
| Type0 — elongated/dispersed endpoint (AR≈4.29, compactness≈0.41, CE≈2.04) | 3,422 (23.5%) |
| Type1 — compact/aggregated endpoint (AR≈2.01, compactness≈0.65, CE≈1.68) | 11,160 (76.5%) |
| K-means silhouette (k=2) | 0.401 |
| Alpha compactness NaN (alpha shape failure) | 995 clusters |
| Acquisition year | 2023 |
| CRS of source data | EPSG:3067 (ETRS89/TM35FIN) |

**Note (post-gate verdict):** the two "types" do NOT survive the four stability gates
(eps-instability, AR-cap dependence, no anisotropy, detector confound — see
`report/REVISION_PLAN.md`). The defensible result is a continuous compactness–elongation
gradient; "Type0/Type1" are labels for its two ends, not disturbance agents.

---

## Key literature (full references)

| Short ref | Full citation |
|---|---|
| Senf & Seidl 2021 | Senf, C. & Seidl, R. (2021). Mapping the forest disturbance regimes of Europe. *Nature Sustainability* 4, 63–70. |
| Seidl et al. 2017 | Seidl, R. et al. (2017). Forest disturbances under climate change. *Nature Climate Change* 7, 395–402. |
| Seidl et al. 2011 | Seidl, R. et al. (2011). Unraveling the drivers of intensifying forest disturbance regimes in Europe. *Global Change Biology* 17, 2842–2852. |
| Raffa et al. 2008 | Raffa, K.F. et al. (2008). Cross-scale drivers of natural disturbances prone to anthropogenic amplification. *BioScience* 58, 501–517. |
| Jönsson et al. 2012 | Jönsson, A.M. et al. (2012). Epidemics of the spruce bark beetle *Ips typographus* in Sweden. *Agricultural and Forest Management* 255, 75–82. |
| Hlásny et al. 2021 | Hlásny, T. et al. (2021). Bark beetle outbreaks in Europe: state of knowledge and ways forward. *Current Forestry Reports* 7, 138–165. |
| Kautz et al. 2011 | Kautz, M. et al. (2011). Simulating bark beetle outbreaks using high-resolution empirical data. *Forest Ecology and Management* 261, 767–777. |
| Wichmann & Ravn 2001 | Wichmann, L. & Ravn, H.P. (2001). The spread of *Ips typographus* in Denmark following wind-throw. *Agricultural and Forest Management* 148, 31–39. |
| Gardiner et al. 2013 | Gardiner, B. et al. (2013). Wind damage to forests and the influence of silviculture. *Forestry* 86, 1–5. |
| Gregow et al. 2011 | Gregow, H. et al. (2011). Combined effects of wind and snow loading on tree damage. *Silva Fennica* 45, 35–51. |
| Schelhaas et al. 2003 | Schelhaas, M.J. et al. (2003). Natural disturbances in the European forests in the 19th and 20th centuries. *Global Change Biology* 9, 1620–1633. |
| Wiegand & Moloney 2014 | Wiegand, T. & Moloney, K.A. (2014). *Handbook of Spatial Point-Pattern Analysis in Ecology*. CRC Press. |
| Reed & McKelvey 2002 | Reed, W.J. & McKelvey, K.S. (2002). Power-law behaviour and parametric models for the size-distribution of forest fires. *Ecological Modelling* 150, 239–254. |
| Clauset et al. 2009 | Clauset, A. et al. (2009). Power-law distributions in empirical data. *SIAM Review* 51, 661–703. |
| Urban & Keitt 2001 | Urban, D. & Keitt, T. (2001). Landscape connectivity: a graph-theoretic perspective. *Ecology* 82, 1205–1218. |
| Anselin 1995 | Anselin, L. (1995). Local indicators of spatial association — LISA. *Geographical Analysis* 27, 93–115. |
| Lundberg & Lee 2017 | Lundberg, S.M. & Lee, S.-I. (2017). A unified approach to interpreting model predictions. *NeurIPS* 30. |

---

## Future directions / open questions

- **Multi-year data**: `alku_aika` in 2023 only; acquiring 2021–2022 imagery from the
  same pipeline would enable true inter-annual wave-front analysis (Senf et al. 2017).
- **Agent attribution validation**: ground-truth a sample of clusters (field survey or
  existing disturbance maps from Metsäkeskus) to convert the discriminant model from
  unsupervised typology to supervised agent classification.
- **Environmental covariates**: overlay with topography (DEM), forest stand age, soil
  moisture index, and storm track data to explain spatial variation in Type0/Type1 ratios.
- **PCF interpretation**: if g(r) shows the characteristic bark-beetle double-peak
  (0–50m aggregation, 100–500m suppression), cite this against Wiegand & Moloney 2014
  Ch. 9 (marked point processes in forestry).
- **Betweenness bridges as monitoring priority**: the top-1% betweenness clusters from E9
  could be cross-referenced with Metsäkeskus monitoring plots to test whether they
  predict future outbreak expansion.
- **`pycircstat` note**: requires the deprecated `nose` package as a transitive dependency.
  If `nose` is unavailable, implement Rayleigh test directly:
  `z = n * r_bar**2; p = np.exp(-z)` (valid for n > 10).
