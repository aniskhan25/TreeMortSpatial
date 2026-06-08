# TreeMortSpatial

Reproducible analysis code for:

> Anis Ur Rahman (2026). *Ecological Drivers of Tree Mortality: A Spatial Analysis of Cluster Morphology in Finland's Forests.*

This repo is the paper's reproducibility artifact. It takes the output of the
[TreeClusters](https://github.com/aniskhan25/TreeClusters) data pipeline and
produces every figure and table in the paper.

## Quick start

```bash
pip install -r requirements.txt
# Copy clusters_with_distance.csv from TreeClusters into data/ (see data/README.md)
jupyter lab notebooks/
```

Run notebooks **in order** (01 → 05). Each notebook is self-contained for its paper section.

## Notebook → paper section map

| Notebook | Paper section | Outputs |
|---|---|---|
| `01_clustering.ipynb` | §3.2 Clustering | `filtered_clusters.csv`, Fig 1 |
| `02_shape_metrics.ipynb` | §3.3 + §4.1 Metrics | `clusters_with_metrics.csv`, Figs 2–4 |
| `03_typology.ipynb` | §4.2 Typology | `clusters_typed.csv`, Table 1, Figs 5–6 |
| `04_spatial_patterns.ipynb` | §4.3 Spatial patterns | Figs 7a–c, 8 |
| `05_advanced_analyses.ipynb` | §4.4 Advanced analyses | Figs 9–11 |

## Source modules

| Module | Purpose |
|---|---|
| `src/clustering.py` | DBSCAN wrapper, NND, cluster filtering |
| `src/metrics/shape.py` | Fractal dimension, alpha compactness, elongation, compactness |
| `src/metrics/spatial.py` | Clark-Evans index, NND stats, Ripley's K, MCF, Moran's I |
| `src/metrics/topology.py` | Euler number (Delaunay V-E+F) |
| `src/severity.py` | Composite Severity Index (CSI) |

## Data

See `data/README.md` for how to produce the input CSV from the TreeClusters pipeline.

## Dependencies

```
numpy pandas scipy matplotlib seaborn tqdm
geopandas shapely alphashape
scikit-learn
libpysal esda pointpats
joblib
```
