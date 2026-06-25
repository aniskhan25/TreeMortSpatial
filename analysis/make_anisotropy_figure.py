"""
fig_anisotropy.png (house style): the G3 negative result — once matched on hull aspect
ratio, the two gradient endpoints have indistinguishable second-order linearity, so the
elongated endpoint carries no directional (wind-swath) signal beyond hull shape.

(a) point-scatter linearity distribution by endpoint (overall — confounded by hull AR);
(b) AR-matched: median linearity vs hull-AR bin, Type0 vs Type1 (the decisive panel).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from figstyle import apply_style, PALETTE, despine, panel_label, save

apply_style()
P = PALETTE
TC = {"Type0": P["secondary"], "Type1": P["primary"]}
TL = {"Type0": "Type 0 (elongated/dispersed)", "Type1": "Type 1 (compact/aggregated)"}

an = pd.read_csv("analysis/out/gate3_anisotropy.csv")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.4), dpi=150)

# (a) overall linearity distributions (confounded by hull AR)
a = ax[0]; despine(a)
bins = np.linspace(an.linearity.min(), 1.0, 40)
for t in ("Type1", "Type0"):
    a.hist(an[an.morph_type == t].linearity, bins=bins, color=TC[t], alpha=0.6,
           edgecolor="white", linewidth=0.3, label=TL[t], zorder=3, density=True)
a.set_xlabel("point-scatter linearity  $\\lambda_1/(\\lambda_1+\\lambda_2)$")
a.set_ylabel("density")
a.legend(loc="upper left")
panel_label(a, "a")
a.set_title("Overall (confounded by hull AR)", fontweight="bold", color="#222", pad=16)

# (b) AR-matched medians
b = ax[1]; despine(b)
edges = [1, 1.5, 2, 2.5, 3, 4, 10]
an["arbin"] = pd.cut(an.hull_ar, edges)
cx, m0, m1, ps = [], [], [], []
for bb, g in an.groupby("arbin", observed=True):
    g0, g1 = g[g.morph_type == "Type0"].linearity, g[g.morph_type == "Type1"].linearity
    if len(g0) >= 20 and len(g1) >= 20:
        cx.append(bb.mid); m0.append(g0.median()); m1.append(g1.median())
        ps.append(mannwhitneyu(g0, g1)[1])
b.plot(cx, m0, "^-", color=TC["Type0"], ms=8, lw=1.8, mfc="white",
       mec=TC["Type0"], mew=1.6, label=TL["Type0"], zorder=3)
b.plot(cx, m1, "o-", color=TC["Type1"], ms=7, lw=1.8, mfc="white",
       mec=TC["Type1"], mew=1.6, label=TL["Type1"], zorder=3)
for x, y0, y1, p in zip(cx, m0, m1, ps):
    if p >= 0.05:
        b.annotate("n.s.", (x, max(y0, y1)), textcoords="offset points", xytext=(0, 6),
                   ha="center", fontsize=7.5, color="#777")
b.set_xlabel("hull aspect-ratio bin (matched)")
b.set_ylabel("median point-scatter linearity")
b.legend(loc="lower right")
panel_label(b, "b")
b.set_title("AR-matched: endpoints are indistinguishable", fontweight="bold", color="#222", pad=16)

fig.tight_layout(w_pad=2.4)
save(fig, "figures/fig_anisotropy")
print("saved figures/fig_anisotropy.png + .pdf")
