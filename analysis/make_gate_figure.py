"""
fig_dbscan_sensitivity.png — stability gates G1/G2, house style (see figstyle.py).
(a) G1 eps-instability: Type-0 prevalence + centroid drift vs eps.
(b) G2 prevalence dependence on the AR cap.
(c) G2 within-Type-0 AR tail dissolving as the cap is released (single mode).
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from figstyle import apply_style, PALETTE, despine, panel_label, save

apply_style()
P = PALETTE

g1 = json.load(open("analysis/out/gate1_eps_stability.json"))
g2 = json.load(open("analysis/out/gate2_ar_sweep.json"))

eps = sorted(int(k) for k in g1)
t0_eps = [g1[str(e)]["type0_pct"] for e in eps]
ar0_eps = [g1[str(e)]["centroids"]["Type0"]["aspect_ratio"] for e in eps]

cap_keys = ["5", "7", "10", "15", "20", "inf"]
cap_x = [5, 7, 10, 15, 20, 25]
cap_lab = ["5", "7", "10", "15", "20", r"$\infty$"]
t0_cap = [g2[k]["type0_pct"] for k in cap_keys]
armax = [float(g2[k]["type0_ar_max"]) for k in cap_keys]
arp95 = [float(g2[k]["type0_ar_p95"]) for k in cap_keys]

fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), dpi=150)

# (a) G1 dual axis
a = ax[0]
despine(a)
ln1 = a.plot(eps, t0_eps, "o-", color=P["primary"], lw=1.8, ms=6,
             mfc="white", mec=P["primary"], mew=1.4, label="Type-0 prevalence", zorder=3)
a.axvline(20, color=P["neutral"], ls=(0, (4, 2)), lw=1, zorder=2)
a.set_xlabel("DBSCAN eps (m)")
a.set_ylabel("Type-0 prevalence (%)", color=P["primary"])
a.tick_params(axis="y", colors=P["primary"])
a.spines["left"].set_color(P["primary"])
ab = a.twinx()
ab.spines["top"].set_visible(False)
ab.tick_params(length=3, width=0.8)
ln2 = ab.plot(eps, ar0_eps, "s--", color=P["secondary"], lw=1.6, ms=5,
              mfc="white", mec=P["secondary"], mew=1.4, label="Type-0 centroid AR", zorder=3)
ab.set_ylabel("Type-0 centroid aspect ratio", color=P["secondary"])
ab.tick_params(axis="y", colors=P["secondary"])
ab.spines["right"].set_color(P["secondary"])
a.legend(ln1 + ln2, [l.get_label() for l in ln1 + ln2], loc="upper center")
panel_label(a, "a")
a.set_title("Not eps-invariant", fontweight="bold", color="#222222", pad=18)

# (b) G2 prevalence vs cap
b = ax[1]
despine(b)
b.plot(cap_x, t0_cap, "o-", color=P["tertiary"], lw=1.8, ms=6,
       mfc="white", mec=P["tertiary"], mew=1.4, zorder=3)
b.axvline(10, color=P["neutral"], ls=(0, (4, 2)), lw=1, zorder=2)
b.annotate("published\n23.5%", xy=(10, 23.5), xytext=(13.5, 29),
           fontsize=8.5, color=P["neutral"], ha="left",
           arrowprops=dict(arrowstyle="->", color="#888", lw=1))
b.set_xticks(cap_x); b.set_xticklabels(cap_lab)
b.set_xlabel("aspect-ratio cap")
b.set_ylabel("Type-0 prevalence (%)")
panel_label(b, "b")
b.set_title("Prevalence is an artefact of the AR cap", fontweight="bold", color="#222222", pad=18)

# (c) G2 AR tail
c = ax[2]
despine(c)
c.plot(cap_x, [5, 7, 10, 15, 20, np.nan], ":", color="#B9C2C8", lw=1.4,
       label="cap value (the wall)", zorder=2)
c.plot(cap_x, armax, "o-", color=P["secondary"], lw=1.8, ms=6,
       mfc="white", mec=P["secondary"], mew=1.4, label="within-Type-0 AR max", zorder=3)
c.plot(cap_x, arp95, "s--", color=P["primary"], lw=1.6, ms=5,
       mfc="white", mec=P["primary"], mew=1.4, label="within-Type-0 AR 95th pct", zorder=3)
c.set_xticks(cap_x); c.set_xticklabels(cap_lab)
c.set_xlabel("aspect-ratio cap")
c.set_ylabel("within-Type-0 aspect ratio")
c.legend(loc="upper left")
panel_label(c, "c")
c.set_title("A single continuous tail, no second mode", fontweight="bold", color="#222222", pad=18)

fig.tight_layout(w_pad=2.4)
save(fig, "figures/fig_dbscan_sensitivity")
print("saved figures/fig_dbscan_sensitivity.png + .pdf")
