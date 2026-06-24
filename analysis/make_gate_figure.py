"""
Regenerate fig_dbscan_sensitivity.png from the gate results so the figure matches its
revised caption: (a) G1 eps-instability (Type-0 prevalence + centroid drift vs eps),
(b) G2 AR-cap prevalence dependence, (c) G2 within-Type-0 AR tail dissolving as the cap
is released.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

g1 = json.load(open("analysis/out/gate1_eps_stability.json"))
g2 = json.load(open("analysis/out/gate2_ar_sweep.json"))

eps = sorted(int(k) for k in g1)
t0_eps = [g1[str(e)]["type0_pct"] for e in eps]
ar0_eps = [g1[str(e)]["centroids"]["Type0"]["aspect_ratio"] for e in eps]

# AR caps in display order; inf shown as a large finite tick
cap_keys = ["5", "7", "10", "15", "20", "inf"]
cap_x = [5, 7, 10, 15, 20, 25]            # 25 = display slot for "inf"
cap_lab = ["5", "7", "10", "15", "20", r"$\infty$"]
t0_cap = [g2[k]["type0_pct"] for k in cap_keys]
armax_cap = [float(g2[k]["type0_ar_max"]) for k in cap_keys]
arp95_cap = [float(g2[k]["type0_ar_p95"]) for k in cap_keys]

fig, ax = plt.subplots(1, 3, figsize=(13.5, 4))

# (a) G1: prevalence + centroid drift vs eps
c0 = "#cb181d"
ax[0].plot(eps, t0_eps, "o-", color=c0, label="Type-0 prevalence (%)")
ax[0].set_xlabel("DBSCAN eps (m)")
ax[0].set_ylabel("Type-0 prevalence (%)", color=c0)
ax[0].tick_params(axis="y", labelcolor=c0)
ax[0].axvline(20, color="#888", ls="--", lw=1)
axb = ax[0].twinx()
axb.plot(eps, ar0_eps, "s--", color="#2171b5", label="Type-0 mean AR")
axb.set_ylabel("Type-0 centroid aspect ratio", color="#2171b5")
axb.tick_params(axis="y", labelcolor="#2171b5")
ax[0].set_title("(a) G1: not eps-invariant\nprevalence slides, centroid drifts")

# (b) G2: prevalence vs AR cap
ax[1].plot(cap_x, t0_cap, "o-", color=c0)
ax[1].set_xticks(cap_x); ax[1].set_xticklabels(cap_lab)
ax[1].axvline(10, color="#888", ls="--", lw=1)
ax[1].annotate("published\n23.5%", xy=(10, 23.5), xytext=(13, 29),
               fontsize=8, color="#555",
               arrowprops=dict(arrowstyle="->", color="#999"))
ax[1].set_xlabel("aspect-ratio cap")
ax[1].set_ylabel("Type-0 prevalence (%)")
ax[1].set_title("(b) G2: prevalence is an\nartefact of the AR cap")

# (c) G2: Type-0 AR tail dissolving. Grey line = the cap itself (the "wall"); the
# AR-max sitting on it at every finite cap shows the tail is pressed against the cap.
cap_val = [5, 7, 10, 15, 20, np.nan]      # inf has no wall to draw
ax[2].plot(cap_x, cap_val, ":", color="#bbb", label="cap value (the wall)")
ax[2].plot(cap_x, armax_cap, "o-", color="#6a51a3", label="within-Type-0 AR max")
ax[2].plot(cap_x, arp95_cap, "s--", color="#9e9ac8", label="within-Type-0 AR 95th pct")
ax[2].set_xticks(cap_x); ax[2].set_xticklabels(cap_lab)
ax[2].set_xlabel("aspect-ratio cap")
ax[2].set_ylabel("within-Type-0 aspect ratio")
ax[2].legend(fontsize=8, frameon=False, loc="upper left")
ax[2].set_title("(c) G2: single continuous tail,\nno second mode")

plt.tight_layout()
out = "figures/fig_dbscan_sensitivity.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print("prevalence vs eps:", list(zip(eps, t0_eps)))
print("prevalence vs cap:", list(zip(cap_lab, t0_cap)))
print("Type-0 AR max vs cap:", list(zip(cap_lab, armax_cap)))
print("saved ->", out)
