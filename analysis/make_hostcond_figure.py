"""
fig_host_conditioned.png (house style): the host-conditioned residual-aggregation result.
(a) L_inhom-r vs the spruce-host envelope; (b) vs the richer spruce+total+edge (ppm) envelope.
Observed above the envelope = residual aggregation beyond host. Plotted to the valid r<=1 km
(beyond which the finite 6 km tile windows destabilise the ring estimator -> spatstat job).
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from figstyle import apply_style, PALETTE, despine, panel_label, save

apply_style(); P = PALETTE
RMAX = 1000  # valid range (r <= ~window/4 for 6 km tiles)

def panel(ax, js, title, color):
    r = np.array(js.get("L_r", js.get("r")))
    Lo = np.array(js["L_obs"]); lo = np.array(js["L_lo"]); hi = np.array(js["L_hi"])
    m = (r > 0) & (r <= RMAX)
    r, Lo, lo, hi = r[m], Lo[m] - r[m], lo[m] - r[m], hi[m] - r[m]
    despine(ax, grid_axis="both")
    ax.fill_between(r, lo, hi, color=P["muted"], alpha=0.5, lw=0, label="host-conditioned null (95% env.)")
    ax.axhline(0, color="#999", lw=0.8, ls=":")
    ax.plot(r, Lo, "-o", color=color, lw=2, ms=5, mfc="white", mec=color, mew=1.4,
            label="observed mortality", zorder=4)
    ax.set_xlabel("distance $r$ (m)"); ax.set_ylabel(r"$L_{\mathrm{inhom}}(r) - r$  (m)")
    ax.set_title(title, fontweight="bold", color="#222", pad=8, fontsize=10.5)
    ax.legend(loc="upper left")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.4), dpi=150)
panel(ax[0], json.load(open("analysis/out/host_pcf_masked.json")),
      "vs spruce-host null", P["primary"])
panel(ax[1], json.load(open("analysis/out/host_pcf_ppm.json")),
      "vs spruce + total-volume + edge null", P["tertiary"])
panel_label(ax[0], "a"); panel_label(ax[1], "b")
fig.suptitle("Residual mortality aggregation beyond forest host (observed above the host-conditioned null)",
             fontsize=11, y=1.02)
fig.tight_layout()
save(fig, "figures/fig_host_conditioned")
print("saved figures/fig_host_conditioned.png + .pdf")
