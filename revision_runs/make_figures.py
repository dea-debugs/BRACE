"""Publication figures from the revision runs -> figures/C3_km_150.png and figures/E_policy_ranking.png (300 dpi).

Style follows the paper's figures/meta_pam_surv_150.png: matplotlib defaults, Followed in red, Not Followed in black,
survival probability 0 to 1, time axis in months to 150.
"""
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter

from brace_paths import REV

OUT = REV / "figures"
OUT.mkdir(exist_ok=True)
RED, BLACK, GREY = "#d62728", "black", "0.5"
TAU = 150
TICKS = np.arange(0, TAU + 1, 30)

def fmt(x):  # typographic minus sign, not a hyphen
    return f"{x:.2f}".replace("-", "\u2212")


# ---------------- 1. C3 Kaplan Meier curves ----------------
label = "C3_pam50_clinical_seed123"
d = pd.read_csv(REV / f"{label}.csv")
ev = json.loads((REV / f"{label}_eval.json").read_text())
assert len(d) == 594
d["followed"] = d.recommendation == d.treatment

panels = [("a", "Pooled", d, "pooled"),
          ("b", "No chemotherapy arm", d[d.treatment == 0], "no_chemo_arm"),
          ("c", "Chemotherapy arm", d[d.treatment == 1], "chemo_arm")]

fig = plt.figure(figsize=(12, 4.5))
gs = fig.add_gridspec(2, 3, height_ratios=[3.6, 0.75], hspace=0.38, wspace=0.3,
                      left=0.1, right=0.985, top=0.93, bottom=0.03)
for j, (letter, name, sub, key) in enumerate(panels):
    ax = fig.add_subplot(gs[0, j])
    tab = fig.add_subplot(gs[1, j])
    groups = [("Followed", sub[sub.followed], RED), ("Not Followed", sub[~sub.followed], BLACK)]
    for i, (gname, g, col) in enumerate(groups):
        km = KaplanMeierFitter().fit(g.time, g.event, label=f"{gname} (n = {len(g)})")
        km.plot_survival_function(ax=ax, ci_show=False, color=col, lw=1.5, show_censors=True,
                                  censor_styles={"ms": 4, "marker": "|", "mew": 0.8})
        at_risk = [(g.time >= t).sum() for t in TICKS]
        for t, n in zip(TICKS, at_risk):
            tab.text(t, 1 - i, str(n), ha="center", va="center", fontsize=8.5, color=col)
    e = ev[key]
    ax.text(0.03, 0.05, f"RMST difference at 150 months\n{fmt(e['diff'])} (95% CI {fmt(e['lower'])} to {fmt(e['upper'])})\n"
            f"p = {e['p']:.3f}", transform=ax.transAxes, fontsize=8.5, va="bottom")
    ax.set_xlim(0, TAU); ax.set_ylim(0, 1.0)
    ax.set_xticks(TICKS)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability" if j == 0 else "")
    ax.legend(loc="upper right", fontsize=8.5, frameon=False)
    ax.text(-0.02, 1.03, f"({letter}) {name}", transform=ax.transAxes, fontsize=11, ha="left", va="bottom")
    tab.set_xlim(0, TAU); tab.set_ylim(-0.6, 1.6); tab.axis("off")
    if j == 0:
        for i, (gname, col) in enumerate((("Followed", RED), ("Not Followed", BLACK))):
            tab.text(-0.08, 1 - i, gname, transform=tab.get_yaxis_transform(), ha="right", va="center",
                     fontsize=8.5, color=col)
    tab.set_title("Number at risk", fontsize=8.5, loc="left", pad=2)
fig.savefig(OUT / "C3_km_150.png", dpi=300)
plt.close(fig)

# ---------------- 2. E policy ranking ----------------
R = pd.read_csv(REV / "E_policy_ranking.csv")  # keep CSV order (ranked by mean V minus observed)
ref_names = {"reference_none": "No chemo for anyone (reference)",
             "reference_all": "Chemo for everyone (reference)",
             "reference_observed": "Observed care (reference)"}


def pretty(row):
    if row.key in ref_names:
        return ref_names[row.key]
    s = row.input.replace("PAM50 clean + clinical (C3 spec)", "PAM50 clean plus clinical")
    s = s.replace("Pearson corr. coeff. (100)", "Pearson correlation (100)")
    return re.sub(r"K_(in|out)", r"K$_{\1}$", s)


fig, ax = plt.subplots(figsize=(7, 7))
n = len(R)
y = np.arange(n)[::-1]
for yi, (_, r) in zip(y, R.iterrows()):
    ref = r.key.startswith("reference")
    hi = r.key == "pam50clean_clinical"
    col = GREY if ref else (RED if hi else BLACK)
    if hi:
        ax.axhspan(yi - 0.45, yi + 0.45, color=RED, alpha=0.08, lw=0)
    ax.hlines(yi, r.V_minus_observed_min, r.V_minus_observed_max, color=col, lw=2.2 if hi else 1.4)
    ax.plot(r.V_minus_observed_mean, yi, "o", color=col, ms=6.5 if hi else 5,
            mfc="white" if ref else col, mew=1.3)
ax.axvline(0, ls="--", color="0.3", lw=1)
ax.set_yticks(y)
ax.set_yticklabels([pretty(r) for _, r in R.iterrows()], fontsize=9)
for tl, (_, r) in zip(ax.get_yticklabels(), R.iterrows()):
    if r.key.startswith("reference"):
        tl.set_color(GREY); tl.set_fontstyle("italic")
    elif r.key == "pam50clean_clinical":
        tl.set_color(RED); tl.set_fontweight("bold")
ax.set_ylim(-0.7, n - 0.3)
ax.set_xlabel("V(BRACE) minus V(observed care) (months)")
ax.tick_params(axis="y", length=0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.text(0.99, 0.005, "Dot: mean over seven splits. Line: minimum to maximum.", ha="right", va="bottom", fontsize=7.5,
         color="0.35")
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(OUT / "E_policy_ranking.png", dpi=300)
plt.close(fig)
print("wrote", OUT / "C3_km_150.png", OUT / "E_policy_ranking.png")
