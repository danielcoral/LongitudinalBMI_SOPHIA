"""Recreate the graphical abstract as a 1500 × 600 PNG."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon


plt.rcParams["font.family"] = "Liberation Sans"

ROOT = Path(__file__).resolve().parents[1]
INK = "#203a52"
BODY = "#354b60"
MUTED = "#526377"
RED = "#d35b55"
TEAL = "#008495"
GRAY = "#70808c"


def text(ax, x, y, words, size, color=BODY, weight="normal", **kwargs):
    ax.text(x, y, words, fontsize=size, color=color, weight=weight,
            ha=kwargs.pop("ha", "center"), va="center", **kwargs)


def arrow(ax, x, y):
    ax.add_patch(Polygon(
        [(x, y - 13), (x + 45, y - 13), (x + 45, y - 25),
         (x + 70, y), (x + 45, y + 25), (x + 45, y + 13),
         (x, y + 13)],
        closed=True, facecolor="#9ab3c4", edgecolor="none"))


def bullet(ax, y, color, words):
    ax.plot(1037, y, "o", color=color, markersize=11)
    text(ax, 1060, y - 2, words, 17, ha="left", linespacing=1.25)


fig, ax = plt.subplots(figsize=(15, 6), dpi=100)
fig.subplots_adjust(0, 0, 1, 1)
ax.set_xlim(0, 1500)
ax.set_ylim(600, 0)
ax.axis("off")

text(ax, 750, 66,
     "Recent BMI trajectory distinguishes metabolic profiles across life stages",
     30, INK, "bold")
ax.hlines(110, 58, 1442, color="#cbdde9", linewidth=1.2)

text(ax, 247, 162, "Within-cohort comparison", 20, INK, "bold")
text(ax, 247, 201, "Lower vs upper slope quartile\nat the same BMI",
     14, MUTED, linespacing=1.2)
for label, y, color, start in (("Lower quartile", 285, RED, 155),
                               ("Upper quartile", 375, TEAL, 145)):
    text(ax, 66, y, label, 15, color, "bold", ha="left")
    ax.plot([start, 349], [y, 330], color=color, linewidth=4,
            solid_capstyle="butt")
ax.plot(350, 330, "o", color=INK, markersize=15)
text(ax, 346, 378, "Same BMI", 14, MUTED)

arrow(ax, 450, 311)
text(ax, 730, 165, "Three complementary life stages", 21, INK, "bold")
text(ax, 730, 270,
     "ALSPAC: adolescent/young adult\nUK Biobank: midlife\n"
     "Rotterdam: older adult",
     17, linespacing=1.16)
box = FancyBboxPatch(
    (580, 343), 300, 82,
    boxstyle="round,pad=0,rounding_size=16",
    facecolor="#f1fafc", edgecolor="#78b8c2", linewidth=1.5
)
ax.add_patch(box)
text(ax, 730, 384, "17 shared traits\nat cohort-specific\nBMI percentiles",
     14, "#175e6d", "bold", linespacing=1.1)

arrow(ax, 936, 311)
text(ax, 1025, 166, "Upper-quartile slopes aligned with",
     16, INK, "bold", ha="left", linespacing=1.15)
bullet(ax, 232, RED,
       "Higher ALT, triglycerides, LDL\nand systolic blood pressure")
bullet(ax, 313, RED,
       "Higher liver-fat ratios at higher\nBMI percentiles in both sexes")

ax.hlines(494, 58, 1442, color="#cbdde9", linewidth=1.2)
text(ax, 750, 532,
     "Stable-reference asymmetry agreed in only 64 of 102 UK Biobank–Rotterdam comparisons",
     16, MUTED)
text(ax, 750, 568,
     "Repeated routine weights add practical longitudinal context to metabolic assessment",
     17, "#086b75", "bold")

output = ROOT / "plots" / "graphical_abstract.png"
fig.savefig(output, dpi=100, facecolor="white")
plt.close(fig)
