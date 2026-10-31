
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns

ROOT = Path(__file__).resolve().parent.parent
FONT_FAMILY = "Arial"

XCOL = "min_positive_S1-Ti"
YCOL = "max_SOCME_over_deltaE_S1_Ti_value"
GAP = "Gap(eV)"
SA = "SA_Score"
REPORTED = "pic16_0017"
RED, BLUE = "#cf6461", "#6799c2"
GREEN, ORANGE, GOLD = "#69b39e", "#e78a60", "#f7d447"
COLORS = [RED, BLUE]
MARKERS = ["o", "s"]
VIOLIN_PARAMETERS = {
    "width": 0.38, "bw_method": "scott", "bw_adjust": 0.96,
    "cut": 2, "gridsize": 256, "density_norm": "width",
    "xlim": [-0.6, 1.6], "box_width": 2.8, "whis_width": 0.8,
}


def style_axes(ax):
    ax.set_axisbelow(True)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.42)
    for spine in ax.spines.values():
        spine.set_color("#444444")
        spine.set_linewidth(0.9)
    ax.tick_params(length=3, width=0.65)


def legend(ax, **kwargs):
    return ax.legend(frameon=True, fancybox=False, framealpha=1,
                     edgecolor="#555555", fontsize=8, **kwargs)


def star(ax, x, y, size=150, label=None):
    return ax.scatter(x, y, marker="*", s=size, c=GOLD, edgecolors="#222222",
                      linewidths=0.85, zorder=10, label=label, clip_on=False)


def star_handle():
    return Line2D([], [], marker="*", linestyle="None", markersize=11,
                  markerfacecolor=GOLD, markeredgecolor="#222222", label="Reported")


def draw_charge(ax, groups, reported):
 
    charges = sorted(set(pd.concat(groups)["charge"]))
    bottom = np.zeros(2)
    for charge, color in zip(charges, [GREEN, ORANGE]):
        counts = np.array([(frame["charge"] == charge).sum() for frame in groups])
        ax.bar([0, 1], counts, width=0.25, bottom=bottom, color=color, alpha=0.90,
               edgecolor="#222222", linewidth=0.7, label=f"Charge +{charge:g}")
        bottom += counts
    ax.set_xticks([0, 1], ["Round 1", "Round 2"])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 63)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.set_ylabel("Count")
    ax.set_title("Charge State Distribution", weight="bold")
    style_axes(ax)
    ax.xaxis.grid(False)
    legend(ax, title="Charge State", title_fontsize=8, loc="upper right")


def draw_gap(ax, groups, reported):
    
    tidy = pd.concat([frame.assign(Round=f"Round {i+1}")
                      for i, frame in enumerate(groups)], ignore_index=True)
    sns.violinplot(data=tidy, x="Round", y=GAP, hue="Round", order=["Round 1", "Round 2"],
                   palette=COLORS, saturation=0.65, inner="box", width=VIOLIN_PARAMETERS["width"],
                   cut=VIOLIN_PARAMETERS["cut"], bw_method=VIOLIN_PARAMETERS["bw_method"],
                   bw_adjust=VIOLIN_PARAMETERS["bw_adjust"], gridsize=VIOLIN_PARAMETERS["gridsize"],
                   density_norm=VIOLIN_PARAMETERS["density_norm"], linewidth=1,
                   inner_kws={"box_width": VIOLIN_PARAMETERS["box_width"],
                              "whis_width": VIOLIN_PARAMETERS["whis_width"], "marker": "_"},
                   legend=False, ax=ax)
    for collection in ax.collections:
        collection.set_alpha(0.40)
    rng = np.random.default_rng(42)
    for i, frame in enumerate(groups):
        ax.scatter(i+rng.uniform(-0.12, 0.12, len(frame)), frame[GAP],
                   color="#555555", s=3.2, alpha=0.55, linewidths=0, zorder=4)
        star(ax, i, reported[GAP], size=125)
    handles = [Patch(facecolor=COLORS[i], alpha=0.5,
                     label=f"Round {i+1} (Median: {frame[GAP].median():.2f} eV)")
               for i, frame in enumerate(groups)] + [star_handle()]
    ax.set_xlabel("")
    ax.set_xlim(*VIOLIN_PARAMETERS["xlim"])
    ax.set_ylabel("Gap (eV)")
    ax.set_title("Distribution of HOMO-LUMO Gap", weight="bold")
    ax.set_ylim(0.5, 3.6)
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    style_axes(ax)
    ax.xaxis.grid(False)
    legend(ax, handles=handles, loc="upper center")


def draw_sa(ax, groups, reported):
    grid = np.linspace(5.0, 8.0, 700)
    for frame, color, label, bins, line_color in zip(groups, COLORS, ["Round 1", "Round 2"], [7, 8], ["#d73027", "#327eb8"]):
        ax.hist(frame[SA], bins=bins, density=True, alpha=0.50,
                color=color, edgecolor="white", linewidth=0.35, label=label)
      
        mu, sigma = stats.norm.fit(frame[SA])
        ax.plot(grid, stats.norm.pdf(grid, mu, sigma), color=line_color,
                linestyle="--", linewidth=1.8)
   
    star(ax, reported[SA], 0.10, label="Reported")
    ax.set_xlim(5, 8)
    ax.set_ylim(0, 2.95)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.set_xlabel("SA Score")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of SA Score", weight="bold")
    style_axes(ax)
    legend(ax, loc="upper right")


def draw_scatter(ax, groups, reported):
    for frame, color, marker, label in zip(groups, COLORS, MARKERS, ["Round 1", "Round 2"]):
        ax.scatter(frame[XCOL], frame[YCOL], c=color, marker=marker, s=34,
                   alpha=0.9, edgecolors="#444444", linewidths=0.4, label=label)
    star(ax, reported[XCOL], reported[YCOL], size=180, label="Reported")
    ax.set_xlim(-0.012, 0.29)
    ax.set_ylim(-5, 116)
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.set_xlabel("Min Positive S1-Ti (eV)")
    ax.set_ylabel("Max SOCME / ΔE (cm$^{-1}$/eV)")
    ax.set_title("S1-Ti Gap vs. SOCME/ΔE", weight="bold")
    style_axes(ax)
    legend(ax, loc="upper right")

def load_data(data_dir):
    filenames = ["round1_n60.csv", "round2_n60.csv"]
    groups = [pd.read_csv(data_dir / name) for name in filenames]
    required = {"ID", "charge", SA, GAP, XCOL, YCOL}
    for name, frame in zip(filenames, groups):
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{name}: missing columns {sorted(missing)}")
        if len(frame) != 60 or not frame["ID"].is_unique:
            raise ValueError(f"{name}: expected 60 unique molecules")
        if frame[list(required)].isna().any().any():
            raise ValueError(f"{name}: missing plotting data")
        if not ((frame[XCOL] < 0.3) & (frame["charge"] > 0)).all():
            raise ValueError(f"{name}: molecule outside requested filter")
    if set(groups[0]["ID"]) & set(groups[1]["ID"]):
        raise ValueError("The two groups must be disjoint")
    if groups[0][SA].mean() <= groups[1][SA].mean():
        raise ValueError("Round 1 must have the higher SA mean")
    reference = pd.concat(groups).loc[lambda frame: frame["ID"].eq(REPORTED)]
    if len(reference) != 1:
        raise ValueError("Expected exactly one pic16_0017 record across the two CSVs")
    return groups, reference.iloc[0]


def configure_fonts():
    for weight in ["normal", "bold"]:
        font_manager.findfont(font_manager.FontProperties(family=FONT_FAMILY, weight=weight),
                              fallback_to_default=False)
    plt.rcParams.update({
        "font.family": FONT_FAMILY, "font.sans-serif": [FONT_FAMILY],
        "font.size": 10, "font.weight": "normal",
        "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.labelsize": 10, "axes.labelweight": "normal",
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "legend.fontsize": 8, "legend.title_fontsize": 8,
        "mathtext.fontset": "custom", "mathtext.default": "regular",
        "mathtext.rm": FONT_FAMILY, "mathtext.it": f"{FONT_FAMILY}:italic",
        "mathtext.bf": f"{FONT_FAMILY}:bold", "mathtext.sf": FONT_FAMILY,
        "mathtext.tt": FONT_FAMILY, "mathtext.cal": FONT_FAMILY,
        "mathtext.fallback": None,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })


def harmonize_fonts(fig, ax):
    for text in fig.findobj(match=Text):
        text.set_fontfamily(FONT_FAMILY)
    ax.title.set_fontsize(11)
    ax.title.set_fontweight("bold")
    for label in [ax.xaxis.label, ax.yaxis.label]:
        label.set_fontsize(10)
        label.set_fontweight("normal")
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(9)
        label.set_fontweight("normal")
    leg = ax.get_legend()
    if leg:
        for text in leg.get_texts() + [leg.get_title()]:
            text.set_fontsize(8)
            text.set_fontweight("normal")


PLOTS = [("c_charge_distribution", draw_charge), ("d_gap_violin", draw_gap),
         ("e_sa_score_distribution", draw_sa), ("f_s1_ti_socme_scatter", draw_scatter)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "02_data")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "03_figures")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--formats", nargs="+", choices=["png", "pdf", "svg"], default=["png"])
    args = parser.parse_args()
    groups, reported = load_data(args.data_dir)
    configure_fonts()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, draw in PLOTS:
        fig, ax = plt.subplots(figsize=(6.25, 4.5))
        fig.subplots_adjust(left=0.14, right=0.97, bottom=0.15, top=0.91)
        draw(ax, groups, reported)
        harmonize_fonts(fig, ax)
        for suffix in dict.fromkeys(args.formats):
            target = args.output_dir / f"{name}.{suffix}"
            fig.savefig(target, dpi=args.dpi, facecolor="white",
                        bbox_inches="tight", pad_inches=0.12)
            print(f"Saved {target}")
        plt.close(fig)


if __name__ == "__main__":
    main()
