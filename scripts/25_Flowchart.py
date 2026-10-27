"""Draw the study flowchart as a high-resolution PNG file.

Run with the Python executable in the project's RLang conda environment.
The axes use inches as coordinates: x increases rightward and y upward.
Box positions and arrow routes are fixed so the layout is reproducible.
Edit BOXES for wording, and the constants below for layout or colors.
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib

# Use a noninteractive backend so the script also works on compute nodes.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as MplPath


# Canvas and box positions are in inches. Each row has one shared height,
# which keeps boxes and vertical arrows aligned across the three columns.
WIDTH = 13.9
HEIGHT = 7.94
BOX_WIDTH = 3.36
COL_X = (0.47, 5.28, 10.09)  # Left edges: left, middle, right.
ROW_TOP = (7.46, 5.60, 3.57)  # Top edges: upper, middle, lower.
ROW_HEIGHT = (1.15, 1.33, 2.83)

# Each long arrow uses one vertical lane in the gap between columns.
# Its last segment enters the left edge of the destination box.
CD_GUTTER_X = 4.38
FG_GUTTER_X = 9.19

# Each box color is (fill, border). Text and arrows use neutral colors.
BLUE = ("#EAF2FB", "#3973AC")
GREEN = ("#EAF7EF", "#33855A")
PINK = ("#F9EDF5", "#9B4A79")
INK = "#233246"
ARROW = "#4A5969"


@dataclass(frozen=True)
class Box:
    """Content, style, and grid address of one flowchart box.

    ``col`` and ``row`` start at zero. ``bold_lines`` contains zero-based
    positions in ``lines``; an empty string adds vertical space.
    """

    name: str
    col: int
    row: int
    lines: tuple[str, ...]
    colors: tuple[str, str]
    size: float = 11.0
    align: str = "center"
    bold_lines: tuple[int, ...] = ()

    @property
    def left(self):
        return COL_X[self.col]

    @property
    def right(self):
        return self.left + BOX_WIDTH

    @property
    def top(self):
        return ROW_TOP[self.row]

    @property
    def bottom(self):
        return self.top - ROW_HEIGHT[self.row]

    @property
    def center_x(self):
        return self.left + BOX_WIDTH / 2

    @property
    def center_y(self):
        return (self.top + self.bottom) / 2


# The box IDs are for connections only; the visible labels are in ``lines``.
# Counts A/B/C were checked on 2026-10-20 against data/UKB/results/
# descriptext.txt (script 03) and recounted from the raw extracts: A is all
# rows in Weight_Recruitment_participant.csv; B is distinct IDs with any GP
# weight record (22A..) in OnlyWeights_gp.tsv, before validity filters.
BOXES = (
    # Left column: recruitment, longitudinal data, analytic sample.
    Box(
        "A", 0, 0,
        (
            "UK Biobank participants at recruitment",
            "(N = 501,962)",
        ),
        BLUE,
        bold_lines=(0,),
    ),
    Box(
        "B", 0, 1,
        (
            "Participants with longitudinal BMI data",
            "from clinical records",
            "(N = 212,991)",
        ),
        BLUE,
        bold_lines=(0,),
    ),
    Box(
        "C", 0, 2,
        (
            "Inclusion criteria", "",
            # Criteria follow scripts/02_DataPrep.R (wrec, wgpc, bmilongrec).
            "• Age at recruitment: 40–70 years",
            "• Height at recruitment ≥150 cm",
            "• Weight 50–200 kg (all measurements)",
            "• BMI at recruitment: 20–50 kg/m²",
            "• ≥1 valid BMI measurement 6 months to",
            "  5 years before recruitment",
            "• Earlier BMI: 20–50 kg/m² and 50%–150%",
            "  of recruitment BMI",
            "",
            "Final analytic sample:",
            "(N = 124,575)",
        ),
        GREEN,
        size=10.0,
        align="left",
        bold_lines=(0, 11, 12),
    ),
    # Middle column: trajectory models and research questions.
    Box(
        "D", 1, 0,
        (
            "Estimate individual trajectories",
            "using sex-specific linear mixed models",
        ),
        PINK,
        bold_lines=(0,),
    ),
    Box(
        "E", 1, 1,
        (
            "Flexible association models between",
            "trajectories and risk biomarkers",
            "at recruitment",
        ),
        PINK,
        bold_lines=(0,),
    ),
    Box(
        "F", 1, 2,
        (
            "Q1. Does prior BMI trajectory explain",
            "significant additional variation in risk",
            "biomarkers beyond current BMI?", "",
            "Q2. Does prior BMI trajectory modify the",
            "relationship between current BMI and",
            "risk biomarkers?",
        ),
        PINK,
        size=9.5,
        align="left",
        bold_lines=(0, 4),
    ),
    # Right column: inference and validation.
    Box(
        "G", 2, 0,
        (
            "Derive plausible BMI change categories",
            "across BMI categories at recruitment",
            "for inference",
        ),
        PINK,
        size=11.0,
        bold_lines=(),
    ),
    Box(
        "H", 2, 1,
        (
            "Validation",
            "Directional consistency of difference in",
            "biomarkers at recruitment between rapid",
            "gain vs loss across BMI categories",
        ),
        PINK,
        size=11.0,
        bold_lines=(0,),
    ),
    Box(
        "I", 2, 2,
        (
            "Internal:",
            "- Bootstrap (B = 1000, 80% of original data)",
            "- Sensitivity analyses:",
            "  Stricter follow-up inclusion criteria and",
            "  exclusion of weight-gain/wasting conditions",
            "",
            "External:",
            "- Rotterdam study (N = 8,087);",
            "- ALSPAC (N = 5,051)",
        ),
        PINK,
        size=10.0,
        align="left",
        bold_lines=(0, 6),
    ),
)


def draw_box(ax, box):
    """Draw one rounded rectangle and center its text block inside it."""

    patch = FancyBboxPatch(
        (box.left, box.bottom),
        BOX_WIDTH,
        ROW_HEIGHT[box.row],
        boxstyle="round,pad=0,rounding_size=0.08",
        facecolor=box.colors[0],
        edgecolor=box.colors[1],
        linewidth=2.1,
        zorder=2,
    )
    ax.add_patch(patch)

    # Font sizes are points; 72 points equal one inch. The multiplier adds
    # leading, while empty lines are shorter than lines containing text.
    line_step = box.size / 72 * 1.47
    heights = [line_step * (0.57 if not line else 1)
               for line in box.lines]
    # Start above the center, then advance downward one line at a time.
    y = box.center_y + sum(heights) / 2
    x = box.center_x if box.align == "center" else box.left + 0.2

    for index, (line, height) in enumerate(zip(box.lines, heights)):
        y -= height / 2
        if line:
            ax.text(
                x, y, line,
                ha=box.align,
                va="center",
                fontsize=box.size,
                fontweight="semibold" if index in box.bold_lines
                else "regular",
                color=INK,
                zorder=3,
            )
        y -= height / 2


def draw_arrow(ax, points):
    """Join the supplied coordinates with straight segments and an arrow."""

    # A Path preserves each right-angle turn in the long connectors.
    vertices = list(points)
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (
        len(vertices) - 1
    )
    path = MplPath(vertices, codes)
    arrow = FancyArrowPatch(
        path=path,
        arrowstyle="-|>",
        mutation_scale=20,
        linewidth=2.15,
        color=ARROW,
        capstyle="round",
        joinstyle="round",
        shrinkA=0,
        shrinkB=0,
        zorder=1,
    )
    ax.add_patch(arrow)


def make_figure():
    """Assemble all boxes and arrows on a borderless, inch-scaled canvas."""

    fig, ax = plt.subplots(figsize=(WIDTH, HEIGHT), facecolor="white")
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(0, HEIGHT)
    ax.set_aspect("equal")
    ax.axis("off")

    # Draw the six short arrows between adjacent boxes in each column.
    boxes = {box.name: box for box in BOXES}
    for upper, lower in (
        ("A", "B"), ("B", "C"),
        ("D", "E"), ("E", "F"),
        ("G", "H"), ("H", "I"),
    ):
        first = boxes[upper]
        second = boxes[lower]
        draw_arrow(
            ax,
            ((first.center_x, first.bottom),
             (second.center_x, second.top)),
        )

    # C -> D: right from C, up the first gutter, right into D's left side.
    draw_arrow(
        ax,
        ((boxes["C"].right, boxes["C"].center_y),
         (CD_GUTTER_X, boxes["C"].center_y),
         (CD_GUTTER_X, boxes["D"].center_y),
         (boxes["D"].left, boxes["D"].center_y)),
    )
    # F -> G follows the same route in the second gutter.
    draw_arrow(
        ax,
        ((boxes["F"].right, boxes["F"].center_y),
         (FG_GUTTER_X, boxes["F"].center_y),
         (FG_GUTTER_X, boxes["G"].center_y),
         (boxes["G"].left, boxes["G"].center_y)),
    )

    # Boxes sit above the lines so arrow endpoints stop at their borders.
    for box in BOXES:
        draw_box(ax, box)

    # The outer margins are already included in the axis coordinates.
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    return fig


def main():
    """Write the finished figure to the plots directory as a PNG."""

    output_dir = Path(__file__).resolve().parents[1] / "plots"
    stem = output_dir / "study_flowchart"
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
    })
    fig = make_figure()
    fig.savefig(stem.with_suffix(".png"), dpi=200,
                facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
