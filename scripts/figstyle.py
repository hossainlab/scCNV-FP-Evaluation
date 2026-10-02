"""Shared figure style for fig1.py ... fig4.py.

The original figures were rendered in an environment that injected these helpers
into the script namespace; they were not part of the supplement. This module is a
reconstruction matched by eye to the published PNG/PDF files (Matplotlib 3.11.1,
DejaVu Sans). Data, layout and colours come from the figure scripts themselves, so
regenerated figures carry identical content; spacing and font sizes may differ
slightly from the published versions.
"""
import matplotlib as mpl
from matplotlib.transforms import Bbox

META_GREY = "#888888"


def apply_figure_style(sizes=(8, 7, 6)):
    """sizes = (axis labels, legend/body text, tick labels), in points."""
    label, body, tick = sizes
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": body,
        "axes.labelsize": label,
        "axes.titlesize": label,
        "legend.fontsize": body,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "pdf.fonttype": 42,
        "savefig.dpi": 330,
    })


def panel_letter(ax, letter, dx=None, dy=0.04, size=9):
    """Bold lower-case panel label above the top-left corner of ``ax``.

    dx is an x offset in axes fraction. By default the label sits just left of the
    outermost decoration (tick labels, y-label) of every axes sharing this axes'
    left edge, so letters in one grid column line up. Resolved at draw time, so
    axes created after this call are taken into account.
    """
    if dx is not None:
        ax.text(dx, 1.0 + dy, letter, transform=ax.transAxes, fontsize=size,
                fontweight="bold", ha="left", va="bottom")
        return
    fig = ax.figure

    def anchor(renderer):
        x0 = ax.get_position().x0
        col = [a for a in fig.axes if abs(a.get_position().x0 - x0) < 1e-3]
        boxes = [b for a in col for b in (a.get_window_extent(renderer),
                                          a.xaxis.get_tightbbox(renderer),
                                          a.yaxis.get_tightbbox(renderer)) if b is not None]
        left = min(b.x0 for b in boxes)
        ext = ax.get_window_extent(renderer)
        return Bbox.from_extents(left, ext.y0, ext.x1, ext.y1)

    ax.annotate(letter, xy=(0, 1.0 + dy), xycoords=anchor, xytext=(-5, 0),
                textcoords="offset points", fontsize=size, fontweight="bold",
                ha="left", va="bottom", annotation_clip=False)
