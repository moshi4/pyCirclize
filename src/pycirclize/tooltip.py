from __future__ import annotations

import textwrap
import uuid
from typing import TYPE_CHECKING, Any

import numpy as np
from matplotlib.colors import is_color_like, to_rgb
from matplotlib.patches import Patch, Rectangle

from pycirclize import config
from pycirclize.utils.plot import degrees, is_lower_loc, is_right_loc, select_textcolor

if TYPE_CHECKING:
    from collections.abc import Sequence

    from Bio.Phylo.BaseTree import Clade
    from Bio.SeqFeature import SeqFeature
    from matplotlib.backend_bases import MouseEvent
    from matplotlib.collections import Collection
    from matplotlib.projections import PolarAxes

    from pycirclize.parser.bed import BedRecord


def gen_gid(prefix: str | None = None) -> str:
    """Generate Group ID for tooltip

    Generate unique id using `uuid.uuid4()` function

    Parameters
    ----------
    prefix : str | None, optional
        Group ID prefix

    Returns
    -------
    gid : str
        Group ID
    """
    gid = uuid.uuid4().hex
    return f"{prefix}-{gid}" if prefix else gid


def to_feature_tooltip(feature: SeqFeature) -> str:
    """Convert feature to tooltip text"""
    start, end = int(feature.location.start), int(feature.location.end)  # type: ignore
    strand = "-" if feature.location.strand == -1 else "+"
    type = feature.type
    gene = str(feature.qualifiers.get("gene", ["na"])[0])
    product = str(feature.qualifiers.get("product", ["na"])[0])
    product = "\n".join(textwrap.wrap(product, width=40))
    return "\n".join(
        [
            f"location: {start:,} - {end:,} ({strand})",
            f"length: {end - start:,}",
            f"type: {type}",
            f"gene: {gene}",
            f"product: {product}",
        ]
    )


def to_link_tooltip(
    sector_region1: tuple[str, float, float],
    sector_region2: tuple[str, float, float],
    direction: int,
) -> str:
    """Convert link data to tooltip text"""
    name1, start1, end1 = sector_region1
    size1 = abs(end1 - start1)
    name2, start2, end2 = sector_region2
    size2 = abs(end2 - start2)
    direction2symbol = {0: "-", 1: "-->", -1: "<--", 2: "<-->"}
    return textwrap.dedent(
        f"""
        Source
          {name1}: {size1:,}
        Target
          {name2}: {size2:,}
        Relation
          {name1} {direction2symbol[direction]} {name2}
        """
    )[1:-1]


def to_node_tooltip(node: Clade) -> str:
    """Convert tree node to tooltip text"""
    tooltip = f"Node: {node.name}\nLength: {node.branch_length}"
    if node.confidence is not None:
        tooltip += f"\nBootstrap: {node.confidence}"
    return tooltip


def to_cytoband_tooltip(rec: BedRecord) -> str:
    """Convert cytoband bed record to tooltip text"""
    return textwrap.dedent(
        f"""
        Cytoband: {rec.name} ({rec.chr})
        Start: {rec.start:,}
        End: {rec.end:,}
        """
    )[1:-1]


def set_patch_tooltip(
    ax: PolarAxes,
    patches: Sequence[Patch],
    gid2tooltip: dict[str, str],
) -> None:
    """Set patch tooltip annotation

    Parameters
    ----------
    ax : PolarAxes
        Target polar axes
    patches : Sequence[Patch]
        List of patches
    gid2tooltip : dict[str, str]
        Patch group id & tooltip dict

    Note
    ----
        To display tooltips, it is necessary to first identify the patch
        at the mouse position. Since it is a heavy load to check all
        the patches against the mouse position every time, the processing
        speed can be improved by dividing the patches into segments in advance
        by a certain radian range.
    """
    gid = gen_gid("patch-tooltip")
    # Hover target patches
    target_patches = [p for p in patches if p.get_gid()]
    seg_patches_list = _segmentize_patches(target_patches)

    def hover(e: MouseEvent) -> None:
        # Remove previous tooltip annotation if exists
        for ann in ax.texts:
            if ann.get_gid() == gid:
                ann.remove()

        # Check axes is target or not
        if e.inaxes != ax:
            return

        # Get hover patch
        hover_patch: Patch | None = None
        for patch in seg_patches_list[_get_mouse_segment_idx(e)]:
            contains, _ = patch.contains(e)
            if contains:
                hover_patch = patch
                break

        # Plot tooltip annotation for hover patch
        if hover_patch is not None:
            tooltip = gid2tooltip.get(str(hover_patch.get_gid()))
            if tooltip is None:
                return
            color = to_rgb(hover_patch.get_facecolor())
            ax.annotate(**_build_tooltip_ann_kwargs(tooltip, gid, e, color, ax))

    # Set hover event handler
    fig = ax.get_figure()
    fig.canvas.mpl_connect("motion_notify_event", hover)  # type: ignore
    for key, value in dict(
        toolbar_visible=False,
        header_visible=False,
        footer_visible=False,
        resizable=True,
    ).items():
        if hasattr(fig.canvas, key):
            setattr(fig.canvas, key, value)


SEG_NUM = 360
SEG_DEG = 1


def _segmentize_patches(patches: list[Patch]) -> list[list[Patch]]:
    """Segmentize patches by patch radian position"""
    seg_patches_list: list[list[Patch]] = [[] for _ in range(SEG_NUM)]
    for p in patches:
        rads = np.array(p.get_path().vertices)[:, 0]
        if isinstance(p, Rectangle):  # For track.bar() Rectangle Patch
            degs = [degrees(p.xy[0]), degrees(p.xy[0] + p.get_width())]
        else:
            degs = list(map(degrees, rads))
        degs = np.append(np.arange(min(degs), max(degs), 0.01), max(degs))
        indices = set(map(lambda v: int(v // SEG_DEG), degs))
        indices = filter(lambda idx: idx if idx < SEG_NUM else SEG_NUM - 1, indices)
        for idx in indices:
            seg_patches_list[idx].append(p)
    return seg_patches_list


def _get_mouse_segment_idx(e: MouseEvent) -> int:
    """Get segment index by mouse radian position"""
    deg = degrees(float(e.xdata))  # type: ignore
    return int(deg // SEG_DEG)


def set_collection_tooltip(ax: PolarAxes, coll: Collection, labels: list[str]) -> None:
    """Set collection tooltip annotation

    Parameters
    ----------
    ax : PolarAxes
        Target polar axes
    coll : Collection
        Plot collection
    labels : list[str]
        Tooltip labels
    """
    gid = gen_gid("collection-tooltip")

    def hover(e: MouseEvent) -> None:
        # Remove previous tooltip annotation if exists
        for ann in ax.texts:
            if ann.get_gid() == gid:
                ann.remove()

        # Check axes is target or not
        if e.inaxes != ax:
            return

        # Get hover target
        contains, ctx = coll.contains(e)
        label, color = None, None
        if contains:
            hover_coll_idx = ctx["ind"][0]
            label = labels[hover_coll_idx]
            color = coll.get_facecolor()
            if is_color_like(color):
                color = to_rgb(color)  # type: ignore
            else:
                color = to_rgb(color[hover_coll_idx])

        if label is not None and color is not None:
            ax.annotate(**_build_tooltip_ann_kwargs(label, gid, e, color, ax))

    # Set hover event handler
    fig = ax.get_figure()
    fig.canvas.mpl_connect("motion_notify_event", hover)  # type: ignore


def _build_tooltip_ann_kwargs(
    text: str,
    gid: str,
    e: MouseEvent,
    fc: tuple[float, float, float],
    ax: PolarAxes,
) -> dict[str, Any]:
    """Build tooltip annotation keyword arguments"""
    tooltip_fc = fc if config.tooltip.fc is None else to_rgb(config.tooltip.fc)
    tooltip_textcolor = select_textcolor(tooltip_fc)
    rad, _ = float(e.xdata), float(e.ydata)  # type: ignore
    xtext = -config.tooltip.margin if is_right_loc(rad) else config.tooltip.margin
    return dict(
        text=text,
        xy=(e.x, e.y),
        xytext=(xtext, 0),
        xycoords=ax.transScale,
        textcoords="offset pixels",
        size=config.tooltip.fontsize,
        gid=gid,
        ha="right" if is_right_loc(rad) else "left",
        va="bottom" if is_lower_loc(rad) else "top",
        ma="left",
        color=tooltip_textcolor,
        zorder=999,
        bbox=dict(
            fc=tooltip_fc,
            ec=tooltip_textcolor,
            lw=config.tooltip.lw,
            boxstyle=config.tooltip.boxstyle,
        ),
    )
