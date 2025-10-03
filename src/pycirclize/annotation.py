from __future__ import annotations

import warnings
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
from matplotlib.text import Annotation, Text

from pycirclize import config, utils
from pycirclize.utils.plot import degrees

if TYPE_CHECKING:
    from matplotlib.projections.polar import PolarAxes
    from matplotlib.transforms import Bbox
    from numpy.typing import NDArray


def adjust_annotation(ax: PolarAxes) -> None:
    """Adjust annotation text position"""
    # Get sorted annotation list for position adjustment
    ann_list = _get_sorted_ann_list(ax)
    if len(ann_list) == 0 or config.ann_adjust.max_iter <= 0:
        return
    if len(ann_list) > config.ann_adjust.limit:
        warnings.warn(
            f"Too many annotations(={len(ann_list)}). Annotation position adjustment is not done.",  # noqa: E501
            stacklevel=2,
        )
        return

    def get_ann_window_extent(ann: Annotation) -> Bbox:
        return Text.get_window_extent(ann).expanded(*config.ann_adjust.expand)

    # Iterate annotation position adjustment
    ax.figure.draw_without_rendering()  # type: ignore
    ann2rad_shift_candidates = _get_ann2rad_shift_candidates(ann_list)
    for idx, ann in enumerate(ann_list[1:], 1):
        orig_rad, orig_r = ann.xyann
        ann_bbox = get_ann_window_extent(ann)
        adj_ann_list = ann_list[:idx]
        adj_ann_bboxes = [get_ann_window_extent(ann) for ann in adj_ann_list]

        # Adjust radian position
        iter, max_iter = 0, config.ann_adjust.max_iter
        if _is_ann_rad_shift_target_loc(orig_rad):
            for rad_shift_candidate in ann2rad_shift_candidates[str(ann)]:
                ann.xyann = (rad_shift_candidate, orig_r)
                ann_bbox = get_ann_window_extent(ann)
                if ann_bbox.count_overlaps(adj_ann_bboxes) == 0 or iter > max_iter:
                    break
                else:
                    ann.xyann = (orig_rad, orig_r)
                iter += 1

        # Adjust radius position
        while ann_bbox.count_overlaps(adj_ann_bboxes) > 0 and iter <= max_iter:
            rad, r = ann.xyann
            ann.xyann = (rad, r + config.ann_adjust.dr)
            ann_bbox = get_ann_window_extent(ann)
            iter += 1

    # Plot annotation text bbox for developer check
    # for ann in ann_list:
    #     utils.plot.plot_bbox(get_ann_window_extent(ann), ax)


def _get_sorted_ann_list(ax: PolarAxes) -> list[Annotation]:
    """Sorted annotation list

    Sorting per 4 sections for adjusting annotation text position
    """
    ann_list = [t for t in ax.texts if isinstance(t, Annotation)]
    loc2ann_list: dict[str, list[Annotation]] = defaultdict(list)
    for ann in ann_list:
        loc = utils.plot.get_loc(ann.xyann[0])
        loc2ann_list[loc].append(ann)

    def sort_by_ann_rad(ann: Annotation) -> float:
        return utils.plot.degrees(ann.xyann[0])

    return (
        sorted(loc2ann_list["upper-right"], key=sort_by_ann_rad, reverse=True)
        + sorted(loc2ann_list["lower-right"], key=sort_by_ann_rad, reverse=False)
        + sorted(loc2ann_list["lower-left"], key=sort_by_ann_rad, reverse=True)
        + sorted(loc2ann_list["upper-left"], key=sort_by_ann_rad, reverse=False)
    )


def _is_ann_rad_shift_target_loc(rad: float) -> bool:
    """Check radian is annotation radian shift target or not

    Parameters
    ----------
    rad : float
        Annotation radian position

    Returns
    -------
    result : bool
        Target or not
    """
    deg = degrees(rad)
    return 30 <= deg <= 150 or 210 <= deg <= 330


def _get_ann2rad_shift_candidates(
    ann_list: list[Annotation],
) -> dict[str, NDArray[np.float64]]:
    """Get candidate radian shift position of annotation text

    Get the candidate radian position to shift of the target annotation
    based on the radian positions of the previous and next annotations and
    the maximum radian shift value.

    Parameters
    ----------
    ann_list : list[Annotation]
        Annotation list

    Returns
    -------
    ann2shift_rad_candidates : dict[str, NDArray[np.float64]]
        Annotation & candidate radian shift position dict
    """
    ann_list = sorted(ann_list, key=lambda a: utils.plot.degrees(a.xyann[0]))
    ann2rad_shift_candidates: dict[str, NDArray[np.float64]] = {}
    for idx, curr_ann in enumerate(ann_list):
        # Get current, prev, next annotation info
        curr_ann_rad = curr_ann.xyann[0]
        prev_ann = curr_ann if idx == 0 else ann_list[idx - 1]
        next_ann = curr_ann if idx == len(ann_list) - 1 else ann_list[idx + 1]
        prev_ann_rad, next_ann_rad = prev_ann.xyann[0], next_ann.xyann[0]
        # Get min-max radian shift position
        if abs(curr_ann_rad - prev_ann_rad) > config.ann_adjust.max_rad_shift:
            min_rad_shift = curr_ann_rad - config.ann_adjust.max_rad_shift
        else:
            min_rad_shift = prev_ann_rad
        if abs(next_ann_rad - curr_ann_rad) > config.ann_adjust.max_rad_shift:
            max_rad_shift = curr_ann_rad + config.ann_adjust.max_rad_shift
        else:
            max_rad_shift = next_ann_rad
        # Calculate candidate radian positions between min-max radian shift position
        # Sort candidate list in order of nearest to current annotation radian
        drad = config.ann_adjust.drad
        candidates = np.arange(min_rad_shift, max_rad_shift + drad, drad)
        candidates = np.append(candidates, curr_ann_rad)
        candidates = candidates[np.argsort(np.abs(candidates - curr_ann_rad))]
        ann2rad_shift_candidates[str(curr_ann)] = candidates
    return ann2rad_shift_candidates
