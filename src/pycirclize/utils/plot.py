from __future__ import annotations

import math
from typing import Any, Literal

from matplotlib.colors import Normalize
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Bbox


def degrees(rad: float) -> float:
    """Convert radian to positive degree (`0 - 360`)

    Parameters
    ----------
    rad : float
        Target radian

    Returns
    -------
    deg : float
        Positive degree (`0 - 360`)
    """
    # Radian to degree
    deg = math.degrees(rad)
    # Normalize degree in 0 - 360 range
    deg = deg % 360
    # Negative to positive
    if deg < 0:
        deg += 360
    return deg


def is_lower_loc(rad: float) -> bool:
    """Check target radian is lower location or not

    Parameters
    ----------
    rad : float
        Target radian

    Returns
    -------
    result : bool
        Lower location or not
    """
    deg = math.degrees(rad)
    return -270 <= deg < -90 or 90 <= deg < 270


def is_right_loc(rad: float) -> bool:
    """Check target radian is right location or not

    Parameters
    ----------
    rad : float
        Target radian

    Returns
    -------
    result : bool
        Right location or not
    """
    deg = math.degrees(rad)
    return -360 <= deg < -180 or 0 <= deg < 180


def is_ann_rad_shift_target_loc(rad: float) -> bool:
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


def get_loc(
    rad: float,
) -> Literal["upper-right", "lower-right", "lower-left", "upper-left"]:
    """Get location of 4 sections

    Returns
    -------
    loc : str
        Location (`upper-right`|`lower-right`|`lower-left`|`upper-left`)
    """
    deg = degrees(rad)
    if 0 <= deg < 90:
        return "upper-right"
    elif 90 <= deg < 180:
        return "lower-right"
    elif 180 <= deg < 270:
        return "lower-left"
    else:
        return "upper-left"


def get_ann_relpos(rad: float) -> tuple[float, float]:
    """Get relative position for annotate by radian text position

    Parameters
    ----------
    rad : float
        Radian text position

    Returns
    -------
    relpos : tuple[float, float]
        Relative position
    """
    deg = degrees(rad)
    if 0 <= deg <= 180:
        return 0.0, Normalize(0, 180)(deg)
    else:
        return 1.0, 1.0 - Normalize(180, 360)(deg)


def plot_bbox(bbox: Bbox, ax: PolarAxes, **kwargs) -> None:
    """Plot bbox to check bbox area for development

    Parameters
    ----------
    bbox : Bbox
        Bounding box
    ax : PolarAxes
        Polar axes
    **kwargs : dict, optional
        Axes.plot properties (e.g. `color="red", lw=0.5, ls="--", ...`)
        <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html>
    """
    trans_bbox = bbox.transformed(ax.transAxes.inverted())
    kwargs.setdefault("clip_on", False)
    x0, y0, x1, y1 = trans_bbox.x0, trans_bbox.y0, trans_bbox.x1, trans_bbox.y1
    x, y = [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]
    ax.plot(x, y, transform=ax.transAxes, **kwargs)


def get_label_params_by_rad(
    rad: float,
    orientation: str,
    outer: bool = True,
    only_rotation: bool = False,
) -> dict[str, Any]:
    """Get proper label parameters by radian position

    Parameters
    ----------
    rad : float
        Radian coordinate
    orientation : str
        Label orientation (`horizontal` or `vertical`)
    outer : bool, optional
        If True, show on `outer` style. Else, show on `inner` style.
    only_rotation : bool, optional
        If True, Only return rotation parameter

    Returns
    -------
    dict_param : dict[str, Any]
        `va`, `ha`, `rotation`, `rotation_mode` dict
    """
    # Get position degree & location info
    deg = math.degrees(rad)
    is_lower, is_right = is_lower_loc(rad), is_right_loc(rad)
    # Get parameters
    if orientation == "horizontal":
        rotation = 180 - deg if is_lower else -deg
        ha = "center"
        if outer:
            va = "top" if is_lower else "bottom"
        else:
            va = "bottom" if is_lower else "top"
    elif orientation == "vertical":
        rotation = 90 - deg if is_right else 270 - deg
        va = "center_baseline"
        if outer:
            ha = "left" if is_right else "right"
        else:
            ha = "right" if is_right else "left"
    else:
        err_msg = f"'{orientation=} is invalid ('horizontal' or 'vertical')"
        raise ValueError(err_msg)

    if only_rotation:
        return dict(rotation=rotation)
    else:
        return dict(va=va, ha=ha, rotation=rotation, rotation_mode="anchor")


def set_axis_default_kwargs(**kwargs) -> dict[str, Any]:
    """Set axis default keyword arguments

    Set simple black axis params (`fc="none", ec="black", lw=0.5`) as default.

    Returns
    -------
    kwargs : dict[str, Any]
        Keyword arguments
    """
    if "fc" not in kwargs and "facecolor" not in kwargs:
        kwargs.update({"fc": "none"})
    if "ec" not in kwargs and "edgecolor" not in kwargs:
        kwargs.update({"ec": "black"})
    if "lw" not in kwargs and "linewidth" not in kwargs:
        kwargs.update({"lw": 0.5})
    return kwargs
