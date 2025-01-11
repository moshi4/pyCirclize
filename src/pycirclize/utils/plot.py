from __future__ import annotations

import math
from typing import Any


def is_lower_loc(deg: float) -> bool:
    """Check target deg is lower location or not

    Parameters
    ----------
    deg : float
        Target degree

    Returns
    -------
    result : bool
        Lower location or not
    """
    return -270 <= deg < -90 or 90 <= deg < 270


def is_right_loc(deg: float) -> bool:
    """Check target deg is right location or not

    Parameters
    ----------
    deg : float
        Target degree

    Returns
    -------
    result : bool
        Right location or not
    """
    return -360 <= deg < -180 or 0 <= deg < 180


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
    # Get parameters
    if orientation == "horizontal":
        rotation = 180 - deg if is_lower_loc(deg) else -deg
        ha = "center"
        if outer:
            va = "top" if is_lower_loc(deg) else "bottom"
        else:
            va = "bottom" if is_lower_loc(deg) else "top"
    elif orientation == "vertical":
        rotation = 90 - deg if is_right_loc(deg) else 270 - deg
        va = "center_baseline"
        if outer:
            ha = "left" if is_right_loc(deg) else "right"
        else:
            ha = "right" if is_right_loc(deg) else "left"
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
