from __future__ import annotations

import math
from typing import Any

import matplotlib as mpl
import numpy as np
from matplotlib.colors import Colormap, to_hex


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
    is_lower_loc = True if -270 <= deg < -90 or 90 <= deg < 270 else False
    is_right_loc = True if -360 <= deg < -180 or 0 <= deg < 180 else False
    # Get parameters
    if orientation == "horizontal":
        rotation = 180 - deg if is_lower_loc else -deg
        ha = "center"
        if outer:
            va = "top" if is_lower_loc else "bottom"
        else:
            va = "bottom" if is_lower_loc else "top"
    elif orientation == "vertical":
        rotation = 90 - deg if is_right_loc else 270 - deg
        va = "center_baseline"
        if outer:
            ha = "left" if is_right_loc else "right"
        else:
            ha = "right" if is_right_loc else "left"
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


class ColorCycler:
    """Color Cycler Class"""

    counter = 0
    cmap: Colormap = mpl.colormaps["tab10"]

    def __new__(cls, n: int | None = None) -> str:
        """Get hexcolor cyclically from cmap by counter or user specified number

        `ColorCycler()` works same as `ColorCycler.get_color()`

        Parameters
        ----------
        n : int | None, optional
            Number for color cycle. If None, counter class variable is used.

        Returns
        -------
        hexcolor : str
            Cyclic hexcolor string
        """
        return cls.get_color(n)

    @classmethod
    def reset_cycle(cls) -> None:
        """Reset cycle counter"""
        cls.counter = 0

    @classmethod
    def set_cmap(cls, name: str) -> None:
        """Set colormap (Default: `tab10`)"""
        cls.cmap = mpl.colormaps[name]
        cls.counter = 0

    @classmethod
    def get_color(cls, n: int | None = None) -> str:
        """Get hexcolor cyclically from cmap by counter or user specified number

        Parameters
        ----------
        n : int | None, optional
            Number for color cycle. If None, counter class variable is used.

        Returns
        -------
        hexcolor : str
            Cyclic hexcolor string
        """
        if n is None:
            n = cls.counter
            cls.counter += 1
        return to_hex(cls.cmap(n % cls.cmap.N), keep_alpha=True)

    @classmethod
    def get_color_list(cls, n: int | None = None) -> list[str]:
        """Get hexcolor list of colormap

        Parameters
        ----------
        n : int | None, optional
            If n is None, all(=cmap.N) hexcolors are extracted from colormap.
            If n is specified, hexcolors are extracted from n equally divided colormap.

        Returns
        -------
        hexcolor_list : list[str]
            Hexcolor list
        """
        if n is None:
            cmap_idx_list = list(range(0, cls.cmap.N))
        elif n > 0:
            cmap_idx_list = [int(i) for i in np.linspace(0, cls.cmap.N, n)]
        else:
            raise ValueError(f"{n=} is invalid number (Must be 'n > 0').")

        return [to_hex(cls.cmap(i), keep_alpha=True) for i in cmap_idx_list]
