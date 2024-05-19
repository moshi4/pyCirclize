from __future__ import annotations

import matplotlib as mpl
import numpy as np
from Bio.SeqFeature import SeqFeature
from matplotlib.colors import Colormap, to_hex


class ColorCycler:
    """Color Cycler Class"""

    counter = 0
    cmap: Colormap = mpl.colormaps["tab10"]  # type: ignore

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
        cls.cmap = mpl.colormaps[name]  # type: ignore
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
        return to_hex(cls.cmap(n % cls.cmap.N), keep_alpha=True)  # type: ignore

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
            cmap_idx_list = list(range(0, cls.cmap.N))  # type: ignore
        elif n > 0:
            cmap_idx_list = [int(i) for i in np.linspace(0, cls.cmap.N, n)]  # type: ignore
        else:
            raise ValueError(f"{n=} is invalid number (Must be 'n > 0').")

        return [to_hex(cls.cmap(i), keep_alpha=True) for i in cmap_idx_list]  # type: ignore


def calc_group_spaces(
    groups: list[int],
    *,
    space_bw_group: float = 15,
    space_in_group: float = 2,
    endspace: bool = True,
) -> list[float]:
    """Calculate spaces between/within groups

    This function can be used to easily calculate the space size
    when data is separated into multiple groups for display.
    For example, to set up a space to divide `[A, B, C, D, E, F, G, H, I, J]`
    into three groups such as `[(A, B, C, D), (E, F, G), (H, I, J)]`,
    set `groups=[4, 3, 3]`.

    Parameters
    ----------
    groups : list[int]
        List of each group number (e.g. `[4, 3, 3]`)
    space_bw_group : float, optional
        Space size between group
    space_in_group : float, optional
        Space size within group
    endspace : bool, optional
        If True, insert space after the end group

    Returns
    -------
    spaces : list[float]
        Spaces between/within groups
    """
    if len(groups) == 0:
        raise ValueError(f"{len(groups)=} is invalid.")
    elif len(groups) == 1:
        group_num = groups[0]
        spaces = [space_in_group] * group_num
    else:
        spaces: list[float] = []
        for group_num in groups:
            group_spaces = [space_in_group] * (group_num - 1)
            group_spaces.extend([space_bw_group])
            spaces.extend(group_spaces)
    if endspace:
        return spaces
    else:
        return spaces[:-1]


def is_pseudo_feature(feature: SeqFeature) -> bool:
    """Check target feature is pseudo or not from qualifiers tag

    Parameters
    ----------
    feature : SeqFeature
        Target feature

    Returns
    -------
    check_result : bool
        pseudo check result
    """
    quals = feature.qualifiers
    return True if "pseudo" in quals or "pseudogene" in quals else False
