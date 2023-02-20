from __future__ import annotations

import itertools
import math
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from matplotlib.projections.polar import PolarAxes

from pycirclize import config, utils
from pycirclize.parser import Bed, Matrix
from pycirclize.patches import ArcLine, ArcRectangle, BezierCurve
from pycirclize.sector import Sector
from pycirclize.track import Track


class Circos:
    """Circos Visualization Class"""

    def __init__(
        self,
        sectors: dict[str, int] | dict[str, float],
        start: float = 0,
        end: float = 360,
        *,
        space: float | list[float] = 0,
        endspace: bool = True,
        sectors_start_pos: dict[str, int] | dict[str, float] = {},
        show_axis_for_debug: bool = False,
    ):
        """
        Parameters
        ----------
        sectors : dict[str, int] | dict[str, float]
            Sector name & size dict
        start : float, optional
            Plot start degree (`-360 <= start < end <= 360`)
        end : float, optional
            Plot end degree (`-360 <= start < end <= 360`)
        space : float | list[float], optional
            Space degree(s) between sector
        endspace : bool, optional
            If True, insert space after the end sector
        sectors_start_pos : dict[str, int] | dict[str, float], optional
            Sector name & start position dict. By default, `start_pos=0`.
        show_axis_for_debug : bool, optional
            Show axis for position check debugging (Developer option)
        """
        # Check start-end degree range
        self._check_degree_range(start, end)

        # Calculate sector region & add sector
        whole_deg_size = end - start
        space_num = len(sectors) if endspace else len(sectors) - 1
        if isinstance(space, (list, tuple)):
            if len(space) != space_num:
                err_msg = f"{space=} is invalid.\n"
                err_msg += f"Length of space list must be {space_num}."
                raise ValueError(err_msg)
            space_list = list(space) + [0]
            space_deg_size = sum(space)
        else:
            space_deg_size = space * space_num
            space_list = [space] * space_num + [0]
        whole_deg_size_without_space = whole_deg_size - space_deg_size
        sector_total_size = sum(sectors.values())

        rad_pos = math.radians(start)
        self._sectors: list[Sector] = []
        for idx, (sector_name, sector_size) in enumerate(sectors.items()):
            sector_size_ratio = sector_size / sector_total_size
            deg_size = whole_deg_size_without_space * sector_size_ratio
            rad_size = math.radians(deg_size)
            rad_lim = (rad_pos, rad_pos + rad_size)
            rad_pos += rad_size + math.radians(space_list[idx])
            start_pos = sectors_start_pos.get(sector_name, 0)
            sector = Sector(sector_name, sector_size, rad_lim, start_pos)
            self._sectors.append(sector)

        self._deg_lim = (start, end)
        self._rad_lim = (math.radians(start), math.radians(end))
        self._patches: list[Patch] = []
        self._plot_funcs: list[Callable[[PolarAxes], None]] = []
        self._show_axis_for_debug = show_axis_for_debug

    ############################################################
    # Property
    ############################################################

    @property
    def rad_size(self) -> float:
        """Circos radian size"""
        return max(self.rad_lim) - min(self.rad_lim)

    @property
    def rad_lim(self) -> tuple[float, float]:
        """Circos radian limit"""
        return self._rad_lim

    @property
    def deg_size(self) -> float:
        """Circos degree size"""
        return max(self.deg_lim) - min(self.deg_lim)

    @property
    def deg_lim(self) -> tuple[float, float]:
        """Circos degree limit"""
        return self._deg_lim

    @property
    def sectors(self) -> list[Sector]:
        """Sectors"""
        return self._sectors

    @property
    def tracks(self) -> list[Track]:
        """Tracks (from sectors)"""
        tracks = []
        for sector in self.sectors:
            for track in sector.tracks:
                tracks.append(track)
        return tracks

    ############################################################
    # Public Method
    ############################################################

    @staticmethod
    def initialize_from_matrix(
        matrix: str | Path | pd.DataFrame | Matrix,
        *,
        start: float = 0,
        end: float = 360,
        space: float | list[float] = 0,
        endspace: bool = True,
        r_lim: tuple[float, float] = (97, 100),
        cmap: str | dict[str, str] = "viridis",
        link_cmap: list[tuple[str, str, str]] = [],
        ticks_interval: int | None = None,
        order: str | list[str] | None = None,
        label_kws: dict[str, Any] = {},
        ticks_kws: dict[str, Any] = {},
        link_kws: dict[str, Any] = {},
    ) -> Circos:
        """Initialize Circos instance from Matrix

        Circos tracks and links are auto-defined by Matrix (For plotting Chord Diagram)

        Parameters
        ----------
        matrix : str | Path | pd.DataFrame | Matrix
            Matrix file or Matrix dataframe or Matrix instance
        start : float, optional
            Plot start degree (-360 <= start < end <= 360)
        end : float, optional
            Plot end degree (-360 <= start < end <= 360)
        space : float | list[float], optional
            Space degree(s) between sector
        endspace : bool, optional
            If True, insert space after the end sector
        r_lim : tuple[float, float], optional
            Outer track radius limit region (0 - 100)
        cmap : str | dict[str, str], optional
            Colormap assigned to each outer track and link.
            User can set matplotlib's colormap (e.g. `viridis`, `jet`, `tab10`) or
            label_name -> color dict (e.g. `dict(A="red", B="blue", C="green", ...)`)
        link_cmap : list[tuple[str, str, str]], optional
            Link colormap to overwrite link colors automatically set by cmap.
            User can set list of `from_label`, `to_label`, `color` tuple
            (e.g. `[("A", "B", "red"), ("A", "C", "#ffff00), ...]`)
        ticks_interval : int | None, optional
            Ticks interval. If None, ticks are not plotted.
        order : str | list[str] | None, optional
            Sort order of matrix for plotting Chord Diagram. If `None`, no sorting.
            If `asc`|`desc`, sort in ascending(or descending) order by node size.
            If node name list is set, sort in user specified node order.
        label_kws : dict[str, Any], optional
            Keyword arguments passed to `sector.text()` method
            (e.g. `dict(r=110, orientation="vertical", size=15, ...)`)
        ticks_kws : dict[str, Any], optional
            Keyword arguments passed to `track.xticks_by_interval()` method
            (e.g. `dict(label_size=10, label_orientation="vertical", ...)`)
        link_kws : dict[str, Any], optional
            Keyword arguments passed to `circos.link()` method
            (e.g. `dict(direction=1, ec="black", lw=0.5, alpha=0.8, ...)`)

        Returns
        -------
        circos : Circos
            Circos instance initialized from Matrix
        """
        # If input matrix is file path, convert to Matrix instance
        if isinstance(matrix, (str, Path, pd.DataFrame)):
            matrix = Matrix(matrix)

        # Sort matrix if order is set
        if order is not None:
            matrix = matrix.sort(order)

        # Get name2color dict from user-specified colormap
        names = matrix.all_names
        if isinstance(cmap, str):
            utils.ColorCycler.set_cmap(cmap)
            colors = utils.ColorCycler.get_color_list(len(names))
            name2color = dict(zip(names, colors))
        else:
            if type(cmap) == defaultdict:
                name2color = cmap
            else:
                name2color: dict[str, str] = defaultdict(lambda: "grey")
                name2color.update(cmap)

        # Initialize circos sectors
        circos = Circos(matrix.to_sectors(), start, end, space=space, endspace=endspace)
        for sector in circos.sectors:
            # Plot label, outer track axis & xticks
            sector.text(sector.name, **label_kws)
            outer_track = sector.add_track(r_lim)
            color = name2color[sector.name]
            outer_track.axis(fc=color)
            if ticks_interval is not None:
                outer_track.xticks_by_interval(ticks_interval, **ticks_kws)

        # Plot links
        fromto_label2color = {f"{t[0]}{t[1]}": t[2] for t in link_cmap}
        for link in matrix.to_links():
            from_label, to_label = link[0][0], link[1][0]
            fromto_label = f"{from_label}{to_label}"
            if fromto_label in fromto_label2color:
                color = fromto_label2color[fromto_label]
            else:
                color = name2color[from_label]
            circos.link(*link, fc=color, **link_kws)

        return circos

    @staticmethod
    def initialize_from_bed(
        bed_file: str | Path,
        start: float = 0,
        end: float = 360,
        *,
        space: float | list[float] = 0,
        endspace: bool = True,
    ) -> Circos:
        """Initialize Circos instance from BED file

        Circos sectors are auto-defined by BED chromosomes

        Parameters
        ----------
        bed_file : str | Path
            Chromosome BED format file (zero-based coordinate)
        start : float, optional
            Plot start degree (-360 <= start < end <= 360)
        end : float, optional
            Plot end degree (-360 <= start < end <= 360)
        space : float | list[float], optional
            Space degree(s) between sector
        endspace : bool, optional
            If True, insert space after the end sector

        Returns
        -------
        circos : Circos
            Circos instance initialized from BED file
        """
        records = Bed(bed_file).records
        sectors = {rec.chr: rec.size for rec in records}
        sectors_start_pos = {rec.chr: rec.start for rec in records}
        return Circos(
            sectors,
            start,
            end,
            space=space,
            endspace=endspace,
            sectors_start_pos=sectors_start_pos,
        )

    def add_cytoband_tracks(
        self,
        r_lim: tuple[float, float],
        cytoband_file: str | Path,
        *,
        track_name: str = "cytoband",
        cytoband_cmap: dict[str, str] | None = None,
    ) -> None:
        """Add track & plot chromosome cytoband on each sector

        Parameters
        ----------
        r_lim : tuple[float, float]
            Radius limit region (0 - 100)
        cytoband_file : str | Path
            Cytoband tsv file (UCSC format)
        track_name : str, optional
            Cytoband track name. By default, `cytoband`.
        cytoband_cmap : dict[str, str] | None, optional
            User-defined cytoband colormap. If None, use Circos style colormap.
            (e.g. `{"gpos100": "#000000", "gneg": "#FFFFFF", ...}`)
        """
        if cytoband_cmap is None:
            cytoband_cmap = config.CYTOBAND_COLORMAP
        cytoband_records = Bed(cytoband_file).records
        for sector in self.sectors:
            track = sector.add_track(r_lim, name=track_name)
            track.axis()
            for rec in cytoband_records:
                if sector.name == rec.chr:
                    color = cytoband_cmap.get(str(rec.score), "white")
                    track.rect(rec.start, rec.end, fc=color)

    def get_sector(self, name: str) -> Sector:
        """Get sector by name

        Parameters
        ----------
        name : str
            Sector name

        Returns
        -------
        sector : Sector
            Sector
        """
        name2sector = {s.name: s for s in self.sectors}
        if name not in name2sector:
            raise ValueError(f"{name=} sector not found.")
        return name2sector[name]

    def axis(self, **kwargs) -> None:
        """Plot axis

        By default, simple black axis params(`fc="none", ec="black", lw=0.5`) are set.

        Parameters
        ----------
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="blue", lw=0.5, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Set default params
        kwargs = utils.plot.set_axis_default_kwargs(**kwargs)

        # Axis facecolor placed behind other patches (zorder=0.99)
        fc_behind_kwargs = {**kwargs, **config.AXIS_FACE_PARAM}
        self.rect(**fc_behind_kwargs)

        # Axis edgecolor placed in front of other patches (zorder=1.01)
        ec_front_kwargs = {**kwargs, **config.AXIS_EDGE_PARAM}
        self.rect(**ec_front_kwargs)

    def text(
        self,
        text: str,
        r: float = 0,
        deg: float = 0,
        **kwargs,
    ) -> None:
        """Plot text

        Parameters
        ----------
        text : str
            Text content
        r : float
            Radius position
        deg : float
            Degree position (0 - 360)
        **kwargs : dict, optional
            Text properties (e.g. `size=12, color="red", rotation=90, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        if "va" not in kwargs and "verticalalignment" not in kwargs:
            kwargs.update(dict(va="center"))
        if "ha" not in kwargs and "horizontalalignment" not in kwargs:
            kwargs.update(dict(ha="center"))

        def plot_text(ax: PolarAxes) -> None:
            ax.text(math.radians(deg), r, text, **kwargs)

        self._plot_funcs.append(plot_text)

    def line(
        self,
        r: float,
        deg_lim: tuple[float, float] | None = None,
        **kwargs,
    ) -> None:
        """Plot line

        Parameters
        ----------
        r : float
            Radius position (0 - 100)
        deg_lim : tuple[float, float]
            Degree limit region (-360 - 360). If None, `circos.deg_lim` is set.
        **kwargs : dict, optional
            Patch properties (e.g. `color="red", lw=3, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        deg_lim = self.deg_lim if deg_lim is None else deg_lim
        rad_lim = (math.radians(min(deg_lim)), math.radians(max(deg_lim)))
        self._patches.append(ArcLine(rad_lim, (r, r), **kwargs))

    def rect(
        self,
        r_lim: tuple[float, float] = (0, 100),
        deg_lim: tuple[float, float] | None = None,
        **kwargs,
    ) -> None:
        """Plot rectangle

        Parameters
        ----------
        r_lim : tuple[float, float]
            Radius limit region (0 - 100)
        deg_lim : tuple[float, float]
            Degree limit region (-360 - 360). If None, `circos.deg_lim` is set.
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="black", lw=1, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        deg_lim = self.deg_lim if deg_lim is None else deg_lim
        rad_lim = (math.radians(min(deg_lim)), math.radians(max(deg_lim)))

        radr = (min(rad_lim), min(r_lim))
        width = max(rad_lim) - min(rad_lim)
        height = max(r_lim) - min(r_lim)
        self._patches.append(ArcRectangle(radr, width, height, **kwargs))

    def link(
        self,
        sector_region1: tuple[str, float, float],
        sector_region2: tuple[str, float, float],
        r1: float | None = None,
        r2: float | None = None,
        *,
        color: str = "grey",
        alpha: float = 0.5,
        height_ratio: float = 0.5,
        direction: int = 0,
        arrow_length_ratio: float = 0.05,
        allow_twist: bool = True,
        **kwargs,
    ) -> None:
        """Plot link to specified region within or between sectors

        Parameters
        ----------
        sector_region1 : tuple[str, float, float]
            Link sector region1 (name, start, end)
        sector_region2 : tuple[str, float, float]
            Link sector region2 (name, start, end)
        r1 : float | None, optional
            Link radius end position for sector_region1.
            If None, lowest radius position of track in target sector is set.
        r2 : float | None, optional
            Link radius end position for sector_region2.
            If None, lowest radius position of track in target sector is set.
        color : str, optional
            Link color
        alpha : float, optional
            Link color alpha (transparency) value
        height_ratio : float, optional
            Bezier curve height ratio
        direction : int, optional
            `0`: No direction edge shape (Default)
            `1`: Forward direction arrow edge shape (region1 -> region2)
            `-1`: Reverse direction arrow edge shape (region1 <- region2)
            `2`: Bidirectional arrow edge shape (region1 <-> region2)
        arrow_length_ratio : float, optional
            Direction arrow length ratio
        allow_twist : bool, optional
            If False, twisted link is automatically resolved.
            <http://circos.ca/documentation/tutorials/links/twists/images>
        **kwargs : dict, optional
            Patch properties (e.g. `ec="red", lw=1.0, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Set data for plot link
        name1, start1, end1 = sector_region1
        name2, start2, end2 = sector_region2
        sector1, sector2 = self.get_sector(name1), self.get_sector(name2)
        r1 = sector1.get_lowest_r() if r1 is None else r1
        r2 = sector2.get_lowest_r() if r2 is None else r2
        rad_start1, rad_end1 = sector1.x_to_rad(start1), sector1.x_to_rad(end1)
        rad_start2, rad_end2 = sector2.x_to_rad(start2), sector2.x_to_rad(end2)

        # Set patch kwargs & default linewidth as 0.1
        # If linewidth=0 is set, twisted part is almost invisible
        kwargs.update(dict(color=color, alpha=alpha))
        if "lw" not in kwargs and "linewidth" not in kwargs:
            kwargs.update(dict(lw=0.1))

        if not allow_twist:
            # Resolve twist
            if (rad_end1 - rad_start1) * (rad_end2 - rad_start2) > 0:
                rad_start2, rad_end2 = rad_end2, rad_start2

        # Create bezier curve path patch
        bezier_curve = BezierCurve(
            rad_start1,
            rad_end1,
            r1,
            rad_start2,
            rad_end2,
            r2,
            height_ratio,
            direction,
            arrow_length_ratio,
            **kwargs,
        )
        self._patches.append(bezier_curve)

    def plotfig(
        self,
        dpi: int = 100,
        *,
        ax: PolarAxes | None = None,
    ) -> Figure:
        """Plot figure

        Parameters
        ----------
        dpi : int, optional
            Figure DPI
        ax : PolarAxes | None
            If None, figure and axes are newly created.

        Returns
        -------
        figure : Figure
            Circos matplotlib figure
        """
        if ax is None:
            # Initialize Figure & PolarAxes
            fig, ax = self._initialize_figure(dpi=dpi)
        else:
            # Check PolarAxes or not
            if not isinstance(ax, PolarAxes):
                ax_class_name = type(ax).__name__
                err_msg = f"Input ax is not PolarAxes (={ax_class_name})."
                raise ValueError(err_msg)
            fig = ax.get_figure()
        self._initialize_polar_axes(ax)

        # Plot all patches
        patches = []
        for patch in self._get_all_patches():
            # Collection cannot handle `zorder`, `hatch`
            # Separate default or user-defined `zorder`, `hatch` property patch
            if patch.get_zorder() == 1 and patch.get_hatch() is None:
                patches.append(patch)
            else:
                ax.add_patch(patch)
        ax.add_collection(PatchCollection(patches, match_original=True))

        # Execute all plot functions
        for plot_func in self._get_all_plot_funcs():
            plot_func(ax)

        return fig  # type: ignore

    def savefig(
        self,
        savefile: str | Path,
        dpi: int = 100,
        pad_inches: float = 0.5,
    ) -> None:
        """Save figure to file

        Parameters
        ----------
        savefile : str | Path
            Save file
        dpi : int, optional
            DPI
        pad_inches : float, optional
            Padding inches
        """
        figure = self.plotfig(dpi=dpi)
        figure.savefig(
            fname=savefile,
            dpi=dpi,
            pad_inches=pad_inches,
            bbox_inches="tight",
        )

    ############################################################
    # Private Method
    ############################################################

    def _check_degree_range(self, start: float, end: float) -> None:
        """Check start-end degree range (`-360 <= start < end <= 360`)

        Parameters
        ----------
        start : float
            Start degree range
        end : float
            End degree range
        """
        min_deg, max_deg = -360, 360
        if not min_deg <= start < end <= max_deg:
            err_msg = "start-end must be "
            err_msg += f"'{min_deg} <= start < end <= {max_deg}' ({start=}, {end=})"
            raise ValueError(err_msg)
        if end - start > max_deg:
            err_msg = f"'end - start' must be less than {max_deg} ({start=}, {end=})"
            raise ValueError(err_msg)

    def _initialize_figure(
        self,
        figsize: tuple[float, float] = (8, 8),
        dpi: int = 100,
    ) -> tuple[Figure, PolarAxes]:
        """Initialize figure

        Parameters
        ----------
        figsize : tuple[float, float], optional
            Figure size
        dpi : int, optional
            Figure DPI

        Returns
        -------
        fig, ax : tuple[Figure, PolarAxes]
            Figure, PolarAxes
        """
        fig = plt.figure(figsize=figsize, dpi=dpi, tight_layout=True)
        ax = fig.add_subplot(projection="polar")
        return fig, ax  # type: ignore

    def _initialize_polar_axes(self, ax: PolarAxes) -> None:
        """Initialize polar axes params

        Parameters
        ----------
        ax : PolarAxes
            PolarAxes
        """
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        # Reason for setting the max radius limit at MAX_R(100) + R_PLOT_MARGIN
        # Because a portion of the patch at the upper boundary of 100 may be missed.
        ax.set_rlim(bottom=config.MIN_R, top=config.MAX_R + config.R_PLOT_MARGIN)

        show_axis = "on" if self._show_axis_for_debug else "off"
        ax.axis(show_axis)

    def _get_all_patches(self) -> list[Patch]:
        """Get all patches from `circos, sector, track`

        Returns
        -------
        all_patches : list[Patch]
            All patches
        """
        circos_patches = self._patches
        sector_patches = list(itertools.chain(*[s.patches for s in self.sectors]))
        track_patches = list(itertools.chain(*[t.patches for t in self.tracks]))
        all_patches = circos_patches + sector_patches + track_patches
        # deepcopy to avoid putting original patch to figure
        return deepcopy(all_patches)

    def _get_all_plot_funcs(self) -> list[Callable[[PolarAxes], None]]:
        """Get all plot functions from `circos, sector, track`

        Returns
        -------
        all_plot_funcs : list[Callable[[PolarAxes], None]]
            All plot functions
        """
        circos_plot_funcs = self._plot_funcs
        sector_plot_funcs = list(itertools.chain(*[s.plot_funcs for s in self.sectors]))
        track_plot_funcs = list(itertools.chain(*[t.plot_funcs for t in self.tracks]))
        all_plot_funcs = circos_plot_funcs + sector_plot_funcs + track_plot_funcs
        return all_plot_funcs
