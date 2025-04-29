from __future__ import annotations

import itertools
import math
import textwrap
import warnings
from collections import defaultdict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np

import plotly.graph_objects as go
from plotly.graph_objs.layout._shape import Shape
from plotly.graph_objs.layout._annotation import Annotation 

from pycirclize_TEST import config, utils
from pycirclize_TEST.parser import Bed, Matrix
from pycirclize_TEST.sector import Sector
from pycirclize_TEST.track import Track
from pycirclize_TEST.patches import PolarSVGPatchBuilder

class Circos:
    """Circos Visualization Class"""

    def __init__(
        self,
        sectors: Mapping[str, int | float | tuple[float, float]],
        start: float = 0,
        end: float = 360,
        *,
        space: float | list[float] = 0,
        endspace: bool = True,
        sector2clockwise: dict[str, bool] | None = None,
        show_axis_for_debug: bool = False,
    ):
        """
        Parameters
        ----------
        sectors : Mapping[str, int | float | tuple[float, float]]
            Sector name & size (or range) dict
        start : float, optional
            Plot start degree (`-360 <= start < end <= 360`)
        end : float, optional
            Plot end degree (`-360 <= start < end <= 360`)
        space : float | list[float], optional
            Space degree(s) between sector
        endspace : bool, optional
            If True, insert space after the end sector
        sector2clockwise : dict[str, bool] | None, optional
            Sector name & clockwise bool dict. By default, `clockwise=True`.
        show_axis_for_debug : bool, optional
            Show axis for position check debugging (Developer option)
        """
        sector2clockwise = {} if sector2clockwise is None else sector2clockwise

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
            space_list = [space] * space_num + [0]
            space_deg_size = space * space_num
        whole_deg_size_without_space = whole_deg_size - space_deg_size
        if whole_deg_size_without_space < 0:
            err_msg = textwrap.dedent(
                f"""
                Too large sector space size is set!!
                Circos Degree Size = {whole_deg_size} ({start} - {end})
                """
                # Total Sector Space Size = {space_deg_size}
                # List of Sector Space Size = {space_list}
                
            )[1:-1]
            raise ValueError(err_msg)

        sector2range = self._to_sector2range(sectors)
        sector_total_size = sum([max(r) - min(r) for r in sector2range.values()])

        rad_pos = math.radians(start)
        
        self._sectors: list[Sector] = []
        for idx, (sector_name, sector_range) in enumerate(sector2range.items()):
            sector_size = max(sector_range) - min(sector_range)
            sector_size_ratio = sector_size / sector_total_size
            deg_size = whole_deg_size_without_space * sector_size_ratio
            rad_size = math.radians(deg_size)
            rad_lim = (rad_pos, rad_pos + rad_size)
            rad_pos += rad_size + math.radians(space_list[idx])
            clockwise = sector2clockwise.get(sector_name, True)
            sector = Sector(sector_name, sector_range, rad_lim, clockwise)
            self._sectors.append(sector)

        self._deg_lim = (start, end)
        self._rad_lim = (math.radians(start), math.radians(end))
        self._show_axis_for_debug = show_axis_for_debug

        # Shapes and annotations for Layout
        self._shapes: list[Shape] = []
        self._annotations: list[Annotation] = []
        # self._ax: PolarAxes | None = None

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
    def initialize_from_bed(
        bed_file: str | Path,
        start: float = 0,
        end: float = 360,
        *,
        space: float | list[float] = 0,
        endspace: bool = True,
        sector2clockwise: dict[str, bool] | None = None,
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
        sector2clockwise : dict[str, bool] | None, optional
            Sector name & clockwise bool dict. By default, `clockwise=True`.

        Returns
        -------
        circos : Circos
            Circos instance initialized from BED file
        """
        records = Bed(bed_file).records
        sectors = {rec.chr: (rec.start, rec.end) for rec in records}
        return Circos(
            sectors,
            start,
            end,
            space=space,
            endspace=endspace,
            sector2clockwise=sector2clockwise,
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


    def axis(self, **kwargs) -> None:
        """Plot axis

        Parameters
        ----------
        **kwargs : dict, optional
            Shape properties (e.g. `fillcolor="red", line=dict(color="darkgreen", width=2, dash="dash", ... ) ...`)
            <https://plotly.com/python/reference/layout/shapes/>
        """
        kwargs = {} if kwargs is None else kwargs

        # Background shape placed behind other shapes (layer="below")
        fc_behind_kwargs = deepcopy(kwargs)
        fc_behind_kwargs.update(config.AXIS_FACE_PARAM)
        self.rect(**fc_behind_kwargs)

        # Edge shape placed in front of other shapes (layer="above")
        ec_front_kwargs = deepcopy(kwargs)
        ec_front_kwargs.update(config.AXIS_EDGE_PARAM)
        self.rect(**ec_front_kwargs)

    def text(
        self,
        text: str,
        *,
        r: float = 0,
        deg: float = 0,
        adjust_rotation: bool = False,
        orientation: str = "horizontal",
        outer: bool = True,
        **kwargs,
    ) -> None:
        """Plot text on the entire circos plot. Uses angular positioning (0-360°).
        Angle is adjusted to Plotly's coordinate system:
            - 0° points upward (Plotly's default)
            - Angles increase clockwise

        Parameters
        ----------
        text : str
            Text content
        r : float, optional
            Radius position (default: 0, centered).
        deg : float, optional
            Degree position (0-360). 0° points upward.
        adjust_rotation : bool, optional
            If True, text rotation is auto set based on `deg` and `orientation`.
        orientation : str, optional
            Text orientation (`horizontal` or `vertical`).
        outer : bool, optional
            If True, text aligns outward from center (for horizontal orientation).
        **kwargs : dict, optional
            Annotation properties (e.g. `font=dict(size=12, color='red')`).
            See: <https://plotly.com/python/reference/layout/annotations/>
        """
        rad = np.radians(deg)
        plotly_rad = -(rad - np.pi/2)  # Convert to Plotly's polar coordinates
        x_pos = r * np.cos(plotly_rad)
        y_pos = r * np.sin(plotly_rad)

        annotation = utils.plot.get_plotly_label_params(rad, adjust_rotation, orientation, outer, **kwargs)
        
        annotation.update({
            "x": x_pos,
            "y": y_pos,
            "text": text,
        })

        self._annotations.append(annotation)

    # def plot_text(ax: PolarAxes) -> None:
    #     ax.text(math.radians(deg), r, text, **kwargs)

    # self._plot_funcs.append(plot_text)

    def line(
        self,
        *,
        r: float | tuple[float, float],
        deg_lim: tuple[float, float] | None = None,
        arc: bool = True,
        **kwargs,
    ) -> None:
        """Plot line

        Parameters
        ----------
        r : float | tuple[float, float]
            Line radius position (0 - 100). If r is float, (r, r) is set.
        deg_lim : tuple[float, float] | None, optional
            Degree limit region (-360 - 360). If None, `circos.deg_lim` is set.
        arc : bool, optional
            If True, plot arc style line for polar projection.
            If False, simply plot linear style line.
        **kwargs : dict, optional
            Patch properties (e.g. `color="red", lw=3, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        deg_lim = self.deg_lim if deg_lim is None else deg_lim
        rad_lim = (math.radians(min(deg_lim)), math.radians(max(deg_lim)))
        r_lim = r if isinstance(r, (tuple, list)) else (r, r)
        LinePatch = ArcLine if arc else Line
        # self._patches.append(LinePatch(rad_lim, r_lim, **kwargs))

    def rect(
        self,
        r_lim: tuple[float, float] = (0, 100),
        deg_lim: tuple[float, float] | None = None,
        **kwargs,
    ) -> None:
        """Plot a rectangle spanning angular and radial ranges.
        Angle is adjusted to Plotly's coordinate system:
            - 0° points upward (Plotly's default)
            - Angles increase clockwise

        Parameters
        ----------
        r_lim : tuple[float, float]
            Radial limits (min, max) between 0 (center) and 100 (outer edge).
        deg_lim : tuple[float, float] | None, optional
            Angular limits in degrees (-360 to 360). If None, uses `circos.deg_lim`.
        **kwargs : dict, optional
            Shape properties (e.g. `fillcolor="red", line=dict(color="blue", width=2)`)
            See: <https://plotly.com/python/reference/layout/shapes/>
        """
        deg_lim = self.deg_lim if deg_lim is None else deg_lim
        rad_lim = (np.radians(min(deg_lim)), np.radians(max(deg_lim)))
        
        # Convert to Plotly's coordinate system
        plotly_rad_lim = [-(rad - np.pi/2) for rad in rad_lim]
        min_rad, max_rad = min(plotly_rad_lim), max(plotly_rad_lim)
        
        # Build rectangle path
        radr = (min_rad, min(r_lim))
        width = max_rad - min_rad
        height = max(r_lim) - min(r_lim)
        
        path = PolarSVGPatchBuilder.arc_rectangle(radr, width, height)
        shape = utils.plot.build_plotly_shape(path, **kwargs)
        self._shapes.append(shape)

    # def link(
    #     self,
    #     sector_region1: tuple[str, float, float],
    #     sector_region2: tuple[str, float, float],
    #     r1: float | None = None,
    #     r2: float | None = None,
    #     *,
    #     color: str = "grey",
    #     alpha: float = 0.5,
    #     height_ratio: float = 0.5,
    #     direction: int = 0,
    #     arrow_length_ratio: float = 0.05,
    #     allow_twist: bool = True,
    #     **kwargs,
    # ) -> None:
    #     """Plot link to specified region within or between sectors

    #     Parameters
    #     ----------
    #     sector_region1 : tuple[str, float, float]
    #         Link sector region1 (name, start, end)
    #     sector_region2 : tuple[str, float, float]
    #         Link sector region2 (name, start, end)
    #     r1 : float | None, optional
    #         Link radius end position for sector_region1.
    #         If None, lowest radius position of track in target sector is set.
    #     r2 : float | None, optional
    #         Link radius end position for sector_region2.
    #         If None, lowest radius position of track in target sector is set.
    #     color : str, optional
    #         Link color
    #     alpha : float, optional
    #         Link color alpha (transparency) value
    #     height_ratio : float, optional
    #         Bezier curve height ratio
    #     direction : int, optional
    #         `0`: No direction edge shape (Default)
    #         `1`: Forward direction arrow edge shape (region1 -> region2)
    #         `-1`: Reverse direction arrow edge shape (region1 <- region2)
    #         `2`: Bidirectional arrow edge shape (region1 <-> region2)
    #     arrow_length_ratio : float, optional
    #         Direction arrow length ratio
    #     allow_twist : bool, optional
    #         If False, twisted link is automatically resolved.
    #         <http://circos.ca/documentation/tutorials/links/twists/images>
    #     **kwargs : dict, optional
    #         Patch properties (e.g. `ec="red", lw=1.0, hatch="//", ...`)
    #         <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
    #     """
    #     # Set data for plot link
    #     name1, start1, end1 = sector_region1
    #     name2, start2, end2 = sector_region2
    #     sector1, sector2 = self.get_sector(name1), self.get_sector(name2)
    #     r1 = sector1.get_lowest_r() if r1 is None else r1
    #     r2 = sector2.get_lowest_r() if r2 is None else r2
    #     rad_start1, rad_end1 = sector1.x_to_rad(start1), sector1.x_to_rad(end1)
    #     rad_start2, rad_end2 = sector2.x_to_rad(start2), sector2.x_to_rad(end2)

    #     # Set patch kwargs & default linewidth as 0.1
    #     # If linewidth=0 is set, twisted part is almost invisible
    #     kwargs.update(dict(color=color, alpha=alpha))
    #     if "lw" not in kwargs and "linewidth" not in kwargs:
    #         kwargs.update(dict(lw=0.1))

    #     if not allow_twist:
    #         # Resolve twist
    #         if (rad_end1 - rad_start1) * (rad_end2 - rad_start2) > 0:
    #             rad_start2, rad_end2 = rad_end2, rad_start2

    #     # Create bezier curve path patch
    #     bezier_curve_link = BezierCurveLink(
    #         rad_start1,
    #         rad_end1,
    #         r1,
    #         rad_start2,
    #         rad_end2,
    #         r2,
    #         height_ratio,
    #         direction,
    #         arrow_length_ratio,
    #         **kwargs,
    #     )
    #     self._patches.append(bezier_curve_link)

    # def link_line(
    #     self,
    #     sector_pos1: tuple[str, float],
    #     sector_pos2: tuple[str, float],
    #     r1: float | None = None,
    #     r2: float | None = None,
    #     *,
    #     color: str = "black",
    #     height_ratio: float = 0.5,
    #     direction: int = 0,
    #     arrow_height: float = 3.0,
    #     arrow_width: float = 2.0,
    #     **kwargs,
    # ) -> None:
    #     """Plot link line to specified position within or between sectors

    #     Parameters
    #     ----------
    #     sector_pos1 : tuple[str, float]
    #         Link line sector position1 (name, position)
    #     sector_pos2 : tuple[str, float]
    #         Link line sector position2 (name, position)
    #     r1 : float | None, optional
    #         Link line radius end position for sector_pos1.
    #         If None, lowest radius position of track in target sector is set.
    #     r2 : float | None, optional
    #         Link line radius end position for sector_pos2.
    #         If None, lowest radius position of track in target sector is set.
    #     color : str, optional
    #         Link line color
    #     height_ratio : float, optional
    #         Bezier curve height ratio
    #     direction : int, optional
    #         `0`: No direction edge shape (Default)
    #         `1`: Forward direction arrow edge shape (pos1 -> pos2)
    #         `-1`: Reverse direction arrow edge shape (pos1 <- pos2)
    #         `2`: Bidirectional arrow edge shape (pos1 <-> pos2)
    #     arrow_height : float, optional
    #         Arrow height size (Radius unit)
    #     arrow_width : float, optional
    #         Arrow width size (Degree unit)
    #     **kwargs : dict, optional
    #         Patch properties (e.g. `lw=1.0, ls="dashed", ...`)
    #         <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
    #     """
    #     # Set data for plot link
    #     name1, pos1 = sector_pos1
    #     name2, pos2 = sector_pos2
    #     sector1, sector2 = self.get_sector(name1), self.get_sector(name2)
    #     r1 = sector1.get_lowest_r() if r1 is None else r1
    #     r2 = sector2.get_lowest_r() if r2 is None else r2
    #     rad_pos1, rad_pos2 = sector1.x_to_rad(pos1), sector2.x_to_rad(pos2)

    #     kwargs.update(color=color)

    #     bezier_curve_line = BezierCurveLine(
    #         rad_pos1,
    #         r1,
    #         rad_pos2,
    #         r2,
    #         height_ratio,
    #         direction,
    #         arrow_height,
    #         arrow_width,
    #         **kwargs,
    #     )
    #     self._patches.append(bezier_curve_line)

    def colorbar(
        self,
        bounds: tuple[float, float, float, float] = (1.02, 0.3, 0.02, 0.4),
        *,
        vmin: float = 0,
        vmax: float = 1,
        cmap: str | Colormap = "bwr",
        orientation: str = "vertical",
        label: str | None = None,
        colorbar_kws: dict[str, Any] | None = None,
        label_kws: dict[str, Any] | None = None,
        tick_kws: dict[str, Any] | None = None,
    ) -> None:
        """Plot colorbar

        Parameters
        ----------
        bounds : tuple[float, float, float, float], optional
            Colorbar bounds tuple (`x`, `y`, `width`, `height`)
        vmin : float, optional
            Colorbar min value
        vmax : float, optional
            Colorbar max value
        cmap : str | Colormap, optional
            Colormap (e.g. `viridis`, `Spectral`, `Reds`, `Greys`)
            <https://matplotlib.org/stable/tutorials/colors/colormaps.html>
        orientation : str, optional
            Colorbar orientation (`vertical`|`horizontal`)
        label : str | None, optional
            Colorbar label. If None, no label shown.
        colorbar_kws : dict[str, Any] | None, optional
            Colorbar properties (e.g. `dict(format="%.1f", ...)`)
            <https://matplotlib.org/stable/api/colorbar_api.html>
        label_kws : dict[str, Any] | None, optional
            Text properties (e.g. `dict(size=15, color="red", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        tick_kws : dict[str, Any] | None, optional
            Axes.tick_params properties (e.g. `dict(labelsize=12, colors="red", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.tick_params.html>
        """
        colorbar_kws = {} if colorbar_kws is None else deepcopy(colorbar_kws)
        label_kws = {} if label_kws is None else deepcopy(label_kws)
        tick_kws = {} if tick_kws is None else deepcopy(tick_kws)

        def plot_colorbar(ax: PolarAxes) -> None:
            axin: Axes = ax.inset_axes(bounds)
            norm = Normalize(vmin=vmin, vmax=vmax)
            cb = Colorbar(
                axin,
                cmap=cmap,  # type: ignore
                norm=norm,
                orientation=orientation,  # type: ignore
                **colorbar_kws,
            )
            axin.tick_params(**tick_kws)
            if label:
                cb.set_label(label, **label_kws)

        # self._plot_funcs.append(plot_colorbar)

    def plotfig(
        self,
        dpi: int = 100,
        figsize: tuple[float, float] = (8, 8),
        **kwargs,
    ) -> go.Figure:
        """Create the Plotly Circos-style figure.

        Parameters
        ----------
        dpi : int, optional
            Dots per inch (used to scale figsize)
        figsize : tuple[float, float], optional
            Size of figure in inches (width, height)
        **kwargs : dict
            Additional layout settings to override defaults

        Returns
        -------
        fig : go.Figure
            Plotly figure object
        """
        layout_dict = self._initialize_plotly_layout(figsize=figsize, dpi=dpi)
        layout_dict.update(kwargs)

        layout_dict['shapes'] = self._get_all_shapes()
        layout_dict['annotations'] = self._get_all_annotations()

        return go.Figure(layout=go.Layout(layout_dict))

    def savefig(
        self,
        savefile: str | Path,
        *,
        dpi: int = 100,
        figsize: tuple[float, float] = (8, 8),
        pad_inches: float = 0.5,
    ) -> None:
        """Save figure to file

        Parameters
        ----------
        savefile : str | Path
            Save file (`*.png`|`*.jpg`|`*.svg`|`*.pdf`)
        dpi : int, optional
            DPI
        figsize : tuple[float, float], optional
            Figure size
        pad_inches : float, optional
            Padding inches

        Warnings
        --------
        To plot a figure that settings a user-defined legend, subtracks, or annotations,
        call `fig.savefig()` instead of `gv.savefig()`.
        """
        # fig = self.plotfig(dpi=dpi, figsize=figsize)
        # fig.savefig(
        #     fname=savefile,  # type: ignore
        #     dpi=dpi,
        #     pad_inches=pad_inches,
        #     bbox_inches="tight",
        # )
        # # Clear & close figure to suppress memory leak
        # if config.clear_savefig:
        #     fig.clear()
        #     plt.close(fig)

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

    def _to_sector2range(
        self,
        sectors: Mapping[str, int | float | tuple[float, float]],
    ) -> dict[str, tuple[float, float]]:
        """Convert sectors to sector2range"""
        sector2range: dict[str, tuple[float, float]] = {}
        for name, value in sectors.items():
            if isinstance(value, (tuple, list)):
                sector_start, sector_end = value
                if not sector_start < sector_end:
                    err_msg = f"{sector_end=} must be larger than {sector_start=}."
                    raise ValueError(err_msg)
                sector2range[name] = (sector_start, sector_end)
            else:
                sector2range[name] = (0, value)
        return sector2range

    @staticmethod
    def _initialize_plotly_layout(
        figsize: tuple[float, float] = (8, 8),
        dpi: int = 100,
    ) -> dict:
        """Initialize default Plotly layout based on config and figure size."""
        width = int(figsize[0] * dpi)
        height = int(figsize[1] * dpi)

        layout = deepcopy(config.plotly_layout_defaults)

        layout['width'] = width
        layout['height'] = height
        layout['xaxis']['range'] = config.AXIS_RANGE
        layout['yaxis']['range'] = config.AXIS_RANGE

        return layout

    def _get_all_shapes(self) -> list[dict]:
        """Gather all shape dictionaries from self, sectors, and tracks."""
        circos_shapes = self._shapes
        sector_shapes = list(itertools.chain(*[s._shapes for s in self.sectors]))
        track_shapes = list(itertools.chain(*[t._shapes for t in self.tracks]))
        return circos_shapes + sector_shapes + track_shapes

    def _get_all_annotations(self) -> list[dict]:
        """Gather all annotation dictionaries from self, sectors, and tracks."""
        circos_ann = self._annotations
        sector_ann = list(itertools.chain(*[s._annotations for s in self.sectors]))
        track_ann = list(itertools.chain(*[t._annotations for t in self.tracks]))
        return circos_ann + sector_ann + track_ann


    def _adjust_annotation(self) -> None:
        """Adjust annotation text position"""
        # Get sorted annotation list for position adjustment
        ann_list = self._get_sorted_ann_list()
        if len(ann_list) == 0 or config.ann_adjust.max_iter <= 0:
            return
        if len(ann_list) > config.ann_adjust.limit:
            warn_msg = f"Too many annotations(={len(ann_list)}). Annotation position adjustment is not done."  # noqa: E501
            warnings.warn(warn_msg)
            return

        def get_ann_window_extent(ann: Annotation) -> Bbox:
            return Text.get_window_extent(ann).expanded(*config.ann_adjust.expand)

        # Iterate annotation position adjustment
        self.ax.figure.draw_without_rendering()  # type: ignore
        ann2rad_shift_candidates = self._get_ann2rad_shift_candidates(ann_list)
        for idx, ann in enumerate(ann_list[1:], 1):
            orig_rad, orig_r = ann.xyann
            ann_bbox = get_ann_window_extent(ann)
            adj_ann_list = ann_list[:idx]
            adj_ann_bboxes = [get_ann_window_extent(ann) for ann in adj_ann_list]

            # Adjust radian position
            iter, max_iter = 0, config.ann_adjust.max_iter
            if utils.plot.is_ann_rad_shift_target_loc(orig_rad):
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
        #     utils.plot.plot_bbox(get_ann_window_extent(ann), self.ax)

    def _get_sorted_ann_list(self) -> list[Annotation]:
        """Sorted annotation list

        Sorting per 4 sections for adjusting annotation text position
        """
        ann_list = [t for t in self.ax.texts if isinstance(t, Annotation)]
        loc2ann_list: dict[str, list[Annotation]] = defaultdict(list)
        for ann in ann_list:
            loc = utils.plot.get_loc(ann.xyann[0])
            loc2ann_list[loc].append(ann)

        def sort_by_ann_rad(ann: Annotation):
            return utils.plot.degrees(ann.xyann[0])

        return (
            sorted(loc2ann_list["upper-right"], key=sort_by_ann_rad, reverse=True)
            + sorted(loc2ann_list["lower-right"], key=sort_by_ann_rad, reverse=False)
            + sorted(loc2ann_list["lower-left"], key=sort_by_ann_rad, reverse=True)
            + sorted(loc2ann_list["upper-left"], key=sort_by_ann_rad, reverse=False)
        )

    def _get_ann2rad_shift_candidates(
        self, ann_list: list[Annotation]
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
