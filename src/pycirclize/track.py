from __future__ import annotations

import math
from copy import deepcopy
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

import matplotlib as mpl
import numpy as np
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.SeqFeature import SeqFeature
from matplotlib.colors import Colormap, Normalize
from matplotlib.patches import Patch
from matplotlib.projections.polar import PolarAxes

from pycirclize import config, utils
from pycirclize.patches import ArcArrow, ArcLine, ArcRectangle

if TYPE_CHECKING:
    # Avoid Sector <-> Track circular import error at runtime
    from pycirclize.sector import Sector


class Track:
    """Circos Track Class"""

    def __init__(
        self,
        name: str,
        r_lim: tuple[float, float],
        r_pad_ratio: float,
        parent_sector: Sector,
    ):
        """
        Parameters
        ----------
        name : str
            Track name
        r_lim : tuple[float, float]
            Track radius limit region
        r_pad_ratio : float
            Track padding ratio for plot data
        parent_sector : Sector
            Parent sector of track
        """
        # Track params
        self._name = name
        self._r_lim = r_lim
        self._r_pad_ratio = r_pad_ratio
        # Inherited from parent sector
        self._parent_sector = parent_sector
        self._rad_lim = parent_sector.rad_lim
        self._start = parent_sector.start
        self._end = parent_sector.end

        # Plot data and functions
        self._patches: list[Patch] = []
        self._plot_funcs: list[Callable[[PolarAxes], None]] = []

    ############################################################
    # Property
    ############################################################

    @property
    def name(self) -> str:
        """Track name"""
        return self._name

    @property
    def size(self) -> float:
        """Track size (x coordinate)"""
        return self.end - self.start

    @property
    def start(self) -> float:
        """Track start position (x coordinate)"""
        return self._start

    @property
    def end(self) -> float:
        """Track end position (x coordinate)"""
        return self._end

    @property
    def r_size(self) -> float:
        """Track radius size"""
        return max(self.r_lim) - min(self.r_lim)

    @property
    def r_lim(self) -> tuple[float, float]:
        """Track radius limit"""
        return self._r_lim

    @property
    def r_plot_size(self) -> float:
        """Track radius size for plot data (`r_size` with padding)"""
        return max(self.r_plot_lim) - min(self.r_plot_lim)

    @property
    def r_plot_lim(self) -> tuple[float, float]:
        """Track radius limit for plot data (`r_lim` with padding)"""
        edge_pad_size = (self.r_size * self._r_pad_ratio) / 2
        min_plot_r = min(self.r_lim) + edge_pad_size
        max_plot_r = max(self.r_lim) - edge_pad_size
        return (min_plot_r, max_plot_r)

    @property
    def rad_size(self) -> float:
        """Track radian size"""
        return max(self.rad_lim) - min(self.rad_lim)

    @property
    def rad_lim(self) -> tuple[float, float]:
        """Track radian limit"""
        return self._rad_lim

    @property
    def deg_size(self) -> float:
        """Track degree size"""
        return max(self.deg_lim) - min(self.deg_lim)

    @property
    def deg_lim(self) -> tuple[float, float]:
        """Track degree limit"""
        return (math.degrees(min(self.rad_lim)), math.degrees(max(self.rad_lim)))

    @property
    def parent_sector(self) -> Sector:
        """Parent sector"""
        return self._parent_sector

    @property
    def clockwise(self) -> bool:
        """Track coordinate direction"""
        return self.parent_sector.clockwise

    @property
    def patches(self) -> list[Patch]:
        """Plot patches"""
        return self._patches

    @property
    def plot_funcs(self) -> list[Callable[[PolarAxes], None]]:
        """Plot functions"""
        return self._plot_funcs

    ############################################################
    # Public Method
    ############################################################

    def x_to_rad(self, x: float, ignore_range_error: bool = False) -> float:
        """Convert x coordinate to radian in track start-end range

        Parameters
        ----------
        x : float
            X coordinate
        ignore_range_error : bool
            Ignore x coordinate range error

        Returns
        -------
        rad : float
            Radian coordinate
        """
        return self.parent_sector.x_to_rad(x, ignore_range_error)

    def axis(self, **kwargs) -> None:
        """Plot axis

        By default, simple black axis params(`fc="none", ec="black", lw=0.5`) are set.

        Parameters
        ----------
        **kwargs : dict, optional
            Patch properties (e.g. `fc="tomato", ec="blue", hatch="//"`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Set default params
        kwargs = utils.plot.set_axis_default_kwargs(**kwargs)

        # Axis facecolor placed behind other patches (zorder=0.99)
        fc_behind_kwargs = {**kwargs, **config.AXIS_FACE_PARAM}
        self.rect(self.start, self.end, ignore_pad=True, **fc_behind_kwargs)

        # Axis edgecolor placed in front of other patches (zorder=1.01)
        ec_front_kwargs = {**kwargs, **config.AXIS_EDGE_PARAM}
        self.rect(self.start, self.end, ignore_pad=True, **ec_front_kwargs)

    def text(
        self,
        text: str,
        x: float | None = None,
        r: float | None = None,
        *,
        adjust_rotation: bool = True,
        orientation: str = "horizontal",
        ignore_range_error: bool = False,
        **kwargs,
    ) -> None:
        """Plot text

        Parameters
        ----------
        text : str
            Text content
        x : float | None
            X position. If None, track center x position is set.
        r : float | None
            Radius position. If None, track center radius position is set.
        adjust_rotation : bool, optional
            If True, text rotation is auto set based on `x` and `orientation` params.
        orientation : str, optional
            Text orientation (`horizontal` or `vertical`)
            If adjust_rotation=True, orientation is used for rotation calculation.
        ignore_range_error : bool, optional
            If True, ignore x position range error
            (ErrorCase: `not track.start <= x <= track.end`)
        **kwargs : dict, optional
            Text properties (e.g. `size=12, color="red", va="center", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        # If value is None, center position is set.
        x = (self.start + self.end) / 2 if x is None else x
        r = sum(self.r_lim) / 2 if r is None else r

        rad = self.x_to_rad(x, ignore_range_error)
        if adjust_rotation:
            params = utils.plot.get_label_params_by_rad(
                rad, orientation, only_rotation=True
            )
            kwargs.update(params)

        if "ha" not in kwargs and "horizontalalignment" not in kwargs:
            kwargs.update(dict(ha="center"))
        if "va" not in kwargs and "verticalalignment" not in kwargs:
            kwargs.update(dict(va="center"))

        def plot_text(ax: PolarAxes) -> None:
            ax.text(rad, r, text, **kwargs)

        self._plot_funcs.append(plot_text)

    def rect(
        self,
        start: float,
        end: float,
        *,
        r_lim: tuple[float, float] | None = None,
        ignore_pad: bool = False,
        **kwargs,
    ) -> None:
        """Plot rectangle

        Parameters
        ----------
        start : float
            Start position (x coordinate)
        end : float
            End position (x coordinate)
        r_lim : tuple[float, float] | None, optional
            Radius limit range.
            If None, `track.r_lim` (ignore_pad=False) or
            `track.r_plot_lim` (ignore_pad=True) is set.
        ignore_pad : bool, optional
            If True, ignore track padding setting.
            If `r_lim` param is set by user, this option not works.
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="blue", lw=1.0, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        rad_rect_start = self.x_to_rad(start)
        rad_rect_end = self.x_to_rad(end)
        rad = min(rad_rect_start, rad_rect_end)
        width = abs(rad_rect_end - rad_rect_start)
        if r_lim is not None:
            if not min(self.r_lim) <= min(r_lim) < max(r_lim) <= max(self.r_lim):
                raise ValueError(f"{r_lim=} is invalid track range.\n{self}")
            radr, height = (rad, min(r_lim)), max(r_lim) - min(r_lim)
        elif ignore_pad:
            radr, height = (rad, min(self.r_lim)), self.r_size
        else:
            radr, height = (rad, min(self.r_plot_lim)), self.r_plot_size
        arc_rect = ArcRectangle(radr, width, height, **kwargs)
        self._patches.append(arc_rect)

    def arrow(
        self,
        start: float,
        end: float,
        *,
        r_lim: tuple[float, float] | None = None,
        head_length: float = 2,
        shaft_ratio: float = 0.5,
        **kwargs,
    ) -> None:
        """Plot arrow

        Parameters
        ----------
        start : float
            Start position (x coordinate)
        end : float
            End position (x coordinate)
        r_lim : tuple[float, float] | None, optional
            Radius limit range. If None, `track.r_lim` is set.
        head_length : float, optional
            Arrow head length (Degree unit)
        shaft_ratio : float, optional
            Arrow shaft ratio (0 - 1.0)
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="blue", lw=1.0, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        rad_arrow_start = self.x_to_rad(start)
        rad_arrow_end = self.x_to_rad(end)
        if r_lim is None:
            r, dr = min(self.r_plot_lim), self.r_plot_size
        else:
            if not min(self.r_lim) <= min(r_lim) < max(r_lim) <= max(self.r_lim):
                raise ValueError(f"{r_lim=} is invalid track range.\n{self}")
            r, dr = min(r_lim), max(r_lim) - min(r_lim)
        arc_arrow = ArcArrow(
            rad=rad_arrow_start,
            r=r,
            drad=rad_arrow_end - rad_arrow_start,
            dr=dr,
            head_length=math.radians(head_length),
            shaft_ratio=shaft_ratio,
            **kwargs,
        )
        self._patches.append(arc_arrow)

    def xticks(
        self,
        x: list[int] | list[float] | np.ndarray,
        labels: list[str] | None = None,
        *,
        tick_length: float = 2,
        outer: bool = True,
        show_bottom_line: bool = False,
        label_size: float = 8,
        label_margin: float = 0.5,
        label_orientation: str = "horizontal",
        line_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Plot xticks & labels on user-specified position

        If you want to plot xticks and their position labels at regular intervals,
        it is recommended to use `track.xticks_by_interval()` instead.

        Parameters
        ----------
        x : list[int] | list[float] | np.ndarray
            X coordinates
        labels : list[str] | None, optional
            Labels on xticks. If None, only plot ticks line.
        tick_length : float, optional
            Tick length (Radius unit)
        outer : bool, optional
            If True, show ticks on outer. If False, show ticks on inner.
        show_bottom_line : bool, optional
            If True, show bottom line.
        label_size : float, optional
            Label size
        label_margin : float, optional
            Label margin size
        label_orientation : str, optional
            Label orientation (`horizontal` or `vertical`)
        line_kws : dict[str, Any] | None, optional
            Patch properties (e.g. `dict(ec="red", lw=1, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        text_kws : dict[str, Any] | None, optional
            Text properties (e.g. `dict(color="red", alpha=0.5, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        line_kws = {} if line_kws is None else deepcopy(line_kws)
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        # Check list length of x & labels
        labels = [""] * len(x) if labels is None else labels
        if len(x) != len(labels):
            err_msg = f"List length is not match ({len(x)=}, {len(labels)=})"
            raise ValueError(err_msg)

        # Plot xticks & labels
        r = max(self.r_lim) if outer else min(self.r_lim)
        tick_r_lim = (r, r + tick_length) if outer else (r - tick_length, r)
        for x_pos, label in zip(x, labels):
            # Plot xticks
            if tick_length > 0:
                self._simpleline((x_pos, x_pos), tick_r_lim, **line_kws)
            # Plot labels
            if label != "":
                rad = self.x_to_rad(x_pos)
                if outer:
                    adj_r = max(tick_r_lim) + label_margin
                else:
                    adj_r = min(tick_r_lim) - label_margin
                params = utils.plot.get_label_params_by_rad(
                    rad, label_orientation, outer
                )
                text_kws.update({**params, **dict(size=label_size)})
                self.text(label, x_pos, adj_r, adjust_rotation=False, **text_kws)

        # Plot bottom line
        if show_bottom_line:
            self._simpleline((self.start, self.end), (r, r), **line_kws)

    def xticks_by_interval(
        self,
        interval: int | float,
        *,
        tick_length: float = 2,
        outer: bool = True,
        show_bottom_line: bool = False,
        show_label: bool = True,
        label_size: float = 8,
        label_margin: float = 0.5,
        label_orientation: str = "horizontal",
        label_formatter: Callable[[float], str] | None = None,
        line_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Plot xticks & position labels by user-specified interval

        `track.xticks_by_interval()` is high-level API function of `track.xticks()`.
        If you want to plot xticks and their labels in any position you like,
        use `track.xticks()` instead.

        Parameters
        ----------
        interval : int | float
            Xticks interval
        tick_length : float, optional
            Tick length (Radius unit)
        outer : bool, optional
            If True, show ticks on outer. If False, show ticks on inner.
        show_bottom_line : bool, optional
            If True, show bottom line.
        show_label : bool, optional
            If True, show label of xticks interval position.
        label_size : float, optional
            Label size
        label_margin : float, optional
            Label margin size
        label_orientation : str, optional
            Label orientation (`horizontal` or `vertical`)
        label_formatter : Callable[[float], str] | None, optional
            User-defined function for label format. (e.g. `1000 -> '1.0 Kb'`)
        line_kws : dict[str, Any] | None, optional
            Patch properties (e.g. `dict(ec="red", lw=1, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        text_kws : dict[str, Any] | None, optional
            Text properties (e.g. `dict(color="red", alpha=0.5, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        line_kws = {} if line_kws is None else deepcopy(line_kws)
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        # Setup xtick positions
        x_list = []
        start_pos, end_pos = self.start - (self.start % interval), self.end + interval
        for x in np.arange(start_pos, end_pos, interval):
            if self.start <= x <= self.end:
                x = int(x) if isinstance(interval, int) else float(x)
                x_list.append(x)

        # Setup xticks labels
        labels = None
        if show_label:
            map_func = str if label_formatter is None else label_formatter
            labels = list(map(map_func, x_list))

        # Plot xticks by user-specified interval
        self.xticks(
            x=x_list,
            labels=labels,
            tick_length=tick_length,
            outer=outer,
            show_bottom_line=show_bottom_line,
            label_size=label_size,
            label_margin=label_margin,
            label_orientation=label_orientation,
            line_kws=line_kws,
            text_kws=text_kws,
        )

    def yticks(
        self,
        y: list[int] | list[float] | np.ndarray,
        labels: list[str] | None = None,
        *,
        vmin: float = 0,
        vmax: float | None = None,
        side: str = "right",
        tick_length: float = 1,
        label_size: float = 8,
        label_margin: float = 0.5,
        line_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Plot yticks & labels on user-specified position

        Parameters
        ----------
        y : list[int] | list[float] | np.ndarray
            Y coordinates
        labels : list[str] | None, optional
            Labels on yticks. If None, only plot ticks line.
        vmin : float, optional
            Y min value
        vmax : float | None, optional
            Y max value. If None, `max(y)` is set.
        side : str, optional
            Ticks side position (`right` or `left`)
        tick_length : float, optional
            Tick length (Degree unit)
        label_size : float, optional
            Label size
        label_margin : float, optional
            Label margin size
        line_kws : dict[str, Any] | None, optional
            Patch properties (e.g. `dict(ec="red", lw=1, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        text_kws : dict[str, Any] | None, optional
            Text properties (e.g. `dict(color="red", alpha=0.5, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        line_kws = {} if line_kws is None else deepcopy(line_kws)
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        # Check y, labels list length
        labels = [""] * len(y) if labels is None else labels
        if len(y) != len(labels):
            err_msg = f"List length is not match ({len(y)=}, {len(labels)=})"
            raise ValueError(err_msg)
        # Check side value
        if side not in ("right", "left"):
            raise ValueError(f"{side=} is invalid ('right' or 'left').")
        # Set vmax & check if y is in min-max range
        vmax = max(y) if vmax is None else vmax
        self._check_value_min_max(y, vmin, vmax)
        # Temporarily set clockwise=True in this method
        original_clockwise = self.clockwise
        self.parent_sector._clockwise = True

        # Plot yticks & labels
        r = [self._y_to_r(v, vmin, vmax) for v in y]
        for r_pos, label in zip(r, labels):
            # Set plot properties
            x_tick_length = (self.size / self.deg_size) * tick_length
            x_label_margin = (self.size / self.deg_size) * label_margin
            if side == "right":
                x_lim = (self.end, self.end + x_tick_length)
                x_text = self.end + (x_tick_length + x_label_margin)
                deg_text = math.degrees(self.x_to_rad(x_text, True))
                is_lower_loc = -270 <= deg_text < -90 or 90 <= deg_text < 270
                ha = "right" if is_lower_loc else "left"
            elif side == "left":
                x_lim = (self.start, self.start - x_tick_length)
                x_text = self.start - (x_tick_length + x_label_margin)
                deg_text = math.degrees(self.x_to_rad(x_text, True))
                is_lower_loc = -270 <= deg_text < -90 or 90 <= deg_text < 270
                ha = "left" if is_lower_loc else "right"
            # Plot yticks
            if tick_length > 0:
                self._simpleline(x_lim, (r_pos, r_pos), **line_kws)
            # Plot ylabels
            if label != "":
                va = "center_baseline"
                _text_kws = deepcopy(text_kws)
                _text_kws.update(
                    dict(ha=ha, va=va, rotation_mode="anchor", size=label_size)
                )
                self.text(label, x_text, r_pos, ignore_range_error=True, **_text_kws)

        # Restore clockwise to original value
        self.parent_sector._clockwise = original_clockwise

    def grid(
        self,
        y_grid_num: int | None = 6,
        x_grid_interval: float | None = None,
        **kwargs,
    ) -> None:
        """Plot grid

        By default, `color="grey", alpha=0.5, zorder=0` line params are set.

        Parameters
        ----------
        y_grid_num : int | None, optional
            Y-axis grid line number. If None, y-axis grid line is not shown.
        x_grid_interval : float | None, optional
            X-axis grid line interval. If None, x-axis grid line is not shown.
        **kwargs : dict, optional
            Axes.plot properties (e.g. `color="red", lw=0.5, ls="--", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html>
        """
        # Check argument values
        if y_grid_num is not None and not y_grid_num >= 2:
            raise ValueError(f"{y_grid_num=} is invalid (y_grid_num >= 2).")
        if x_grid_interval is not None and not x_grid_interval > 0:
            raise ValueError(f"{x_grid_interval=} is invalid (x_grid_interval > 0).")

        # Set default grid line properties
        default_props = dict(color="grey", alpha=0.5, zorder=0)
        for name, value in default_props.items():
            if name not in kwargs:
                kwargs.update({name: value})

        # Plot y-axis grid line
        if y_grid_num is not None:
            vmin, vmax = 0, y_grid_num - 1
            for y_grid_idx in range(y_grid_num):
                x = [self.start, self.end]
                y: list[float] = [y_grid_idx, y_grid_idx]
                self.line(x, y, vmin=vmin, vmax=vmax, **kwargs)

        # Plot x-axis grid line
        if x_grid_interval is not None:
            vmin, vmax = 0, 1.0
            x_grid_idx = 0
            while True:
                x_pos = self.start + (x_grid_interval * x_grid_idx)
                if x_pos > self.end:
                    break
                x, y = [x_pos, x_pos], [vmin, vmax]
                self.line(x, y, vmin=vmin, vmax=vmax, **kwargs)
                x_grid_idx += 1

    def line(
        self,
        x: list[float] | np.ndarray,
        y: list[float] | np.ndarray,
        *,
        vmin: float = 0,
        vmax: float | None = None,
        arc: bool = True,
        **kwargs,
    ) -> None:
        """Plot line

        Parameters
        ----------
        x : list[float] | np.ndarray
            X coordinates
        y : list[float] | np.ndarray
            Y coordinates
        vmin : float, optional
            Y min value
        vmax : float | None, optional
            Y max value. If None, `max(y)` is set.
        arc : bool, optional
            If True, plot arc style line for polar projection.
            If False, simply plot linear style line.
        **kwargs : dict, optional
            Axes.plot properties (e.g. `color="red", lw=0.5, ls="--", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html>
        """
        # Check x, y list length
        if len(x) != len(y):
            err_msg = f"List length is not match ({len(x)=}, {len(y)=})"
            raise ValueError(err_msg)

        # Convert (x, y) to (rad, r)
        rad = list(map(self.x_to_rad, x))
        vmax = max(y) if vmax is None else vmax
        self._check_value_min_max(y, vmin, vmax)
        r = [self._y_to_r(v, vmin, vmax) for v in y]

        if arc:
            # Convert linear line to arc line (rad, r) points
            plot_rad, plot_r = self._to_arc_radr(rad, r)
        else:
            plot_rad, plot_r = rad, r

        def plot_line(ax: PolarAxes) -> None:
            ax.plot(plot_rad, plot_r, **kwargs)

        self._plot_funcs.append(plot_line)

    def scatter(
        self,
        x: list[float] | np.ndarray,
        y: list[float] | np.ndarray,
        *,
        vmin: float = 0,
        vmax: float | None = None,
        **kwargs,
    ) -> None:
        """Plot scatter

        Parameters
        ----------
        x : list[float] | np.ndarray
            X position list
        y : list[float] | np.ndarray
            Y position list
        vmin : float, optional
            Y min value
        vmax : float | None, optional
            Y max value. If None, `max(y)` is set.
        **kwargs : dict, optional
            Axes.scatter properties (e.g. `ec="black", lw=1.0, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.scatter.html>
        """
        # Check x, y list length
        if len(x) != len(y):
            err_msg = f"List length is not match ({len(x)=}, {len(y)=})"
            raise ValueError(err_msg)

        # Convert (x, y) to (rad, r)
        rad = list(map(self.x_to_rad, x))
        vmax = max(y) if vmax is None else vmax
        self._check_value_min_max(y, vmin, vmax)
        r = [self._y_to_r(v, vmin, vmax) for v in y]

        def plot_scatter(ax: PolarAxes) -> None:
            ax.scatter(rad, r, **kwargs)

        self._plot_funcs.append(plot_scatter)

    def bar(
        self,
        x: list[float] | np.ndarray,
        height: list[float] | np.ndarray,
        width: float = 0.8,
        bottom: float | list[float] | np.ndarray = 0,
        align: str = "center",
        *,
        vmin: float = 0,
        vmax: float | None = None,
        **kwargs,
    ) -> None:
        """Plot bar

        Parameters
        ----------
        x : list[float] | np.ndarray
            Bar x coordinates
        height : list[float] | np.ndarray
            Bar heights
        width : float, optional
            Bar width
        bottom : float | np.ndarray, optional
            Bar bottom(s)
        align : str, optional
            Bar alignment type (`center` or `edge`)
        vmin : float, optional
            Y min value
        vmax : float | None, optional
            Y max value. If None, `np.max(height + bottom)` is set.
        **kwargs : dict, optional
            Axes.bar properties (e.g. `color="tomato", ec="black", lw=0.5, hatch="//"`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.bar.html>
        """
        # Check x, height list length
        if len(x) != len(height):
            err_msg = f"List length is not match ({len(x)=}, {len(height)=})"
            raise ValueError(err_msg)

        # Calculate top & vmax
        if isinstance(bottom, (list, tuple, np.ndarray)):
            bottom = np.array(bottom)
        else:
            bottom = np.array([bottom])
        top = np.array(height) + bottom
        vmax = max(top) if vmax is None else vmax

        # Check if bottom & top(height + bottom) is in valid min-max range
        self._check_value_min_max(bottom, vmin, vmax)
        self._check_value_min_max(top, vmin, vmax)

        # Calculate bar params
        rad = list(map(self.x_to_rad, x))
        r_bottom = np.array([self._y_to_r(v, vmin, vmax) for v in bottom])
        r_height = [self._y_to_r(v, vmin, vmax) for v in top] - r_bottom
        rad_width = self.rad_size * (width / self.size)

        def plot_bar(ax: PolarAxes) -> None:
            ax.bar(rad, r_height, rad_width, r_bottom, align=align, **kwargs)

        self._plot_funcs.append(plot_bar)

    def fill_between(
        self,
        x: list[float] | np.ndarray,
        y1: list[float] | np.ndarray,
        y2: float | list[float] | np.ndarray = 0,
        *,
        vmin: float = 0,
        vmax: float | None = None,
        **kwargs,
    ) -> None:
        """Fill the area between two horizontal(y1, y2) curves

        Parameters
        ----------
        x : list[float] | np.ndarray
            X coordinates
        y1 : list[float] | np.ndarray
            Y coordinates (first curve definition)
        y2 : float | list[float] | np.ndarray, optional
            Y coordinate[s] (second curve definition)
        vmin : float, optional
            Y min value
        vmax : float | None, optional
            Y max value. If None, `max(y1 + y2)` is set.
        **kwargs : dict, optional
            Axes.fill_between properties (e.g. `fc="red", ec="black", lw=0.1, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.fill_between.html>
        """
        rad = list(map(self.x_to_rad, x))
        if isinstance(y2, (list, tuple, np.ndarray)):
            y_all = list(y1) + list(y2)
        else:
            y_all = list(y1) + [y2]
            y2 = [y2] * len(x)
        vmin = min(y_all) if vmin is None else vmin
        vmax = max(y_all) if vmax is None else vmax
        self._check_value_min_max(y_all, vmin, vmax)

        r2 = [self._y_to_r(v, vmin, vmax) for v in y2]
        arc_rad, arc_r2 = self._to_arc_radr(rad, r2)
        r = [self._y_to_r(v, vmin, vmax) for v in y1]
        _, arc_r = self._to_arc_radr(rad, r)

        def plot_fill_between(ax: PolarAxes) -> None:
            ax.fill_between(arc_rad, arc_r, arc_r2, **kwargs)  # type: ignore

        self._plot_funcs.append(plot_fill_between)

    def heatmap(
        self,
        data: np.ndarray,
        *,
        vmin: float | None = None,
        vmax: float | None = None,
        start: float | None = None,
        end: float | None = None,
        cmap: str | Colormap = "bwr",
        show_value: bool = False,
        rect_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Plot heatmap

        Parameters
        ----------
        data : np.ndarray
            Numpy 2d array matrix values
        vmin : float | None, optional
            Min value for heatmap plot. If None, `np.min(data)` is set.
        vmax : float | None, optional
            Max value for heatmap plot. If None, `np.max(data)` is set.
        start : float | None, optional
            Start position for heatmap plot (x coordinate).
            If None, `track.start` is set.
        end : float | None, optional
            End position for heatmap plot (x coordinate).
            If None, `track.end` is set.
        cmap : str | Colormap, optional
            Colormap (e.g. `viridis`, `Spectral`, `Reds`, `Greys`)
            <https://matplotlib.org/stable/tutorials/colors/colormaps.html>
        show_value : bool, optional
            If True, show data value on heatmap rectangle
        rect_kws : dict[str, Any] | None, optional
            Patch properties (e.g. `dict(ec="black", lw=0.5, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        text_kws : dict[str, Any] | None, optional
            Text properties (e.g. `dict(size=6, color="red", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        rect_kws = {} if rect_kws is None else deepcopy(rect_kws)
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        # Set default value for None properties
        vmin = np.min(data) if vmin is None else vmin
        vmax = np.max(data) if vmax is None else vmax
        start = self.start if start is None else start
        end = self.end if end is None else end

        # Calculate radius & x position range list of heatmap rectangle
        data_row_num, data_col_num = data.shape
        unit_r_size = self.r_plot_size / data_row_num
        unit_x_size = (end - start) / data_col_num
        self._check_value_min_max(data, vmin, vmax)

        r_range_list: list[tuple[float, float]] = []
        for i in range(data_row_num):
            max_range = max(self.r_plot_lim) - (unit_r_size * i)
            min_range = max_range - unit_r_size
            r_range_list.append((min_range, max_range))
        x_range_list: list[tuple[float, float]] = []
        for i in range(data_col_num):
            min_range = start + (unit_x_size * i)
            max_range = min_range + unit_x_size
            x_range_list.append((min_range, max_range))

        # Plot heatmap
        colormap = cmap if isinstance(cmap, Colormap) else mpl.colormaps[cmap]
        norm = Normalize(vmin=vmin, vmax=vmax)
        for row_idx, row in enumerate(data):
            for col_idx, v in enumerate(row):
                # Plot heatmap rectangle
                rect_start, rect_end = x_range_list[col_idx]
                rect_r_lim = r_range_list[row_idx]
                color = colormap(norm(v))
                rect_kws.update(dict(fc=color, facecolor=color))
                self.rect(rect_start, rect_end, r_lim=rect_r_lim, **rect_kws)

                if show_value:
                    # Plot value text on heatmap rectangle
                    text_value = f"{v:.2f}" if isinstance(v, float) else str(v)
                    text_x = (rect_end + rect_start) / 2
                    text_r = sum(rect_r_lim) / 2
                    self.text(text_value, text_x, text_r, **text_kws)

    def tree(
        self,
        treedata: str | Path | StringIO | Tree,
        *,
        format: str = "newick",
        outer: bool = True,
        use_branch_length: bool = False,
        leaf_label_size: float = 0,
        innode_label_size: float = 0,
        leaf_label_margin: float = 0.5,
        label_formatter: Callable[[Clade], str] | None = None,
        node_color_list: list[tuple[list[str], str]] | None = None,
        line_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Plot tree

        It is recommended that the track(sector) size be the same as the number of
        leaf nodes in the tree, to make it easier to combine with `bar` and `heatmap`.

        Parameters
        ----------
        treedata : str | Path | StringIO | Tree
            Phylogenetic tree data (`File-like object` or `Tree object`)
        format : str, optional
            Tree format (e.g. `newick`, `nexus`, `phyloxml`, ...)
        outer : bool, optional
            If True, plot tree on outer side. If False, plot tree on inner side.
        use_branch_length : bool, optional
            If True, tree branch length is used for plot tree.
            If False, plot tree with ultrametric tree style.
        leaf_label_size : float, optional
            Leaf node label size. By default, `size=0` (No display).
        innode_label_size : float, optional
            Internal node label size. By default, `size=0` (No display).
        leaf_label_margin : float, optional
            leaf node label margin size (Radius unit)
        label_formatter : Callable[[Clade], str] | None, optional
            User-defined label format function.
            (e.g. `lambda node: node.name.replace("_", " ")`)
        node_color_list : list[tuple[list[str], str]] | None, optional
            Tree node & color setting list.
            If multi nodes are set, MRCA(Most Recent Common Ancestor) node of
            target nodes is set.
            (e.g. `[(["node1"], "red"), (["taxa1", "taxa2"], "blue"), ...]`)
        line_kws : dict[str, Any] | None, optional
            Patch properties (e.g. `dict(ec="red", lw=1, ls="dashed", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        text_kws : dict[str, Any] | None, optional
            Text properties (e.g. `dict(color="red", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        node_color_list = [] if node_color_list is None else deepcopy(node_color_list)
        line_kws = {} if line_kws is None else deepcopy(line_kws)
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        # Load tree data, set node names, set node colors
        tree = utils.TreeUtil.load_tree(treedata, format)
        tree = utils.TreeUtil.set_unique_node_name(tree)
        tree = utils.TreeUtil.set_node_color(tree, node_color_list)
        utils.TreeUtil.check_node_name_dup(tree)
        # If not `use_branch_length` or not branch length exists
        # Convert tree to ultrametric tree
        name2depth = {n.name: d for n, d in tree.depths().items()}
        max_tree_depth = max(name2depth.values())
        if not use_branch_length or max_tree_depth == 0:
            tree = utils.TreeUtil.to_ultrametric_tree(tree)
            name2depth = {n.name: d for n, d in tree.depths().items()}
            max_tree_depth = max(name2depth.values())
        # Calculate x, r unit size of depth
        x_unit_size = self.size / tree.count_terminals()
        r_unit_size = self.r_plot_size / max_tree_depth
        # Calculate leaf node (x, r) coordinate
        name2xr: dict[Any, tuple[float, float]] = {}
        node: Clade
        for idx, node in enumerate(tree.get_terminals()):
            x = self.start + (x_unit_size * idx) + (x_unit_size / 2)
            if outer:
                r = min(self.r_plot_lim) + r_unit_size * name2depth[node.name]
            else:
                r = max(self.r_plot_lim) - r_unit_size * name2depth[node.name]
            name2xr[node.name] = (x, r)
        # Calculate internal node (x, r) coordinate
        for node in tree.get_nonterminals(order="postorder"):
            x = sum([name2xr[n.name][0] for n in node.clades]) / len(node.clades)
            if outer:
                r = min(self.r_plot_lim) + r_unit_size * name2depth[node.name]
            else:
                r = max(self.r_plot_lim) - r_unit_size * name2depth[node.name]
            name2xr[node.name] = (x, r)
        # Plot tree by node (x, r) coordinate
        for node in tree.get_nonterminals():
            parent_x, parent_r = name2xr[node.name]
            child_node: Clade
            for child_node in node.clades:
                child_x, child_r = name2xr[child_node.name]
                # Set node color if exists
                _line_kws = deepcopy(line_kws)
                if child_node.color is not None:
                    _line_kws.update(dict(color=child_node.color))
                # Plot horizontal line
                h_line_points = (parent_x, child_x), (parent_r, parent_r)
                self._simpleline(*h_line_points, **_line_kws)
                # Plot vertical line
                v_line_points = (child_x, child_x), (parent_r, child_r)
                self._simpleline(*v_line_points, **_line_kws)
        # Plot label if required
        for node in tree.find_clades():
            label = str(node.name) if label_formatter is None else label_formatter(node)
            x, r = name2xr[node.name]
            _text_kws = deepcopy(text_kws)
            if node.is_terminal() and leaf_label_size > 0:
                orientation = "vertical"
                r = r + leaf_label_margin if outer else r - leaf_label_margin
                _text_kws.update(dict(size=leaf_label_size))
                if node.color is not None:
                    _text_kws.update(dict(color=node.color))
            elif not node.is_terminal() and innode_label_size > 0:
                orientation = "horizontal"
                _text_kws.update(dict(size=innode_label_size))
            else:
                continue
            rad = self.x_to_rad(x)
            params = utils.plot.get_label_params_by_rad(rad, orientation, outer)
            _text_kws.update(params)
            self.text(label, x, r, orientation=orientation, **_text_kws)

    def genomic_features(
        self,
        features: list[SeqFeature],
        *,
        plotstyle: str = "box",
        r_lim: tuple[float, float] | None = None,
        facecolor_handler: Callable[[SeqFeature], str] | None = None,
        **kwargs,
    ) -> None:
        """Plot genomic features

        Parameters
        ----------
        features : list[SeqFeature]
            Biopython's SeqFeature list
        plotstyle : str, optional
            Plot style (`box` or `arrow`)
        r_lim : tuple[float, float] | None, optional
            Radius limit range. If None, `track.r_plot_lim` is set.
        facecolor_handler : Callable[[SeqFeature], str] | None, optional
            User-defined function to handle facecolor
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="blue", lw=1.0, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        if r_lim is None:
            r_lim = self.r_plot_lim
        else:
            if not min(self.r_lim) <= min(r_lim) < max(r_lim) <= max(self.r_lim):
                raise ValueError(f"{r_lim=} is invalid track range.\n{self}")

        for feature in features:
            # Set qualifier tag facecolor if exists
            tag_color = feature.qualifiers.get("facecolor", [None])[0]
            if tag_color is not None:
                kwargs.update(dict(fc=tag_color, facecolor=tag_color))
            # Set facecolor by user-defined function
            if facecolor_handler is not None:
                color = facecolor_handler(feature)
                kwargs.update(dict(fc=color, facecolor=color))
            # Plot feature
            try:
                start = int(str(feature.location.parts[0].start))
                end = int(str(feature.location.parts[-1].end))
            except ValueError:
                print(f"Failed to parse feature's start-end position.\n{feature}")
                continue
            if feature.strand == -1:
                start, end = end, start
            if plotstyle == "box":
                self.rect(start, end, r_lim=r_lim, **kwargs)
            elif plotstyle == "arrow":
                self.arrow(start, end, r_lim=r_lim, **kwargs)
            else:
                raise ValueError(f"{plotstyle=} is invalid ('box' or 'arrow').")

    ############################################################
    # Private Method
    ############################################################

    def _y_to_r(self, y: float, vmin: float, vmax: float) -> float:
        """Convert y coordinate to radius in track

        Parameters
        ----------
        y : float
            Y coordinate
        vmin : float
            Min y coordinate
        vmax : float
            Max y coordinate

        Returns
        -------
        r : float
            Converted radius position
        """
        norm = Normalize(vmin, vmax)
        r = min(self.r_plot_lim) + (self.r_plot_size * norm(y))
        return r

    def _to_arc_radr(
        self,
        rad: list[float] | np.ndarray,
        r: list[float] | np.ndarray,
    ) -> tuple[list[float], list[float]]:
        """Convert radian & radius to arc radian & arc radius

        Parameters
        ----------
        rad : list[float] | np.ndarray
            Radian list
        r : list[float] | np.ndarray
            Radius list

        Returns
        -------
        arc_rad, arc_r : tuple[list[float], list[float]]
            Arc radian list, Ard radius list
        """
        all_arc_rad, all_arc_r = [], []
        for i in range(len(rad) - 1):
            rad1, rad2, r1, r2 = rad[i], rad[i + 1], r[i], r[i + 1]
            if rad1 == rad2:
                all_arc_rad.extend([rad1, rad2])
                all_arc_r.extend([r1, r2])
            else:
                step = config.ARC_RADIAN_STEP
                if rad1 > rad2:
                    step *= -1
                arc_rad = list(np.arange(rad1, rad2, step)) + [rad2]
                all_arc_rad.extend(arc_rad)
                arc_r = np.linspace(r1, r2, len(arc_rad), endpoint=True)
                all_arc_r.extend(arc_r)
        return all_arc_rad, all_arc_r

    def _simpleline(
        self,
        x_lim: tuple[float, float],
        r_lim: tuple[float, float],
        **kwargs,
    ) -> None:
        """Plot simple patch line between two points (x1, r1), (x2, r2)

        Used to plot simple lines such as ticks internally

        Parameters
        ----------
        x_lim : tuple[float, float]
            X start-end limit region
        r_lim : tuple[float, float]
            Radius start-end limit region
        **kwargs : dict, optional
            Patch properties (e.g. `ec="red", lw=1.0, ...`)
            https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html
        """
        rad_lim = tuple(map(self.x_to_rad, x_lim, (True, True)))
        self._patches.append(ArcLine(rad_lim, r_lim, **kwargs))

    def _check_value_min_max(
        self,
        value: float | list[int] | list[float] | np.ndarray,
        vmin: float,
        vmax: float,
    ) -> None:
        """Check if value(s) is in valid min-max range

        Parameters
        ----------
        value : float | list[int] | list[float] | np.ndarray
            Check value(s)
        vmin : float
            Min value
        vmax : float
            Max value
        """
        if isinstance(value, (list, tuple, np.ndarray)):
            if isinstance(value, np.ndarray):
                value = list(value.flatten())
            for v in value:
                if not vmin <= v <= vmax:
                    err_msg = f"value={v} is not in valid range ({vmin=}, {vmax=})"
                    raise ValueError(err_msg)
        else:
            if not vmin <= value <= vmax:
                err_msg = f"{value=} is not in valid range ({vmin=}, {vmax=})"
                raise ValueError(err_msg)

    def __str__(self):
        return (
            f"# Track = '{self.name}' (Parent Sector = '{self.parent_sector.name}')\n"
            f"# Size = {self.size} ({self.start} - {self.end})\n"
            f"# Radius size = {self.r_size:.2f} "
            f"({min(self.r_lim):.2f} - {max(self.r_lim):.2f})\n"
            f"# Degree size = {self.deg_size:.2f} "
            f"({min(self.deg_lim):.2f} - {max(self.deg_lim):.2f})\n"
        )
