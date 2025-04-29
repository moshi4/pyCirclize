from __future__ import annotations

import math
import textwrap
import warnings
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable

import numpy as np
from matplotlib.patches import Patch
from matplotlib.projections.polar import PolarAxes
from PIL import Image, ImageOps

from pycirclize import config, utils
from pycirclize.patches import ArcLine, ArcRectangle, Line
from pycirclize.track import Track
from pycirclize.utils.plot import get_label_params_by_rad


class Sector:
    """Circos Sector Class"""

    def __init__(
        self,
        name: str,
        size: float | tuple[float, float],
        rad_lim: tuple[float, float],
        clockwise: bool = True,
    ):
        """
        Parameters
        ----------
        name : str
            Sector name
        size : float | tuple[float, float]
            Sector size (or range)
        rad_lim : tuple[float, float]
            Sector radian limit region
        clockwise : bool, optional
            Sector coordinate direction (clockwise or anti-clockwise).
        """
        self._name = name
        if isinstance(size, (tuple, list)):
            start, end = size[0], size[1]
        else:
            start, end = 0, size
        self._start = start
        self._end = end
        self._size = end - start
        self._rad_lim = rad_lim
        self._clockwise = clockwise
        self._tracks: list[Track] = []

        # Plot data and functions
        self._patches: list[Patch] = []
        self._plot_funcs: list[Callable[[PolarAxes], None]] = []

    ############################################################
    # Property
    ############################################################

    @property
    def name(self) -> str:
        """Sector name"""
        return self._name

    @property
    def size(self) -> float:
        """Sector size (x coordinate)"""
        return self._size

    @property
    def start(self) -> float:
        """Sector start position (x coordinate)"""
        return self._start

    @property
    def end(self) -> float:
        """Sector end position (x coordinate)"""
        return self._end

    @property
    def center(self) -> float:
        """Sector center position (x coordinate)"""
        return (self.start + self.end) / 2

    @property
    def rad_size(self) -> float:
        """Sector radian size"""
        return max(self.rad_lim) - min(self.rad_lim)

    @property
    def rad_lim(self) -> tuple[float, float]:
        """Sector radian limit"""
        return self._rad_lim

    @property
    def deg_size(self) -> float:
        """Sector degree size"""
        return max(self.deg_lim) - min(self.deg_lim)

    @property
    def deg_lim(self) -> tuple[float, float]:
        """Sector degree limit"""
        return (math.degrees(self.rad_lim[0]), math.degrees(self.rad_lim[1]))

    @property
    def clockwise(self) -> bool:
        """Sector coordinate direction"""
        return self._clockwise

    @property
    def tracks(self) -> list[Track]:
        """Tracks in sector"""
        return self._tracks

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

    def add_track(
        self,
        r_lim: tuple[float, float],
        *,
        r_pad_ratio: float = 0,
        name: str | None = None,
    ) -> Track:
        """Add track to sector

        Parameters
        ----------
        r_lim : tuple[float, float]
            Radius limit region (0 - 100)
        r_pad_ratio : float
            Radius padding ratio for plot data (0 - 1.0)
        name : str | None, optional
            Track name. If None, `Track{track_idx}` is set.

        Returns
        -------
        track : Track
            Track
        """
        name = f"Track{len(self.tracks) + 1:02d}" if name is None else name
        if name in [t.name for t in self.tracks]:
            raise ValueError(f"{name=} track is already exists.")
        if not 0 <= min(r_lim) <= max(r_lim) <= 100:
            warnings.warn(f"{r_lim=} is unexpected plot range (0 <= r <= 100).")
        track = Track(name, r_lim, r_pad_ratio, self)
        self._tracks.append(track)
        return track

    def get_track(self, name: str) -> Track:
        """Get track by name

        Parameters
        ----------
        name : str
            Track name

        Returns
        -------
        track : Track
            Target name track
        """
        name2track = {t.name: t for t in self.tracks}
        if name not in name2track:
            raise ValueError(f"{name=} track not exists.")
        return name2track[name]

    def get_lowest_r(self) -> float:
        """Get lowest radius position of sector from tracks data

        Returns
        -------
        lowest_r : float
            Lowest radius position. If no tracks found, `lowest_r=100`.
        """
        if len(self.tracks) == 0:
            return config.MAX_R
        return min([min(t.r_lim) for t in self.tracks])

    def x_to_rad(self, x: float, ignore_range_error: bool = False) -> float:
        """Convert x coordinate to radian in sector start-end range

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
        # Check target x is in valid sector range
        if not ignore_range_error:
            # Apply relative torelance value to sector range to avoid
            # unexpected invalid range error due to rounding errors (Issue #27, #67)
            min_range = self.start - config.EPSILON
            max_range = self.end + config.EPSILON
            if not min_range <= x <= max_range:
                raise ValueError(f"{x=} is invalid range of '{self.name}' sector.\n{self}")  # fmt: skip  # noqa: E501

        if not self.clockwise:
            x = (self.start + self.end) - x
        size_ratio = self.rad_size / self.size
        x_from_start = x - self.start
        rad_from_start = x_from_start * size_ratio
        rad = min(self.rad_lim) + rad_from_start
        return rad

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
        self.rect(self.start, self.end, config.R_LIM, **fc_behind_kwargs)

        # Axis edgecolor placed in front of other patches (zorder=1.01)
        ec_front_kwargs = {**kwargs, **config.AXIS_EDGE_PARAM}
        self.rect(self.start, self.end, config.R_LIM, **ec_front_kwargs)

    def text(
        self,
        text: str,
        x: float | None = None,
        r: float = 105,
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
        x: float | None, optional
            X position. If None, sector center x is set.
        r : float, optional
            Radius position. By default, outer position `r=105` is set.
        adjust_rotation : bool, optional
            If True, text rotation is auto set based on `x` and `orientation` params.
        orientation : str, optional
            Text orientation (`horizontal` or `vertical`)
        ignore_range_error : bool, optional
            If True, ignore x position range error
            (ErrorCase: `not track.start <= x <= track.end`)
        **kwargs : dict, optional
            Text properties (e.g. `size=12, color="red", va="center", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        # If value is None, center position is set.
        x = self.center if x is None else x
        rad = self.x_to_rad(x, ignore_range_error)

        if adjust_rotation:
            # Set label proper alignment, rotation parameters by radian
            params = utils.plot.get_label_params_by_rad(rad, orientation)
            kwargs.update(params)

        if "ha" not in kwargs and "horizontalalignment" not in kwargs:
            kwargs.update(dict(ha="center"))
        if "va" not in kwargs and "verticalalignment" not in kwargs:
            kwargs.update(dict(va="center"))

        def plot_text(ax: PolarAxes) -> None:
            ax.text(rad, r, text, **kwargs)

        self._plot_funcs.append(plot_text)

    def line(
        self,
        *,
        r: float | tuple[float, float],
        start: float | None = None,
        end: float | None = None,
        arc: bool = True,
        **kwargs,
    ) -> None:
        """Plot line

        Parameters
        ----------
        r : float | tuple[float, float]
            Line radius position (0 - 100). If r is float, (r, r) is set.
        start : float | None, optional
            Start position (x coordinate). If None, `sector.start` is set.
        end : float | None, optional
            End position (x coordinate). If None, `sector.end` is set.
        arc : bool, optional
            If True, plot arc style line for polar projection.
            If False, simply plot linear style line.
        **kwargs : dict, optional
            Patch properties (e.g. `color="red", lw=3, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        start = self.start if start is None else start
        end = self.end if end is None else end
        rad_lim = (self.x_to_rad(start), self.x_to_rad(end))
        r_lim = r if isinstance(r, (tuple, list)) else (r, r)
        LinePatch = ArcLine if arc else Line
        self._patches.append(LinePatch(rad_lim, r_lim, **kwargs))

    def rect(
        self,
        start: float | None = None,
        end: float | None = None,
        r_lim: tuple[float, float] | None = None,
        **kwargs,
    ) -> None:
        """Plot rectangle

        Parameters
        ----------
        start : float | None, optional
            Start position (x coordinate). If None, `sector.start` is set.
        end : float | None, optional
            End position (x coordinate). If None, `sector.end` is set.
        r_lim : tuple[float, float] | None, optional
            Radius limit region. If None, (0, 100) is set.
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="blue", lw=1.0, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        start = self.start if start is None else start
        end = self.end if end is None else end
        rad_rect_start = self.x_to_rad(start)
        rad_rect_end = self.x_to_rad(end)

        r_lim = config.R_LIM if r_lim is None else r_lim
        min_rad = min(rad_rect_start, rad_rect_end)
        max_rad = max(rad_rect_start, rad_rect_end)

        radr = (min_rad, min(r_lim))
        width = max_rad - min_rad
        height = max(r_lim) - min(r_lim)
        self._patches.append(ArcRectangle(radr, width, height, **kwargs))

    def raster(
        self,
        img: str | Path | Image.Image,
        *,
        size: float = 0.05,
        x: float | None = None,
        r: float = 105,
        rotation: int | float | str | None = None,
        border_width: int = 0,
        label: str | None = None,
        label_pos: str = "bottom",
        label_margin: float = 0.1,
        imshow_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Plot raster image

        This method is experimental. API may change in the future release.

        Parameters
        ----------
        img : str | Path | Image
            Image data (`File Path`|`URL`|`PIL Image`)
        size : float, optional
            Image size (ratio to overall figure size)
        x : float | None, optional
            X position. If None, sector center x position is set.
        r : float, optional
            Radius position
        rotation : int | float | str | None, optional
            Image rotation setting.
            If `None`, no rotate image (default).
            If `auto`, rotate image by auto set rotation.
            If `int` or `float` value, rotate image by user-specified value.
        border_width : int, optional
            Border width in pixel. By default, 0 is set (no border shown).
        label : str | None, optional
            Image label. If None, no label shown.
        label_pos : str, optional
            Label plot position (`bottom` or `top`)
        label_margin : float, optional
            Label margin
        imshow_kws : dict[str, Any] | None, optional
            Axes.imshow properties
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html>
        text_kws : dict[str, Any] | None, optional
            Text properties (e.g. `dict(size=10, color="red", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        imshow_kws = {} if imshow_kws is None else deepcopy(imshow_kws)
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        # Load image data
        im = utils.load_image(img)

        # Draw border on image
        if border_width > 0:
            im = ImageOps.expand(im, border=border_width, fill="black")

        # Rotate image
        x = self.center if x is None else x
        rad = self.x_to_rad(x)
        if isinstance(rotation, (int, float)):
            im = im.rotate(rotation, expand=True)
            rotate_value = rotation
        elif rotation == "auto":
            rotate_value: float = get_label_params_by_rad(rad, "horizontal")["rotation"]
            im = im.rotate(rotate_value, expand=True)
        elif rotation is None:
            rotate_value = 0
        else:
            raise ValueError(f"{rotation=} is invalid.")

        # Calculate x, y image set position
        max_r_lim = config.MAX_R + config.R_PLOT_MARGIN
        im_x: float = np.cos((np.pi / 2) - rad) * (r / max_r_lim)
        im_y: float = np.sin((np.pi / 2) - rad) * (r / max_r_lim)
        # Normalize (-1, 1) to (0, 1) axis range
        im_x = (im_x + 1) / 2
        im_y = (im_y + 1) / 2

        # TODO: Terrible code to be fixed in the future
        # Approximate image size calculation logic, not complete
        scale = 1 - (abs(abs(rotate_value) % 90 - 45) / 45)  # 0 - 1.0
        size_ratio = 1 + (scale * (np.sqrt(2) - 1))
        size = size * size_ratio

        def plot_raster(ax: PolarAxes) -> None:
            # Set inset axes & plot raster image
            bounds = (im_x - (size / 2), im_y - (size / 2), size, size)
            axin = ax.inset_axes(bounds, transform=ax.transAxes)
            axin.axis("off")
            axin.imshow(im, **imshow_kws)  # type: ignore

            # Plot label
            if label is not None:
                text_x = sum(axin.get_xlim()) / 2
                y_size = max(axin.get_ylim()) - min(axin.get_ylim())
                if label_pos == "bottom":
                    text_y = max(axin.get_ylim()) + (y_size * label_margin)
                    va = "top"
                elif label_pos == "top":
                    text_y = min(axin.get_ylim()) - (y_size * label_margin)
                    va = "bottom"
                else:
                    raise ValueError(f"{label_pos=} is invalid ('top' or 'bottom').")
                axin.text(text_x, text_y, label, ha="center", va=va, **text_kws)

        self._plot_funcs.append(plot_raster)

    ############################################################
    # Private Method
    ############################################################

    def __str__(self):
        min_deg_lim, max_deg_lim = min(self.deg_lim), max(self.deg_lim)
        track_names = [t.name for t in self.tracks]
        return textwrap.dedent(
            f"""
            # Sector = '{self.name}'
            # Size = {self.size} ({self.start} - {self.end})
            # Degree Size = {self.deg_size:.2f} ({min_deg_lim:.2f} - {max_deg_lim:.2f})
            # Track List = {track_names}
            """
        )[1:]
