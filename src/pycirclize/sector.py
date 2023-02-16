from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Callable
from urllib.parse import urlparse
from urllib.request import urlopen

import numpy as np
from matplotlib.patches import Patch
from matplotlib.projections.polar import PolarAxes
from PIL import Image

from pycirclize import config, utils
from pycirclize.patches import ArcLine, ArcRectangle
from pycirclize.track import Track


class Sector:
    """Circos Sector Class"""

    def __init__(
        self,
        name: str,
        size: float,
        rad_lim: tuple[float, float],
        start_pos: float = 0,
    ):
        """
        Parameters
        ----------
        name : str
            Sector name
        size : float
            Sector size
        rad_lim : tuple[float, float]
            Sector radian limit region
        start_pos : float, optional
            Sector start position
        """
        self._name = name
        self._size = size
        self._rad_lim = rad_lim
        self._start_pos = start_pos
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
        return self._start_pos

    @property
    def end(self) -> float:
        """Sector end position (x coordinate)"""
        return self._start_pos + self._size

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
        return tuple(map(math.degrees, self.rad_lim))

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
        if not self.start <= x <= self.end and not ignore_range_error:
            err_msg = f"{x=} is invalid range of '{self.name}' sector.\n{self}"
            raise ValueError(err_msg)
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
        orientation: str = "horizontal",
        **kwargs,
    ) -> None:
        """Plot text

        Parameters
        ----------
        text : str
            Text content
        x: float, optional
            X position. If None, sector center x is set.
        r : float, optional
            Radius position. By default, outer position `r=105` is set.
        orientation : str, optional
            Text orientation (`horizontal` or `vertical`)
        **kwargs : dict, optional
            Text properties (e.g. `size=12, color="red", va="center", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        if x is None:
            # Get sector center radian position
            center_x = (self.start + self.end) / 2
            rad = self.x_to_rad(center_x)
        else:
            rad = self.x_to_rad(x)

        # Set label proper alignment, rotation parameters by radian
        params = utils.plot.get_label_params_by_rad(rad, orientation, outer=True)
        kwargs.update(params)

        def plot_text(ax: PolarAxes) -> None:
            ax.text(rad, r, text, **kwargs)

        self._plot_funcs.append(plot_text)

    def line(
        self,
        r: float,
        start: float | None = None,
        end: float | None = None,
        **kwargs,
    ) -> None:
        """Plot line

        Parameters
        ----------
        r : float
            Line radius position (0 - 100)
        start : float, optional
            Start position (x coordinate). If None, `sector.start` is set.
        end : float, optional
            End position (x coordinate). If None, `sector.end` is set.
        **kwargs : dict, optional
            Patch properties (e.g. `color="red", lw=3, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        start = self.start if start is None else start
        end = self.end if end is None else end
        rad_lim = (self.x_to_rad(start), self.x_to_rad(end))
        self._patches.append(ArcLine(rad_lim, (r, r), **kwargs))

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
        start : float, optional
            Start position (x coordinate). If None, `sector.start` is set.
        end : float, optional
            End position (x coordinate). If None, `sector.end` is set.
        r_lim : tuple[float, float] | None
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
        img: str | Path | np.ndarray | Image.Image,
        *,
        size: float = 0.05,
        x: float | None = None,
        r: float = 105,
        label: str | None = None,
        label_pos: str = "bottom",
        label_margin: float = 0.1,
        show_border: bool = True,
        imshow_kws: dict[str, Any] = {},
        text_kws: dict[str, Any] = {},
    ) -> None:
        """Plot raster image

        This method is experimental. API may change in the future release.

        Parameters
        ----------
        img : str | Path | np.ndarray | Image.Image
            Image for plotting (`File Path`|`URL`|`Numpy Array`|`PIL Image`)
        size : float, optional
            Image size (ratio to overall figure size)
        x : float | None, optional
            X position. If None, sector center x position is set.
        r : float, optional
            Radius position
        label : str | None, optional
            Image label. If None, no plotting label.
        label_pos : str, optional
            Label plot position (`bottom` or `top`)
        label_margin : float, optional
            Label margin
        show_border : bool, optional
            If True, show label border.
        imshow_kws : dict[str, Any], optional
            Axes.imshow properties
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html>
        text_kws : dict[str, Any], optional
            Text properties (e.g. `dict(size=10, color="red", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        # Load image data
        if isinstance(img, str) and urlparse(img).scheme in ("http", "https"):
            im = Image.open(urlopen(img))
        elif isinstance(img, (str, Path)):
            im = Image.open(str(img))
        else:
            im = img

        # Calculate x, y image set position
        x = (self.start + self.end) / 2 if x is None else x
        rad = self.x_to_rad(x)
        max_r_lim = config.MAX_R + config.R_PLOT_MARGIN
        im_x = np.cos((np.pi / 2) - rad) * (r / max_r_lim)
        im_y = np.sin((np.pi / 2) - rad) * (r / max_r_lim)
        # Normalize (-1, 1) to (0, 1) axis range
        im_x = (im_x + 1) / 2
        im_y = (im_y + 1) / 2

        def plot_raster(ax: PolarAxes) -> None:
            # Set inset axes
            bounds = [im_x - (size / 2), im_y - (size / 2), size, size]
            axin = ax.inset_axes(bounds, transform=ax.transAxes)
            # Set image border param
            if show_border:
                axin.tick_params(bottom=False, left=False)
                axin.tick_params(labelbottom=False, labelleft=False)
            else:
                axin.axis("off")
            # Plot raster image
            axin.imshow(im, **imshow_kws)

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
        return (
            f"# Sector = '{self.name}'\n"
            f"# Size = {self.size} ({self.start} - {self.end})\n"
            f"# Radian size = {self.rad_size:.2f} "
            f"({min(self.rad_lim):.2f} - {max(self.rad_lim):.2f})\n"
            f"# Degree size = {self.deg_size:.2f} "
            f"({min(self.deg_lim):.2f} - {max(self.deg_lim):.2f})\n"
        )
