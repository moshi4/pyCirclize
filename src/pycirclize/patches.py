from __future__ import annotations

import math

import numpy as np
from matplotlib.patches import PathPatch
from matplotlib.path import Path

from pycirclize import config


class ArcLine(PathPatch):
    """Arc Line Patch"""

    def __init__(
        self,
        rad_lim: tuple[float, float],
        r_lim: tuple[float, float],
        **kwargs,
    ):
        """
        Parameters
        ----------
        rad_lim : tuple[float, float]
            Radian limit region
        r_lim : tuple[float, float]
            Radius limit region
        **kwargs : dict, optional
            Patch properties (e.g. `ec="red", lw=1.0, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Default params: fc='none', color='black', linewidth=0.5
        kwargs.update(dict(fc="none"))
        if "ec" not in kwargs and "edgecolor" not in kwargs and "color" not in kwargs:
            kwargs.update(dict(ec="black"))
        if "lw" not in kwargs and "linewidth" not in kwargs:
            kwargs.update(dict(lw=0.5))

        # Calculate line path vertices
        if rad_lim[1] >= rad_lim[0]:
            rad_start, rad_end = rad_lim
            r_start, r_end = r_lim
        else:
            rad_start, rad_end = rad_lim[::-1]
            r_start, r_end = r_lim[::-1]
        if rad_start == rad_end:
            arc_rads = [rad_start, rad_end]
        else:
            step = config.ARC_RADIAN_STEP
            arc_rads = list(np.arange(rad_start, rad_end, step)) + [rad_end]
        arc_r_list = np.linspace(r_start, r_end, len(arc_rads), endpoint=True)

        # Set line path
        verts = list(zip(arc_rads, arc_r_list))
        super().__init__(Path(verts), **kwargs)


class ArcRectangle(PathPatch):
    """Arc Rectangle PathPatch"""

    def __init__(
        self,
        radr: tuple[float, float],
        width: float,
        height: float,
        **kwargs,
    ):
        """
        Parameters
        ----------
        radr : tuple[float, float]
            Anchor point (rad=`radian`, r=`radius`)
        width : float
            Rectangle radian width
        height : float
            Rectangle radius height
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="blue", lw=2.0, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        min_rad, min_r = radr
        max_rad, max_r = min_rad + width, min_r + height
        arc_rads = np.arange(min_rad, max_rad, config.ARC_RADIAN_STEP)
        arc_rads = np.append(arc_rads, max_rad)
        bottom_arc_path = list(zip(arc_rads, [min_r] * len(arc_rads)))
        upper_arc_path = list(zip(arc_rads[::-1], [max_r] * len(arc_rads)))
        arc_rect_path = Path(
            bottom_arc_path + upper_arc_path + [bottom_arc_path[0]],
            closed=True,
        )
        super().__init__(arc_rect_path, **kwargs)


class ArcArrow(PathPatch):
    """Arc Arrow PathPatch"""

    def __init__(
        self,
        rad: float,
        r: float,
        drad: float,
        dr: float,
        head_length: float = np.pi / 90,
        shaft_ratio: float = 0.5,
        **kwargs,
    ):
        """
        Parameters
        ----------
        rad : float
            Radian base coordinate
        r : float
            Radius base coordinate
        drad : float
            Radian size
        dr : float
            Radius size
        head_length : float, optional
            Arrow head length (Radian unit)
        shaft_ratio : float, optional
            Arrow shaft ratio (0 - 1.0)
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="blue", lw=1.0, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Set position parameters
        shaft_size = dr * shaft_ratio
        y_shaft_bottom = r + ((dr - shaft_size) / 2)
        y_shaft_upper = r + dr - ((dr - shaft_size) / 2)

        is_forward = True if drad >= 0 else False
        drad = abs(drad)
        if head_length > drad:
            head_length = drad
        if is_forward:
            rad_shaft_tip = rad + (drad - head_length)
            rad_arrow_tip = rad + drad
        else:
            rad_shaft_tip = rad - (drad - head_length)
            rad_arrow_tip = rad - drad

        # ArcArrow vertex points
        p1 = rad, y_shaft_bottom
        p2 = rad_shaft_tip, y_shaft_bottom
        p3 = rad_shaft_tip, r  # Arrow bottom tip point
        p4 = rad_arrow_tip, (r + (r + dr)) / 2  # Arrow center tip point
        p5 = rad_shaft_tip, r + dr  # Arrow upper tip point
        p6 = rad_shaft_tip, y_shaft_upper
        p7 = rad, y_shaft_upper

        # Create ArcArrow Path from vertex points
        step = config.ARC_RADIAN_STEP if is_forward else -config.ARC_RADIAN_STEP
        shaft_arc_rads = np.arange(p1[0], p2[0], step)
        bottom_shaft_r_list = [p1[1]] * len(shaft_arc_rads)
        upper_shaft_r_list = [p7[1]] * len(shaft_arc_rads)
        bottom_shaft_arc_path = list(zip(shaft_arc_rads, bottom_shaft_r_list))
        upper_shaft_arc_path = list(zip(shaft_arc_rads[::-1], upper_shaft_r_list))
        arc_arrow_path = Path(
            bottom_shaft_arc_path + [p2, p3, p4, p5, p6] + upper_shaft_arc_path + [p1],
            closed=True,
        )
        super().__init__(arc_arrow_path, **kwargs)


class BezierCurve(PathPatch):
    """Bezier Curve PathPatch"""

    def __init__(
        self,
        rad_start1: float,
        rad_end1: float,
        r1: float,
        rad_start2: float,
        rad_end2: float,
        r2: float,
        height_ratio: float = 0.5,
        direction: int = 0,
        arrow_length_ratio: float = 0.05,
        **kwargs,
    ):
        """
        Parameters
        ----------
        rad_start1 : float
            Radian start1
        rad_end1 : float
            Radian end1
        r1 : float
            Radius position1
        rad_start2 : float
            Radian start2
        rad_end2 : float
            Radian end2
        r2 : float
            Radius position2
        height_ratio : float, optional
            Bezier curve height ratio parameter
        direction : int, optional
            `0`: Circular edge shape (Default)
            `1`: Directional(1 -> 2) arrow edge shape
            `-1`: Directional(1 <- 2) arrow edge shape
            `2`: Bidirectional arrow edge shape
        arrow_length_ratio : float, optional
            Arrow length ratio.
        **kwargs : dict, optional
            Patch properties (e.g. `lw=1.0, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """

        def make_arc_paths(rad1: float, rad2: float, r: float):
            # If rad1 == rad2, return blank list
            arc_paths = []
            step = config.ARC_RADIAN_STEP if rad1 <= rad2 else -config.ARC_RADIAN_STEP
            for rad in np.arange(rad1, rad2, step):
                arc_paths.append((Path.LINETO, (rad, r)))
            return arc_paths

        def make_arrow_paths(rad1: float, rad2: float, r_side: float, r_top: float):
            return [
                (Path.LINETO, (rad1, r_side)),
                (Path.LINETO, ((rad1 + rad2) / 2, r_top)),
                (Path.LINETO, (rad2, r_side)),
            ]

        def make_bezier_paths(
            rad1: float, rad2: float, r1: float, r2: float, height_ratio: float = 0.5
        ):
            if height_ratio >= 0.5:
                # Case1: height_ratio: 0.50 => r_ctl_pos: 0
                # Case2: height_ratio: 0.75 => r_ctl_pos: 25
                # Case3: height_ratio: 1.00 => r_ctl_pos: 50
                r_ctl_pos = config.MAX_R * (height_ratio - 0.5)
                rad_ctl_pos = (rad1 + rad2) / 2 + math.pi
            else:
                # Case1: height_ratio: 0.25 => r_ctl_pos: 25
                # Case2: height_ratio: 0.00 => r_ctl_pos: 50
                r_ctl_pos = config.MAX_R * (0.5 - height_ratio)
                rad_ctl_pos = (rad1 + rad2) / 2
            return [
                (Path.LINETO, (rad1, r1)),
                (Path.CURVE3, (rad_ctl_pos, r_ctl_pos)),
                (Path.LINETO, (rad2, r2)),
            ]

        # Circos style plot order `start1 -> end1 -> end2 -> start2 -> start1`
        # http://circos.ca/documentation/tutorials/links/twists/images
        arrow_r1 = r1 * (1 - arrow_length_ratio)
        arrow_r2 = r2 * (1 - arrow_length_ratio)
        if direction == config.Direction.NONE:
            path_data = [
                (Path.MOVETO, (rad_start1, r1)),
                *make_arc_paths(rad_start1, rad_end1, r1),
                (Path.LINETO, (rad_end1, r1)),
                *make_bezier_paths(rad_end1, rad_end2, r1, r2, height_ratio),
                (Path.LINETO, (rad_end2, r2)),
                *make_arc_paths(rad_end2, rad_start2, r2),
                (Path.LINETO, (rad_start2, r2)),
                *make_bezier_paths(rad_start2, rad_start1, r2, r1, height_ratio),
                (Path.CLOSEPOLY, (rad_start1, r1)),
            ]
        elif direction == config.Direction.FORWARD:
            path_data = [
                (Path.MOVETO, (rad_start1, r1)),
                *make_arc_paths(rad_start1, rad_end1, r1),
                (Path.LINETO, (rad_end1, r1)),
                *make_bezier_paths(rad_end1, rad_end2, r1, arrow_r2, height_ratio),
                (Path.LINETO, (rad_end2, arrow_r2)),
                *make_arrow_paths(rad_end2, rad_start2, arrow_r2, r2),
                (Path.LINETO, (rad_start2, arrow_r2)),
                *make_bezier_paths(rad_start2, rad_start1, arrow_r2, r1, height_ratio),
                (Path.CLOSEPOLY, (rad_start1, r1)),
            ]
        elif direction == config.Direction.REVERSE:
            path_data = [
                (Path.MOVETO, (rad_start1, arrow_r1)),
                *make_arrow_paths(rad_start1, rad_end1, arrow_r1, r1),
                (Path.LINETO, (rad_end1, arrow_r1)),
                *make_bezier_paths(rad_end1, rad_end2, arrow_r1, r2, height_ratio),
                (Path.LINETO, (rad_end2, r2)),
                *make_arc_paths(rad_end2, rad_start2, r2),
                (Path.LINETO, (rad_start2, r2)),
                *make_bezier_paths(rad_start2, rad_start1, r2, arrow_r1, height_ratio),
                (Path.CLOSEPOLY, (rad_start1, arrow_r1)),
            ]
        elif direction == config.Direction.BIDIRECTIONAL:
            path_data = [
                (Path.MOVETO, (rad_start1, arrow_r1)),
                *make_arrow_paths(rad_start1, rad_end1, arrow_r1, r1),
                (Path.LINETO, (rad_end1, arrow_r1)),
                *make_bezier_paths(
                    rad_end1, rad_end2, arrow_r1, arrow_r2, height_ratio
                ),
                (Path.LINETO, (rad_end2, arrow_r2)),
                *make_arrow_paths(rad_end2, rad_start2, arrow_r2, r2),
                (Path.LINETO, (rad_start2, arrow_r2)),
                *make_bezier_paths(
                    rad_start2, rad_start1, arrow_r2, arrow_r1, height_ratio
                ),
                (Path.CLOSEPOLY, (rad_start1, arrow_r1)),
            ]
        else:
            err_msg = f"{direction=} is invalid value (0 or 1 or -1 or 2)."
            raise ValueError(err_msg)

        codes, verts = zip(*path_data)
        bezier_curve_path = Path(verts, codes, closed=True)
        super().__init__(bezier_curve_path, **kwargs)
