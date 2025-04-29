from typing import ClassVar

import numpy as np
from pycirclizely_TEST import config


class PolarSVGPatchBuilder:
    step: ClassVar[float] = config.ARC_RADIAN_STEP

    @staticmethod
    def _polar_to_cart(theta: float, r: float) -> tuple[float, float]:
        """Convert polar coordinates to cartesian with Plotly's orientation:
        - Shift by π/2 to make 0 point upward
        - Negate to make angles increase clockwise
        """
        adjusted_theta = -(theta - np.pi / 2)
        x = r * np.cos(adjusted_theta)
        y = r * np.sin(adjusted_theta)
        return (x, y)

    @staticmethod
    def _svg_path_from_points(
        points: list[tuple[float, float]], closed: bool = False
    ) -> str:
        if not points:
            return ""
        path = [f"M {points[0][0]} {points[0][1]}"]
        path += [f"L {x} {y}" for x, y in points[1:]]
        if closed:
            path.append("Z")
        return " ".join(path)

    @classmethod
    def arc_rectangle(
        cls, radr: tuple[float, float], width: float, height: float
    ) -> str:
        min_rad, min_r = radr
        max_rad = min_rad + width
        max_r = min_r + height

        arc_rads = np.arange(min_rad, max_rad, cls.step)
        arc_rads = np.append(arc_rads, max_rad)

        bottom_arc = [cls._polar_to_cart(theta, min_r) for theta in arc_rads]
        upper_arc = [cls._polar_to_cart(theta, max_r) for theta in arc_rads[::-1]]
        points = bottom_arc + upper_arc + [bottom_arc[0]]
        return cls._svg_path_from_points(points, closed=True)

    @classmethod
    def arc_line(cls, rad_lim: tuple[float, float], r_lim: tuple[float, float]) -> str:
        rad_start, rad_end = rad_lim
        r_start, r_end = r_lim

        if rad_start > rad_end:
            rad_start, rad_end = rad_end, rad_start
            r_start, r_end = r_end, r_start

        if rad_start == rad_end:
            # Special case: straight radial line
            points = [
                cls._polar_to_cart(rad_start, r_start),
                cls._polar_to_cart(rad_end, r_end),
            ]
        else:
            arc_rads = np.arange(rad_start, rad_end, cls.step)
            arc_rads = np.append(arc_rads, rad_end)
            arc_rs = np.linspace(r_start, r_end, len(arc_rads))
            points = [cls._polar_to_cart(t, r) for t, r in zip(arc_rads, arc_rs)]

        return cls._svg_path_from_points(points)

    @classmethod
    def line(cls, rad_lim: tuple[float, float], r_lim: tuple[float, float]) -> str:
        points = [cls._polar_to_cart(rad, r) for rad, r in zip(rad_lim, r_lim)]
        return cls._svg_path_from_points(points)

    @classmethod
    def arc_arrow(
        cls,
        rad: float,
        r: float,
        drad: float,
        dr: float,
        head_length: float = np.pi / 90,
        shaft_ratio: float = 0.5,
    ) -> str:
        shaft_size = dr * shaft_ratio
        y_shaft_bottom = r + ((dr - shaft_size) / 2)
        y_shaft_upper = r + dr - ((dr - shaft_size) / 2)

        is_forward = drad >= 0
        drad = abs(drad)
        head_length = min(head_length, drad)

        rad_shaft_tip = (
            rad + (drad - head_length) if is_forward else rad - (drad - head_length)
        )
        rad_arrow_tip = rad + drad if is_forward else rad - drad

        p1 = rad, y_shaft_bottom
        p2 = rad_shaft_tip, y_shaft_bottom
        p3 = rad_shaft_tip, r
        p4 = rad_arrow_tip, (r + r + dr) / 2
        p5 = rad_shaft_tip, r + dr
        p6 = rad_shaft_tip, y_shaft_upper
        p7 = rad, y_shaft_upper

        shaft_rads = np.arange(p1[0], p2[0], cls.step if is_forward else -cls.step)
        bottom_shaft = [cls._polar_to_cart(t, p1[1]) for t in shaft_rads]
        top_shaft = [cls._polar_to_cart(t, p7[1]) for t in shaft_rads[::-1]]
        head = [cls._polar_to_cart(*pt) for pt in [p2, p3, p4, p5, p6]]
        path_points = bottom_shaft + head + top_shaft + [cls._polar_to_cart(*p1)]

        return cls._svg_path_from_points(path_points, closed=True)
