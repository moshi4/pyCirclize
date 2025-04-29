from __future__ import annotations

from copy import deepcopy
import math
import numpy as np
from typing import Any, Literal

from pycirclize_TEST import config

def degrees(rad: float) -> float:
    """Convert radian to positive degree (`0 - 360`)

    Parameters
    ----------
    rad : float
        Target radian

    Returns
    -------
    deg : float
        Positive degree (`0 - 360`)
    """
    # Radian to degree
    deg = math.degrees(rad)
    # Normalize degree in 0 - 360 range
    deg = deg % 360
    # Negative to positive
    if deg < 0:
        deg += 360
    return deg


def is_lower_loc(rad: float) -> bool:
    """Check target radian is lower location or not

    Parameters
    ----------
    rad : float
        Target radian

    Returns
    -------
    result : bool
        Lower location or not
    """
    deg = math.degrees(rad)
    return -270 <= deg < -90 or 90 <= deg < 270


def is_right_loc(rad: float) -> bool:
    """Check target radian is right location or not

    Parameters
    ----------
    rad : float
        Target radian

    Returns
    -------
    result : bool
        Right location or not
    """
    deg = math.degrees(rad)
    return -360 <= deg < -180 or 0 <= deg < 180


def is_ann_rad_shift_target_loc(rad: float) -> bool:
    """Check radian is annotation radian shift target or not

    Parameters
    ----------
    rad : float
        Annotation radian position

    Returns
    -------
    result : bool
        Target or not
    """
    deg = degrees(rad)
    return 30 <= deg <= 150 or 210 <= deg <= 330


def get_loc(
    rad: float,
) -> Literal["upper-right", "lower-right", "lower-left", "upper-left"]:
    """Get location of 4 sections

    Returns
    -------
    loc : str
        Location (`upper-right`|`lower-right`|`lower-left`|`upper-left`)
    """
    deg = degrees(rad)
    if 0 <= deg < 90:
        return "upper-right"
    elif 90 <= deg < 180:
        return "lower-right"
    elif 180 <= deg < 270:
        return "lower-left"
    else:
        return "upper-left"


def get_ann_relpos(rad: float) -> tuple[float, float]:
    """Get relative position for annotate by radian text position

    Parameters
    ----------
    rad : float
        Radian text position

    Returns
    -------
    relpos : tuple[float, float]
        Relative position
    """
    deg = degrees(rad)
    if 0 <= deg <= 180:
        return 0.0, Normalize(0, 180)(deg)
    else:
        return 1.0, 1.0 - Normalize(180, 360)(deg)


def plot_bbox(bbox: Bbox, ax: PolarAxes, **kwargs) -> None:
    """Plot bbox to check bbox area for development

    Parameters
    ----------
    bbox : Bbox
        Bounding box
    ax : PolarAxes
        Polar axes
    **kwargs : dict, optional
        Axes.plot properties (e.g. `color="red", lw=0.5, ls="--", ...`)
        <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html>
    """
    trans_bbox = bbox.transformed(ax.transAxes.inverted())
    kwargs.setdefault("clip_on", False)
    x0, y0, x1, y1 = trans_bbox.x0, trans_bbox.y0, trans_bbox.x1, trans_bbox.y1
    x, y = [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]
    ax.plot(x, y, transform=ax.transAxes, **kwargs)


def get_plotly_label_params(rad: float, adjust_rotation: bool, orientation: str, outer: bool = True, 
                            only_rotation: bool = False, **kwargs) -> dict:
    # Start with global defaults
    annotation = deepcopy(config.plotly_annotation_defaults)
    
    if not only_rotation:
        # Apply orientation-specific defaults
        orientation_defaults = config.plotly_text_orientation_defaults.get(orientation, {})
        annotation.update(orientation_defaults)
        
        # Handle outer/inner alignment
        if not outer:
            if orientation == "horizontal":
                annotation["yanchor"] = "bottom" if annotation["yanchor"] == "top" else "top"
            elif orientation == "vertical":
                annotation["xanchor"] = "right" if annotation["xanchor"] == "left" else "left"
        else:
            annotation["xanchor"] = "center"
            annotation["yanchor"] = "middle"
    
    # Override with user-provided kwargs
    annotation.update(kwargs)

    if adjust_rotation:
        rotation = np.degrees(rad)
        
        if orientation == "horizontal":
            rotation = rotation % 360
            # Flip if upside-down
            if 90 < rotation <= 270:
                rotation += 180
        elif orientation == "vertical":
            # Point text radially (90° offset from horizontal)
            rotation = (rotation + 90) % 360
            # No flipping needed for vertical text

        annotation.update({"textangle": rotation})

    return annotation

def build_plotly_shape(path: str, **kwargs) -> dict:
    shape_defaults = deepcopy(config.plotly_shape_defaults)
    shape_defaults.update(**kwargs)
    return {
        "type": "path",
        "path": path,
        **shape_defaults
    }
