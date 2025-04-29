from __future__ import annotations

import math
from enum import IntEnum
from typing import ClassVar

import matplotlib as mpl

###########################################################
# Constant Value Config
###########################################################

# Fundamental Plot Parameters
MIN_R = 0
MAX_R = 100
R_PLOT_MARGIN = 5
ARC_RADIAN_STEP = 0.01
R_LIM = (MIN_R, MAX_R)
AXIS_FACE_PARAM = dict(zorder=0.99, ec="none", edgecolor="none")
AXIS_EDGE_PARAM = dict(zorder=1.01, fc="none", facecolor="none")
EPSILON = 1e-10

# Circos Color Scheme
# http://circos.ca/tutorials/lessons/configuration/colors/
CYTOBAND_COLORMAP = {
    "gpos100": "#000000",  # 0,0,0
    "gpos": "#000000",  # 0,0,0
    "gpos75": "#828282",  # 130,130,130
    "gpos66": "#A0A0A0",  # 160,160,160
    "gpos50": "#C8C8C8",  # 200,200,200
    "gpos33": "#D2D2D2",  # 210,210,210
    "gpos25": "#C8C8C8",  # 200,200,200
    "gvar": "#DCDCDC",  # 220,220,220
    "gneg": "#FFFFFF",  # 255,255,255
    "acen": "#D92F27",  # 217,47,39
    "stalk": "#647FA4",  # 100,127,164
}


class Direction(IntEnum):
    """Link BezierCurve Direction Enum"""

    REVERSE = -1
    NONE = 0
    FORWARD = 1
    BIDIRECTIONAL = 2


###########################################################
# Mutable Value Config (Mainly for Developer)
###########################################################


class _AnnotationAdjustConfig:
    """Annotation Position Adjustment Config"""

    enable: ClassVar[bool] = True
    """Enable Annotation position adjustment (default: `True`)"""
    limit: ClassVar[int] = 200
    """Limit of Annotation number for position adjustment (default: `200`)"""
    max_iter: ClassVar[int] = 1000
    """Max iteration number for Annotation position adjustment (default: `1000`)"""
    drad: ClassVar[float] = math.radians(0.1)
    """Delta radian for iterative position adjustment (default: `math.radians(0.1)`)"""
    dr: ClassVar[float] = 0.1
    """Delta radius for iterative position adjustment (default: `0.1`)"""
    expand: ClassVar[tuple[float, float]] = (1.2, 1.2)
    """Expand width & height factor of text bbox (default: `(1.2, 1.2)`)"""
    max_rad_shift: ClassVar[float] = math.radians(3.0)
    """Max radian of Annotation position shift (default: `math.radians(3.0)`)"""


clear_savefig: bool = True
"""
By default, after saving a figure using the `savefig()` method, figure object is
automatically deleted to avoid memory leaks (no display on jupyter notebook)
If you want to display the figure on jupyter notebook using `savefig()` method,
set clear_savefig=False.
"""
ann_adjust = _AnnotationAdjustConfig


###########################################################
# Matplotlib Runtime Configuration
###########################################################

# Setting matplotlib rc(runtime configuration) parameters
# https://matplotlib.org/stable/tutorials/introductory/customizing.html
mpl_rc_params = {
    # Font
    "font.size": 10,  # Default: 10
    # Lines
    "lines.linewidth": 0.5,  # Default: 1.5
    "lines.color": "black",  # Default: 'C0'
    "lines.markersize": 3,  # Default: 6
    # Patch
    "patch.linewidth": 0,  # Default: 1.0
    "patch.facecolor": "C0",  # Default: 'C0'
    "patch.edgecolor": "black",  # Default: 'black'
    # Legend
    "legend.loc": "upper left",  # Default: best
    "legend.frameon": False,  # Default: True
    "legend.handlelength": 1,  # Default: 2.0
    "legend.handleheight": 1,  # Default: 0.7
    # Savefig
    "savefig.bbox": "tight",  # Default: None
    "savefig.pad_inches": 0.5,  # Default: 0.1
    # SVG
    "svg.fonttype": "none",
}
mpl.rcParams.update(mpl_rc_params)

###########################################################
# GitHub Eukaryote & Prokaryote Dataset Config
###########################################################

GITHUB_DATA_URL = "https://raw.githubusercontent.com/moshi4/pycirclize-data/master/"

EUKARYOTE_DATASET = {
    "hg38": [
        "hg38_chr.bed",
        "hg38_cytoband.tsv",
        "hg38_genomic_link.tsv",
    ],
    "hs1": [
        "hs1_chr.bed",
        "hs1_cytoband.tsv",
        "hs1_genomic_link.tsv",
    ],
    "mm10": [
        "mm10_chr.bed",
        "mm10_cytoband.tsv",
        "mm10_genomic_link.tsv",
    ],
    "mm39": [
        "mm39_chr.bed",
        "mm39_cytoband.tsv",
        "mm39_genomic_link.tsv",
    ],
}

PROKARYOTE_FILES = [
    "enterobacteria_phage.gbk",
    "enterobacteria_phage.gff",
    "mycoplasma_alvi.gbk",
    "mycoplasma_alvi.gff",
    "escherichia_coli.gbk.gz",
    "escherichia_coli.gff.gz",
]
