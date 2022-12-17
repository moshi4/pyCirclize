from __future__ import annotations

import warnings
from enum import IntEnum

import matplotlib as mpl

warnings.filterwarnings("ignore")


###########################################################
# Constant Value Config
###########################################################

# Fundamental Plot Parameters
MIN_R = 0
MAX_R = 100
R_PLOT_MARGIN = 5
ARC_RADIAN_STEP = 0.001
R_LIM = (MIN_R, MAX_R)
AXIS_FACE_PARAM = dict(zorder=0.99, ec="none", edgecolor="none")
AXIS_EDGE_PARAM = dict(zorder=1.01, fc="none", facecolor="none")

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
    "mm10": [
        "mm10_chr.bed",
        "mm10_cytoband.tsv",
        "mm10_genomic_link.tsv",
    ],
}

PROKARYOTE_FILES = [
    "enterobacteria_phage.gbk",
    "enterobacteria_phage.gff",
    "escherichia_coli.gbk.gz",
    "escherichia_coli.gff.gz",
]
