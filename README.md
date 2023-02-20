# pyCirclize: Circular visualization in Python

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![OS](https://img.shields.io/badge/OS-_Windows_|_Mac_|_Linux-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)
[![Latest PyPI version](https://img.shields.io/pypi/v/pycirclize.svg)](https://pypi.python.org/pypi/pycirclize)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/pycirclize.svg?color=green)](https://anaconda.org/conda-forge/pycirclize)
[![CI](https://github.com/moshi4/pyCirclize/actions/workflows/ci.yml/badge.svg)](https://github.com/moshi4/pyCirclize/actions/workflows/ci.yml)

## Table of contents

- [Overview](#overview)
- [Installation](#installation)
- [API Usage](#api-usage)
- [Code Example](#code-example)
- [Not Implemented Features](#not-implemented-features)

## Overview

pyCirclize is a circular visualization python package implemented based on matplotlib.
This package is developed for the purpose of easily and beautifully plotting circular figure such as Circos Plot and Chord Diagram in Python.
In addition, useful genome and phylogenetic tree visualization methods for the bioinformatics field are also implemented.
pyCirclize was inspired by [circlize](https://github.com/jokergoo/circlize) and [pyCircos](https://github.com/ponnhide/pyCircos).
More detailed documentation is available [here](https://moshi4.github.io/pyCirclize/).

![pyCirclize_gallery.png](https://raw.githubusercontent.com/moshi4/pyCirclize/main/docs/images/pyCirclize_gallery.png)  
**Fig.1 pyCirclize example plot gallery**

## Installation

`Python 3.8 or later` is required for installation.

**Install PyPI package:**

    pip install pycirclize

**Install conda-forge package:**

    conda install -c conda-forge pycirclize

## API Usage

API usage is described in each of the following sections in the [document](https://moshi4.github.io/pyCirclize/).

- [Getting Started](https://moshi4.github.io/pyCirclize/getting_started/)
- [Plot API Example](https://moshi4.github.io/pyCirclize/plot_api_example/)
- [Chord Diagram](https://moshi4.github.io/pyCirclize/chord_diagram/)
- [Circos Plot (Genomics)](https://moshi4.github.io/pyCirclize/circos_plot/)
- [Phylogenetic Tree](https://moshi4.github.io/pyCirclize/phylogenetic_tree/)

## Code Example

### 1. Circos Plot

```python
from pycirclize import Circos
import numpy as np
np.random.seed(0)

# Initialize Circos sectors
sectors = {"A": 10, "B": 15, "C": 12, "D": 20, "E": 15}
circos = Circos(sectors, space=5)

for sector in circos.sectors:
    # Plot sector name
    sector.text(f"Sector: {sector.name}", r=110, size=15)
    # Create x positions & randomized y values
    x = np.arange(sector.start, sector.end) + 0.5
    y = np.random.randint(0, 100, len(x))
    # Plot line track
    line_track = sector.add_track((80, 100), r_pad_ratio=0.1)
    line_track.xticks_by_interval(interval=1)
    line_track.axis()
    line_track.line(x, y)
    # Plot points track
    points_track = sector.add_track((55, 75), r_pad_ratio=0.1)
    points_track.axis()
    points_track.scatter(x, y)
    # Plot bar track
    bar_track = sector.add_track((30, 50), r_pad_ratio=0.1)
    bar_track.axis()
    bar_track.bar(x, y)

# Plot links 
circos.link(("A", 0, 3), ("B", 15, 12))
circos.link(("B", 0, 3), ("C", 7, 11), color="skyblue")
circos.link(("C", 2, 5), ("E", 15, 12), color="chocolate", direction=1)
circos.link(("D", 3, 5), ("D", 18, 15), color="lime", ec="black", lw=0.5, hatch="//", direction=2)
circos.link(("D", 8, 10), ("E", 2, 8), color="violet", ec="red", lw=1.0, ls="dashed")

circos.savefig("example01.png")
```

![example01.png](https://raw.githubusercontent.com/moshi4/pyCirclize/main/docs/images/example01.png)  

### 2. Chord Diagram

```python
from pycirclize import Circos
import pandas as pd

# Create matrix dataframe (3 x 6)
row_names = ["F1", "F2", "F3"]
col_names = ["T1", "T2", "T3", "T4", "T5", "T6"]
matrix_data = [
    [10, 16, 7, 7, 10, 8],
    [4, 9, 10, 12, 12, 7],
    [17, 13, 7, 4, 20, 4],
]
matrix_df = pd.DataFrame(matrix_data, index=row_names, columns=col_names)

# Initialize Circos from matrix for plotting Chord Diagram
circos = Circos.initialize_from_matrix(
    matrix_df,
    space=5,
    cmap="tab10",
    label_kws=dict(size=12),
    link_kws=dict(ec="black", lw=0.5, direction=1),
)

circos.savefig("example02.png")
```

![example02.png](https://raw.githubusercontent.com/moshi4/pyCirclize/main/docs/images/example02.png)  

## Not Implemented Features

List of features implemented in other Circos plotting tools but not yet implemented in pyCirclize.
I may implement them when I feel like it.

- Plot histogram
- Plot boxplot
- Plot violin
- Plot raster image
- Label position auto adjustment
