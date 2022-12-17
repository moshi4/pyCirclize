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

## Not Implemented Features

List of features implemented in other Circos plotting tools but not yet implemented in pyCirclize.
I may implement them when I feel like it.

- Plot histogram
- Plot boxplot
- Plot violin
- Plot raster image
- Label position auto adjustment
