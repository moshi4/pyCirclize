from __future__ import annotations

import math

import pytest

from pycirclizely.sector import Sector
from pycirclizely.track import Track


@pytest.fixture
def track() -> Track:
    """Track fixture"""
    sector = Sector(
        name="A",
        size=1000,
        rad_lim=(0, math.pi),
    )
    track = Track(
        name="Track01",
        r_lim=(90, 100),
        r_pad_ratio=0.1,
        parent_sector=sector,
    )
    return track


def test_track_property(track: Track):
    """Test track property"""
    assert track.name == "Track01"
    assert track.size == 1000
    assert track.start == 0
    assert track.end == 1000
    assert track.center == 500
    assert track.r_size == 10
    assert track.r_lim == (90, 100)
    assert track.r_center == 95
    assert track.r_plot_size == 9
    assert track.r_plot_lim == (90.5, 99.5)
    assert track.rad_size == math.pi
    assert track.rad_lim == (0, math.pi)
    assert track.deg_size == 180
    assert track.deg_lim == (0, 180)
