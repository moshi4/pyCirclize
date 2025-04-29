from __future__ import annotations

import pytest

from pycirclizely import Circos


def test_circos_init():
    """Test circos initialization"""
    circos = Circos({"A": 10, "B": 20, "C": 15})
    assert [s.name for s in circos.sectors] == ["A", "B", "C"]
    assert [s.size for s in circos.sectors] == [10, 20, 15]

    circos = Circos(dict(D=10, E=(10, 20), F=(30, 50), G=100))
    assert [s.name for s in circos.sectors] == ["D", "E", "F", "G"]
    assert [s.size for s in circos.sectors] == [10, 10, 20, 100]


@pytest.mark.parametrize(
    "start, end",
    [
        (-10, 360),  # End - Start > 360
        (0, -90),  # Start > End
        (-400, -200),  # Start < -360
        (200, 400),  # End > 360
    ],
)
def test_circos_init_range_error(start: float, end: float):
    """Test circos initialization range error"""
    with pytest.raises(ValueError):
        Circos({s: 10 for s in "ABC"}, start=start, end=end)


@pytest.mark.parametrize(
    "space, endspace, success",
    [
        # If `sector num = 3` and `endspace = True`
        # List of space length must be 3
        ([5], True, False),
        ([5, 10], True, False),
        ([5, 10, 15], True, True),
        ([5, 10, 15, 20], True, False),
        # If `sector num = 3` and `endspace = False`
        # List of space length must be 2
        ([5], False, False),
        ([5, 10], False, True),
        ([5, 10, 15], False, False),
    ],
)
def test_circos_init_space_list(space: list[float], endspace: bool, success: bool):
    """Test circos initialization space list length error"""
    sectors = {s: 10 for s in "ABC"}
    if success:
        Circos(sectors, space=space, endspace=endspace)
    else:
        with pytest.raises(ValueError):
            Circos({s: 10 for s in "ABC"}, space=space, endspace=endspace)


def test_get_sector():
    """Test `get_sector()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors)
    # Case1: Successfully get sector
    for sector_name in sectors.keys():
        circos.get_sector(sector_name)
    # Case2: Failed to get sector
    with pytest.raises(ValueError):
        circos.get_sector("error")


def test_get_group_sectors_deg_lim():
    """Test `get_group_sectors_deg_lim()`"""
    sectors = dict(A=10, B=10, C=10, D=10, E=10, F=10, G=10, H=10)

    group1 = list("BCD")
    circos = Circos(sectors)
    group1_deg_lim = circos.get_group_sectors_deg_lim(group1)
    assert tuple(map(round, group1_deg_lim)) == (45, 180)

    group2 = list("HEF")
    circos = Circos(sectors, start=20, end=340)
    group2_deg_lim = circos.get_group_sectors_deg_lim(group2)
    assert tuple(map(round, group2_deg_lim)) == (180, 340)


def test_ax_property():
    """Test `circos.ax` property"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors)
    # Raise error before calling `circos.plotfig()` method
    with pytest.raises(ValueError):
        assert circos.ax
    circos.plotfig()
    assert circos.ax
