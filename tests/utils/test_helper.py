import pytest

from pycirclize.utils import ColorCycler, calc_group_spaces


def test_color_cycler():
    """Test color cycler"""
    # Check get color list length
    ColorCycler.set_cmap("tab10")
    assert len(ColorCycler.get_color_list()) == 10
    assert len(ColorCycler.get_color_list(5)) == 5
    assert len(ColorCycler.get_color_list(20)) == 20

    # Check cycle index, color
    assert ColorCycler(0) != ColorCycler(1)
    assert ColorCycler(0) == ColorCycler(10)
    assert ColorCycler(15) == ColorCycler(25)

    # Check cycle counter
    assert ColorCycler() != ColorCycler()
    assert ColorCycler.counter == 2

    # Check reset cycle
    ColorCycler.reset_cycle()
    assert ColorCycler.counter == 0

    # Check cmap change
    ColorCycler.set_cmap("tab20")
    with pytest.raises(KeyError):
        ColorCycler.set_cmap("invalid name")
    assert len(ColorCycler.get_color_list()) == 20


def test_calc_group_spaces():
    """Test `calc_group_spaces`"""
    # Case1. Blank list (error)
    with pytest.raises(ValueError):
        calc_group_spaces([])

    # Case2. List length = 1 (endspace=True)
    spaces = calc_group_spaces([5])
    expected_spaces = [2, 2, 2, 2, 2]
    assert spaces == expected_spaces

    # Case3. List length = 1 (endspace=False)
    spaces = calc_group_spaces([5], space_in_group=3, endspace=False)
    expected_spaces = [3, 3, 3, 3]
    assert spaces == expected_spaces

    # Case4. List length > 1 (endspace=True)
    spaces = calc_group_spaces([4, 3, 3])
    expected_spaces = [2, 2, 2, 5, 2, 2, 5, 2, 2, 5]
    assert spaces == expected_spaces

    # Case5. List length > 1 (endspace=False)
    spaces = calc_group_spaces(
        [4, 3, 3], space_bw_group=8, space_in_group=1, endspace=False
    )
    expected_spaces = [1, 1, 1, 8, 1, 1, 8, 1, 1]
    assert spaces == expected_spaces
