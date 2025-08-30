from collections.abc import Sequence

from numpy.typing import NDArray

Numeric = int | float
NumericArrayLike = Sequence[Numeric] | NDArray
