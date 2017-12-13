import numpy as np

import pytest

from recalibration import solar_abundance
from recalibration import recalSingleLine


def test_solar_abundance():
    assert solar_abundance(26) == 7.47
    assert solar_abundance(42) == 1.88
    for element in range(1, 96):
        assert isinstance(solar_abundance(element), float)
    with pytest.raises(ValueError):
        solar_abundance('atom')


# def test_recalSingleLine():
#     line = np.array([4532.40, 26.0, 3.65, -1.800, 44.2])
#     params = (5777, 4.44, 0.00, 1.00)
#     loggf = recalSingleLine(line, params)
#     assert isinstance(loggf, float)
#     assert loggf != line[3]
