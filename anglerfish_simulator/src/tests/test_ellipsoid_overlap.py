import math
import numpy as np
import sys
sys.path.insert(1, '/Users/luna/Documents/GitHub/anglerfish_simulator/src/ellipsoid.py')
from ellipsoid import Ellipsoid


def test_matrix_representation():
    e1 = Ellipsoid([0, 0, 0], [3, 3, 3])
    x = np.array([0, 0, 3, 1]).reshape((4, 1))
    return math.isclose(x.T @ e1.to_matrix() @ x, 0, abs_tol=1e-10)


def test_overlapping():
    e1 = Ellipsoid([0, 0, 0], [1, 1, 1])
    e2 = Ellipsoid([0, 0, 3], [3, 3, 3])
    assert Ellipsoid.check_overlap(e1, e2) == True


def test_not_overlapping():
    e1 = Ellipsoid([0, 0, 0], [1, 1, 1])
    e2 = Ellipsoid([0, 0, 3], [1, 1, 1])
    assert Ellipsoid.check_overlap(e1, e2) == False


def test_completely_overlap():
    e1 = Ellipsoid([0, 0, 0.1], [1, 1, 1])
    e2 = Ellipsoid([0, 0, 0], [1, 1, 1])
    assert Ellipsoid.check_overlap(e1, e2) == True


def test_touch_but_not_overlap():
    e1 = Ellipsoid([1, 0, 3], [2, 3, 6], [0, 0, 0])
    e2 = Ellipsoid([1, 5, 3], [2, 2, 2], [0, 0, 0])
    assert Ellipsoid.check_overlap(e1, e2) == False


def test_separate_but_close_to_overlap():
    e1 = Ellipsoid([1, 0, 3], [2, 3, 6], [0, 0, 0])
    e2 = Ellipsoid([0, 5, 3], [2, 2, 2], [0, 0, 0])
    assert Ellipsoid.check_overlap(e1, e2) == False


def test_rotated_90_degrees_touch_but_not_overlap():
    e1 = Ellipsoid([4, 0, 0], [2, 2, 2], [0, 0, 0])
    e2 = Ellipsoid([0, 0, 0], [2, 6, 3], [90, 0, 0])
    assert Ellipsoid.check_overlap(e1, e2) == False


def test_rotated_90_degrees_overlap():
    e1 = Ellipsoid([2, 0, 1], [2, 2, 2], [0, 0, 0])
    e2 = Ellipsoid([0, 0, 0], [2, 6, 3], [90, 0, 0])
    assert Ellipsoid.check_overlap(e1, e2) == True


if __name__ == "__main__":
    heyhey = test_matrix_representation()
    print(heyhey)