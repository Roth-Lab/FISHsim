import numpy as np
import numpy.linalg as LA
from ellipsoid import Ellipsoid


def test_ellipsoid_projection():
    Q = np.array([[1, -0.5, 0, 0], [-0.5, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
    inv_Q = LA.inv(Q)
    P = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    C = LA.inv(P @ inv_Q @ P.T)

    w, v = LA.eig(C[:2, :2])  # eigenvalue, eigenvector


if __name__ == "__main__":
    test_ellipsoid_projection()
