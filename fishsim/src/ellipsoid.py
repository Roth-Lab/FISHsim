from __future__ import annotations  # for self referencing type hints
import math
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation


class Ellipsoid:
    """Class that represents an ellipsoid.
    Provide utility functions that allow checking for overlap between two ellipsoids

    Method adopted from:
    "Random generation of periodic hard ellipsoids based on molecular dynamics: A computationally-efficient algorithm" by Ghossein et al."
    """

    def __init__(self, center: list, axes: list, rot: list = None):
        """Constructs an ellipsoid object

        Args:
            axes (list): principal axes of the ellipsoid [a,b,c]
            center (list): center position of the ellipsoid
            rot (list, optional): euler angles [z, y, x] that represents rotation. If not provided, the angle is randomly generated.
        """

        self.axes = np.array(axes)
        self.center = np.array(center)
        if rot is None:
            rot = Ellipsoid.random_rotation()
        self.R = Rotation.from_euler("zyx", rot, degrees=True).as_matrix()

    def transformation_matrix(self) -> np.ndarray:
        """Creates a transformation matrix that represents rotation and tranaslation applied to an ellipsoid

        Returns:
            np.ndarray: Creates a [4 x 4] transformation matrix
        """

        M = np.eye(4)
        M[:3, :3] = self.R
        M[:3, 3] = self.center
        return M

    def to_matrix(self) -> np.ndarray:
        """Creates a matrix representation of an ellipsoid w/ rotation and translation reflected

        Returns:
            np.ndarray: [4 x 4] matrix representation of an ellipsoid
        """

        M = self.transformation_matrix()

        # Create an ellipsoid matrix w/ transformations and rotation applied
        # TODO: check whether the matrix invertible
        inv_M = np.linalg.inv(M)
        A = np.diag(np.append(1 / self.axes ** 2, -1))
        return inv_M.T @ A @ inv_M

    def random_point(
        self, radial_extent: list, num_points: int = 1, boundary_box_dim: list = None
    ) -> np.ndarray:
        """Generates a random point within the volume of the nucleus.
            One can specify which section of the nucleus to sample from by specifying radial extent


        Args:
            radial_extent (list): [inner, outer] specifies the section along the axes of the ellipsoid
                                  to sample from. Specify a float between 0 - 1
            num_points (int): number of points to generate
            boundary_box_dim (list): bounding box for the randomly generated points

        Returns:
            np.ndarray: [description]
        """
        prev_points = set()
        points = np.empty(shape=(0, 3))

        count = 0
        while count < num_points:
            # Compute the coordinates of the points in the form (theta, phi, radial)
            theta = np.random.uniform(0, 2 * np.pi, size=1)
            phi = np.arccos(1 - 2 * np.random.uniform(0, 1, size=1))
            radial = np.random.uniform(
                radial_extent[0] ** 3, radial_extent[1] ** 3, size=1
            ) ** (
                1 / 3
            )  # radial value is cube rooted to distribute the points according to the volume

            # Convert back to cartesian coordinate system
            x = self.axes[0] * np.sin(phi) * np.cos(theta) * radial
            y = self.axes[1] * np.sin(phi) * np.sin(theta) * radial
            z = self.axes[2] * np.cos(phi) * radial

            # Apply rotation
            coords = self.R @ np.vstack((x, y, z))

            # Apply translation
            coords = coords.T + self.center

            # Check if the cell is being cut off by the boundary condition
            if boundary_box_dim is not None and not (
                boundary_box_dim[0][0] <= coords[0, 0] < boundary_box_dim[0][1]
                and boundary_box_dim[1][0] <= coords[0, 1] < boundary_box_dim[1][1]
                and boundary_box_dim[2][0] <= coords[0, 2] < boundary_box_dim[2][1]
            ):
                continue

            # Checking for duplicates within the cell. Checking on (row, col) for now.
            coords_floored = np.floor(coords).astype(int)  # floor the coordinates
            if (coords_floored[0, 0], coords_floored[0, 1]) in prev_points:
                continue

            points = np.vstack((points, coords))
            prev_points.add((coords_floored[0, 0], coords_floored[0, 1]))
            count += 1

        return points

    def project(self) -> tuple:
        """ Projects the ellipsoid onto the xy-plane. 

            Returns:
                (major axis, minor axis, center) that define the ellipse in the plane.
                e.g. ([0,2],[1,0],[0,0])

            References:
            https://laurentlessard.com/teaching/cs524/slides/11%20-%20quadratic%20forms%20and%20ellipsoids.pdf
            https://astarmathsandphysics.com/university-maths-notes/matrices-and-linear-algebra/4507-lengths-and-directions-of-the-principal-axes-of-an-ellipse.html
            https://math.stackexchange.com/questions/573055/projection-of-ellipsoid?rq=1
            https://math.stackexchange.com/questions/874522/matrix-notation-of-an-ellipse
        """
        # Q: quadric matrix form of the ellipsoid (4x4)
        # P: projection matrix that defines the plane of projection (3x4)
        # C: conic matrix form of the projected ellipse (3x3)
        # C = (P*Q^-1*P.T)^-1
        Q = self.to_matrix()
        Q_inv = LA.inv(Q)
        P = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
        C = LA.inv(P @ Q_inv @ P.T)

        # The 2x2 matrix from the top left corner of C is the symmetric quadratic matrix of the ellipse
        # Note that an ellipse can be defined using the following quadratic: ax^2+2bxy+cy^2=1
        # From this quadratic matrix:
        #   eigen-vector: semi-axes directions
        #   eigen-values: 1/sqrt(lambda) = axes length
        w, v = LA.eig(C[:2, :2])  # eigenvalue, eigenvector
        axes_len = 1 / np.sqrt(w)
        ax1, ax2 = v[:, 0] * axes_len[0], v[:, 1] * axes_len[1]
        axes = sorted([ax1, ax2], key=LA.norm)

        return axes[1], axes[0], self.center[:2]  # major, minor, center


    def generate__surface(self) -> np.ndarray:
        points = np.empty(shape=(2, 3))

        theta = np.linspace(0, 2 * np.pi, num = 200)
        phi = np.arccos(1 - 2 * np.linspace(0, 1, num = 200))
        X,Y = np.meshgrid(theta, phi)
        x = self.axes[0] * np.sin(X) * np.cos(Y)
        # print(x.shape)
        y = self.axes[1] * np.sin(X) * np.sin(Y) 
        # print(y.shape)
        z = self.axes[2] * np.cos(X) 
        # print(z.shape)
        for i in range(200):
            for j in range(200):
            # Apply rotation
                coords = self.R @ np.vstack((x[i,j], y[i,j], z[i,j]))
            # Apply translation
                coords = coords.T + self.center

                points = np.vstack((points, coords))

        return points




    @staticmethod
    def random_rotation() -> list:
        """Creates a random rotation vector

        Returns:
            list: euler angles [x y z]
        """

        # euler's rotation
        angle_x = np.random.uniform(0, 180)
        angle_y = np.random.uniform(0, 180)
        angle_z = np.random.uniform(0, 180)

        return [angle_x, angle_y, angle_z]

       

    @staticmethod
    def __compute_p_coefficients(
        ellipsoid1: Ellipsoid, ellipsoid2: Ellipsoid
    ) -> np.ndarray:
        """
        Computes the coefficients of the characteristic equation of two ellipsoids defined as det(lambda*A+B).
        Refer to Ghossein et al.

        Returns:
            np.ndarray: [p1, p2, p3, p4, p5] coefficients of the 4th order characteristic equation
        """
        M1 = ellipsoid1.transformation_matrix()
        A2 = ellipsoid2.to_matrix()
        C = M1.T @ A2 @ M1

        sig = 1 / ellipsoid1.axes ** 2

        p1 = -np.prod(sig)
        p2 = -(
            sig[1] * sig[2] * C[0, 0]
            + sig[0] * sig[2] * C[1, 1]
            + sig[0] * sig[1] * C[2, 2]
            + p1 * C[3, 3]
        )
        p3 = (
            sig[0] * sig[1] * (C[2, 2] * C[3, 3] - C[2, 3] * C[3, 2])
            + sig[1] * sig[2] * (C[0, 0] * C[3, 3] - C[0, 3] * C[3, 0])
            + sig[0] * sig[2] * (C[1, 1] * C[3, 3] - C[1, 3] * C[3, 1])
            + sig[0] * (C[1, 2] * C[2, 1] - C[1, 1] * C[2, 2])
            + sig[1] * (C[0, 2] * C[2, 0] - C[0, 0] * C[2, 2])
            + sig[2] * (C[0, 1] * C[1, 0] - C[0, 0] * C[1, 1])
        )

        p4 = (
            sig[0]
            * (
                C[1, 1] * C[2, 2] * C[3, 3]
                - C[1, 1] * C[2, 3] * C[3, 2]
                - C[2, 2] * C[3, 1] * C[1, 3]
                - C[3, 3] * C[2, 1] * C[1, 2]
                + C[2, 1] * C[1, 3] * C[3, 2]
                + C[3, 1] * C[1, 2] * C[2, 3]
            )
            + sig[1]
            * (
                C[0, 0] * C[2, 2] * C[3, 3]
                - C[0, 0] * C[2, 3] * C[3, 2]
                - C[2, 2] * C[0, 3] * C[3, 0]
                - C[3, 3] * C[0, 2] * C[2, 0]
                + C[2, 0] * C[0, 3] * C[3, 2]
                + C[3, 0] * C[0, 2] * C[2, 3]
            )
            + sig[2]
            * (
                C[0, 0] * C[1, 1] * C[3, 3]
                - C[0, 0] * C[1, 3] * C[3, 1]
                - C[1, 1] * C[0, 3] * C[3, 0]
                - C[3, 3] * C[0, 1] * C[1, 0]
                + C[1, 0] * C[0, 3] * C[3, 1]
                + C[3, 0] * C[0, 1] * C[1, 3]
            )
            + C[0, 0] * C[1, 2] * C[2, 1]
            + C[1, 1] * C[0, 2] * C[2, 0]
            + C[2, 2] * C[0, 1] * C[1, 0]
            - C[0, 0] * C[1, 1] * C[2, 2]
            - C[1, 0] * C[0, 2] * C[2, 1]
            - C[2, 0] * C[0, 1] * C[1, 2]
        )
        p5 = np.linalg.det(A2)

        return np.array([p1, p2, p3, p4, p5])

    @staticmethod
    def __compute_n_coefficients(p: np.ndarray) -> np.ndarray:
        """
        Modifies the coefficients of the characteristic equation such that one can check for ellipsoid overlap
        Refer to Ghossein et al.

        Args:
            p (np.ndarray): coefficients of the characteristic equation

        Returns:
            np.ndarray: [n1, n2, n3, n4, n5] modified coefficients
        """
        p1 = -p[1] / (4 * p[0])
        p2 = p[2] / (6 * p[0])
        p3 = -p[3] / (4 * p[0])
        p4 = p[4] / p[0]

        b1 = (p4 - p3 * p1) + 3 * (p2 ** 2 - p1 * p3)
        b2 = -p3 * (p3 - p1 * p2) - p4 * (p1 ** 2 - p2) - p2 * (p2 ** 2 - p1 * p3)

        n1 = b1 ** 3 - 27 * b2 ** 2
        n2 = (
            -9 * (p3 - p1 * p2) ** 2
            + 27 * (p1 ** 2 - p2) * (p2 ** 2 - p1 * p3)
            - 3 * (p4 - p1 * p3) * (p1 ** 2 - p2)
        )
        n3 = b1 * (p3 - p1 * p2) - 3 * p1 * b2
        n4 = -(p4 - p1 * p3)
        n5 = p1 ** 2 - p2

        coeffs = [n1, n2, n3, n4, n5]
        # return np.array(coeffs)
        return np.array([n if abs(n) > 1e-8 else 0 for n in coeffs])

    @staticmethod
    def check_overlap(ellipsoid1: Ellipsoid, ellipsoid2: Ellipsoid) -> bool:
        """Checks for overlap between two ellipsoids. Externally touching is considered to non-overlapping.

        Returns:
            bool: true if overlapping, false otherwise
        """
        # Compute the distance between the centre of ellipsoid
        distance = LA.norm(ellipsoid1.center - ellipsoid2.center)

        # If the distance is greater than the sum of the maximum of the principle axes
        # of the ellipsoid then we know it's separated
        if distance >= (max(ellipsoid1.axes) + max(ellipsoid2.axes)):
            return False

        p = Ellipsoid.__compute_p_coefficients(ellipsoid1, ellipsoid2)
        n = Ellipsoid.__compute_n_coefficients(p)

        # Checking for complete seperation (externally tangent considered to be overlaping)
        return not (
            (n[0] == 0 and n[1] > 0 and n[2] > 0 and n[4] > 0)
            or (n[0] > 0 and n[1] > 0 and n[4] > 0)
            or (n[0] == 0 and n[1] > 0 and n[2] < 0 and n[4] > 0)
            or (n[0] == 0 and n[1] == 0 and n[3] < 0 and n[4] > 0)
        )

    @staticmethod
    def volume(axes: list) -> float:
        """Computes the volume of an ellipsoid

        Args:
            axes (list): principal axes of an ellipsoid

        Returns:
            float: volume of an ellipsoid
        """
        return (4 / 3) * np.pi * np.prod(np.array(axes))


if __name__ == "__main__":
    pass