from scipy.spatial.transform import Rotation

import numpy as np
import numpy.linalg as LA


class Ellipsoid(object):
    """Class that represents an ellipsoid.
    Provide utility functions that allow checking for overlap between two ellipsoids

    Method adopted from:
    "Random generation of periodic hard ellipsoids based on molecular dynamics: A computationally-efficient algorithm" by Ghossein et al."
    """

    @staticmethod
    def random_axes(axes_bounds):
        """Creates random principal axes for the ellipsoid based on the bounds specified
            in axes_bounds

        Args:
            axes_bounds (dict): specifies the min, max bounds for each of the axes

        Returns:
            list: randomly generated pricipal axes of an ellipsoid [a,b,c]
        """
        a = np.random.uniform(axes_bounds["a"][0], axes_bounds["a"][1])

        b = np.random.uniform(axes_bounds["b"][0], axes_bounds["b"][1])

        c = np.random.uniform(axes_bounds["c"][0], axes_bounds["c"][1])

        return [a, b, c]

    def __init__(self, axes, center, rot=None):
        """Constructs an ellipsoid object

        Args:
            axes (list): principal axes of the ellipsoid [a,b,c]
            center (list): center position of the ellipsoid
            rot (list, optional): euler angles [z, y, x] that represents rotation. If not provided, the angle is randomly generated.
        """
        self.axes = np.array(axes)

        self.center = np.array(center)

        if rot is None:
            rot = [np.random.uniform(0, 180), np.random.uniform(0, 180), np.random.uniform(0, 180)]

        self.R = Rotation.from_euler("zyx", rot, degrees=True).as_matrix()

    @property
    def projection(self):
        """Projects the ellipsoid onto the xy-plane.

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

    @property
    def volume(self):
        """Computes the volume of an ellipsoid

        Args:
            axes (list): principal axes of an ellipsoid

        Returns:
            float: volume of an ellipsoid
        """
        return (4 / 3) * np.pi * np.prod(np.array(self.axes))

    def check_overlap(self, x):
        """Checks for overlap between this ellipsoid and another.
        Externally touching is considered to non-overlapping.

        Args:
            x (Ellipsoid): Shape to check if overlaps with this one.

        Returns:
            bool: true if overlapping, false otherwise
        """
        # Compute the distance between the centre of ellipsoid
        distance = LA.norm(self.center - x.center)

        # If the distance is greater than the sum of the maximum of the principle axes
        # of the ellipsoid then we know it's separated
        if distance >= (max(self.axes) + max(x.axes)):
            return False

        n = self._compute_n_coefficients(x)

        # Checking for complete seperation (externally tangent considered to be overlaping)
        return not (
            (n[0] == 0 and n[1] > 0 and n[2] > 0 and n[4] > 0)
            or (n[0] > 0 and n[1] > 0 and n[4] > 0)
            or (n[0] == 0 and n[1] > 0 and n[2] < 0 and n[4] > 0)
            or (n[0] == 0 and n[1] == 0 and n[3] < 0 and n[4] > 0)
        )

    def generate_random_points(self, boundary_box, radial_extent, min_dist=0, num_points=1):
        """Generates a random point within the volume of the nucleus.
            One can specify which section of the nucleus to sample from by specifying radial extent


        Args:
            boundary_box (dict): bounding box for the randomly generated points
            radial_extent (list): [inner, outer] specifies the section along the axes of the ellipsoid
                                  to sample from. Specify a float between 0 - 1
            min_dist (float): minimum distance in x-y plane points must be apart
            num_points (int): number of points to generate

        Returns:
            np.ndarray: [description]
        """

        def is_in_boundary_box(coords):
            return (
                boundary_box["y"][0] <= coords[0] < boundary_box["y"][1]
                and boundary_box["x"][0] <= coords[1] < boundary_box["x"][1]
                and boundary_box["z"][0] <= coords[2] < boundary_box["z"][1]
            )

        def is_valid(candidate, accepted_points):
            if min_dist <= 0:
                return True
            cx = int(candidate[0])
            cy = int(candidate[1])
            for pt in accepted_points:
                px = int(pt[0])
                py = int(pt[1])
                dist = ((cx - px) ** 2 + (cy - py) ** 2) ** 0.5
                if dist < min_dist:
                    return False
            return True

        points = []

        while len(points) < num_points:
            # Compute the coordinates of the points in the form (theta, phi, radial)
            theta = np.random.uniform(0, 2 * np.pi, size=1)

            phi = np.arccos(1 - 2 * np.random.uniform(0, 1, size=1))

            radial = np.random.uniform(radial_extent[0] ** 3, radial_extent[1] ** 3, size=1) ** (
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

            coords = np.squeeze(coords)

            if is_in_boundary_box(coords) and is_valid(coords, points):
                points.append(coords)

        return np.array(points)

    def to_matrix(self):
        """Creates a matrix representation of an ellipsoid w/ rotation and translation reflected

        Returns:
            np.ndarray: [4 x 4] matrix representation of an ellipsoid
        """
        M = self.transformation_matrix()

        # Create an ellipsoid matrix w/ transformations and rotation applied
        # TODO: check whether the matrix invertible
        inv_M = np.linalg.inv(M)

        A = np.diag(np.append(1 / self.axes**2, -1))

        return inv_M.T @ A @ inv_M

    def transformation_matrix(self):
        """Creates a transformation matrix that represents rotation and tranaslation applied to an ellipsoid

        Returns:
            np.ndarray: Creates a [4 x 4] transformation matrix
        """
        M = np.eye(4)

        M[:3, :3] = self.R

        M[:3, 3] = self.center

        return M

    def _compute_p_coefficients(self, x):
        """
        Computes the coefficients of the characteristic equation of two ellipsoids defined as det(lambda*A+B).
        Refer to Ghossein et al.

        Returns:
            np.ndarray: [p1, p2, p3, p4, p5] coefficients of the 4th order characteristic equation
        """
        M1 = self.transformation_matrix()

        A2 = x.to_matrix()

        C = M1.T @ A2 @ M1

        sig = 1 / self.axes**2

        p1 = -np.prod(sig)

        p2 = -(sig[1] * sig[2] * C[0, 0] + sig[0] * sig[2] * C[1, 1] + sig[0] * sig[1] * C[2, 2] + p1 * C[3, 3])

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

    def _compute_n_coefficients(self, x):
        """
        Modifies the coefficients of the characteristic equation such that one can check for ellipsoid overlap
        Refer to Ghossein et al.

        Args:
            p (np.ndarray): coefficients of the characteristic equation

        Returns:
            np.ndarray: [n1, n2, n3, n4, n5] modified coefficients
        """
        p = self._compute_p_coefficients(x)

        p1 = -p[1] / (4 * p[0])
        p2 = p[2] / (6 * p[0])
        p3 = -p[3] / (4 * p[0])
        p4 = p[4] / p[0]

        b1 = (p4 - p3 * p1) + 3 * (p2**2 - p1 * p3)
        b2 = -p3 * (p3 - p1 * p2) - p4 * (p1**2 - p2) - p2 * (p2**2 - p1 * p3)

        n1 = b1**3 - 27 * b2**2
        n2 = -9 * (p3 - p1 * p2) ** 2 + 27 * (p1**2 - p2) * (p2**2 - p1 * p3) - 3 * (p4 - p1 * p3) * (p1**2 - p2)
        n3 = b1 * (p3 - p1 * p2) - 3 * p1 * b2
        n4 = -(p4 - p1 * p3)
        n5 = p1**2 - p2

        coeffs = [n1, n2, n3, n4, n5]
        # return np.array(coeffs)
        return np.array([n if abs(n) > 1e-8 else 0 for n in coeffs])


if __name__ == "__main__":
    pass
