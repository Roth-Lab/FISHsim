import math

import numpy as np
import numpy.linalg as LA

from fishsim.ellipsoid import Ellipsoid
from fishsim.utils import make_gaussian_2d


class Cell:
    def __init__(self):
        self.emitters = None


class EllipsoidCell(Cell):
    """
    A class to represent ellipsoid cell inherits from the Cell class, random ellipsoid representating cells
    """

    def __init__(
        self,
        axes_bounds,
        center,
        nucleus_size=0.5,
        nucleus_emitter_density=0.3,
        rot=None,
    ):
        """Creates an EllipsoidCell object

        Args:
            center (list): center of the ellipsoid that defines cell shape
            axes_bounds (dict): min, max bounds on the principal axes of the ellisoid that defines cell shape
            nuc_size (float): fractional length of the nucleus axes compared to the cell axes
            nucleus_emitter_density (float): fractional emitter density of the nucleus compared to the rest of the cell
            rot (list, optional): euler angles [z, y, x] that represents rotation. If not provided, the angle is randomly generated.
        """
        super().__init__()

        if rot is None:
            self.shape = Ellipsoid(Ellipsoid.random_axes(axes_bounds), center)

        else:
            self.shape = Ellipsoid(Ellipsoid.random_axes(axes_bounds), center, rot=rot)

        self.nuc_size = nucleus_size

        self.nucelus_emitter_density = nucleus_emitter_density

    @property
    def axes(self):
        return self.shape.axes

    @property
    def center(self):
        return self.shape.center

    @property
    def R(self):
        return self.shape.R

    @property
    def volume(self):
        return self.shape.volume

    def background(self, background_level):
        """Creates a 2D gaussian image based on the cell shape to serve as cellular background

        Args:
            background_level (float): the gaussian image will be scaled such that max(img) == background_level
        """
        major, minor, _ = self.shape.projection

        len_major, len_minor = LA.norm(major), LA.norm(minor)

        std = [len_major, len_minor]

        shape = (4 * math.ceil(len_major), 4 * math.ceil(len_major))

        # Compute the rotation angle. Ensure that the major axis is in the righ quadrant
        major = -major if major[1] < 0 else major

        rot = math.acos(np.dot(major, np.array([1, 0])) / len_major)  # radians

        # Create the base gaussian image
        gaussian_img = make_gaussian_2d(
            std,
            shape,
            rot,
            inner_extent=self.nuc_size,
            outer_extent=1.2,
            smoothing=int(len_minor / 4),
        )

        gaussian_img *= background_level / np.max(gaussian_img)

        return gaussian_img

    def check_overlap(self, x):
        """Check for overlap between this cell and another cell.

        Args:
            x (EllipsoidCell): Cell to check for overlap with

        Returns:
            bool: Whether the cell overlaps with this one.
        """
        return self.shape.check_overlap(x.shape)

    def generate_emitters(self, boundary_box, num_emitters, sim_nucleus=True):
        """Randomly generates the specified number of emitters within the volume of the cell

        Args:
            boundary_box (dict): external bounding box that limits the extent of the cell volume.
            num_emitters (int): number of emitters to distribute
            is_nucleus (bool, optional): if true, emitters are distributed in the nucleus with a lower emitter density. Defaults to False.
        """
        emitters = []

        if sim_nucleus:
            nuc_vol = Ellipsoid(self.nuc_size * self.axes, self.center).volume

            num_nuc_emitters = np.floor((nuc_vol / self.volume) * num_emitters * self.nucelus_emitter_density)

            emitters.extend(
                self.shape.generate_random_points(boundary_box, [0, self.nuc_size], num_points=num_nuc_emitters)
            )

        else:
            num_nuc_emitters = 0

        # Append the sampled emitters into emitter_x,y,z
        num_cyto_emitters = num_emitters - num_nuc_emitters

        emitters.extend(
            self.shape.generate_random_points(boundary_box, [self.nuc_size, 1], num_points=num_cyto_emitters)
        )

        return np.array(emitters)
