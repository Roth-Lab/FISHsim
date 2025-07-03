import math
import numpy as np
import numpy.linalg as LA
from scipy.spatial.transform import Rotation as R
from .ellipsoid import Ellipsoid
from .utils import make_gaussian_2d


class Cell:
    def __init__(self):
        self.emitters = None


class EllipsoidCell(Cell):
    """
    A class to represent ellipsoid cell inherits from the Cell class, random ellipsoid representating cells
    """

    def __init__(
        self,
        center: list,
        axes_bounds: dict,
        rot: list = None,
        nucleus_size: float = 5 / 10,
        nucleus_emitter_density: float = 3 / 10,
    ):
        """Creates an EllipsoidCell object

        Args:
            center (list): center of the ellipsoid that defines cell shape
            axes_bounds (dict): min, max bounds on the principal axes of the ellisoid that defines cell shape
            nucleus_size (float): fractional length of the nucleus axes compared to the cell axes
            nucleus_emitter_density (float): fractional emitter density of the nucleus compared to the rest of the cell
        """
        super().__init__()
        if rot is None:
            self.shape = Ellipsoid(center, EllipsoidCell.random_axes(axes_bounds))
        else:
            self.shape = Ellipsoid(center, EllipsoidCell.random_axes(axes_bounds), rot)

        self.axes = self.shape.axes
        self.center = self.shape.center
        self.R = self.shape.R

        # Nucleus attributes
        self.nucleus_size = nucleus_size  # fractional length of the axes
        self.nucelus_emitter_density = (
            nucleus_emitter_density  # fractional emitter density
        )

    def generate_emitters(
        self, num_emitters: int, boundary_box_dim: list = None, is_nucleus: bool = False
    ) -> None:
        """Randomly generates the specified number of emitters within the volume of the cell

        Args:
            num_emitters (int): number of emitters to distribute
            boundary_box_dim (list, optional): external bounding box that limits the extent of the cell volume. Defaults to None.
            is_nucleus (bool, optional): if true, emitters are distributed in the nucleus with a lower emitter density. Defaults to False.
        """
        # Ellipsoid parametrization
        nucleus_emitter_count = 0

        if is_nucleus:
            # Compute how many emitters should be in the nucleus
            nucleus_volume = Ellipsoid.volume(self.nucleus_size * self.axes)
            nucleus_emitter_count = np.floor(
                (nucleus_volume / self.volume())
                * num_emitters
                * self.nucelus_emitter_density
            )
            nucleus_emitters = self.shape.random_point(
                [0, self.nucleus_size], nucleus_emitter_count, boundary_box_dim
            )

        # Sample the emitters in nucleus first
        # Append the sampled emitters into emitter_x,y,z
        self.emitters = self.shape.random_point(
            [self.nucleus_size, 1],
            num_emitters - nucleus_emitter_count,
            boundary_box_dim,
        )

        if is_nucleus:
            self.emitters = np.vstack((nucleus_emitters, self.emitters))

    def volume(self) -> float:
        """Computes volume of the cell

        Returns:
            volume (float): volume of the cell
        """
        return Ellipsoid.volume(self.axes)

    def background(self, bg_lvl) -> np.ndarray:
        """ Creates a 2D gaussian image based on the cell shape to serve as cellular background

            Args:
                bg_lvl (float): the gaussian image will be scaled such that max(img) == bg_lvl
        """

        major, minor, _ = self.shape.project()
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
            inner_extent=self.nucleus_size,
            outer_extent=1.2,
            smoothing=int(len_minor / 4),
        )
        gaussian_img *= bg_lvl / np.max(gaussian_img)

        return gaussian_img

    @staticmethod
    def random_axes(axes_bounds: dict) -> list:
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

    class ComplexCell(Cell):
        # TODO: Use a more biologically accurate model of a cell
        pass


if __name__ == "__main__":
  pass

