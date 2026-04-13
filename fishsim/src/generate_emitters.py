import math
import numpy as np
import time
from typing import Tuple
from .cells import EllipsoidCell
from .ellipsoid import Ellipsoid

def random_emitter_position(
    x_dim: tuple, y_dim: tuple, z_dim: tuple, num_emitter: int
) -> np.ndarray:
    """Generate randomly distributed emitter positions in a volume

    Args:
        x_dim (tuple): x dim of the volume
        y_dim (tuple): y dim of the volume
        z_dim (tuple): z dim of the volume
        num_emitter (int): number of emitter
    Returns:
        np.ndarray: [N x 3] array denoting (x, y, z) positions
    """

    # Uniformly distribute them in a 3D rectangular box
    pproc = np.random.uniform(0, 1, size=(num_emitter, 3))

    # rescale the position grid
    pproc[:, 0] = pproc[:, 0] * (x_dim[1] - x_dim[0]) + x_dim[0]
    pproc[:, 1] = pproc[:, 1] * (y_dim[1] - y_dim[0]) + y_dim[0]
    pproc[:, 2] = pproc[:, 2] * (z_dim[1] - z_dim[0]) + z_dim[0]

    return pproc


def cell_emitter_position(
    x_dim: tuple,
    y_dim: tuple,
    z_dim: tuple,
    num_emitter: int,
    cell_count: int,
    cell_axes_bounds: dict,
    is_nucleus: bool,
) -> Tuple[np.ndarray, list]:
    """Generates emitter positions taking cell structure into account

    Args:
        x_dim (tuple): x dim of the volume
        y_dim (tuple): y dim of the volume
        z_dim (tuple): z dim of the volume
        num_emitter (int): number of emitter
        cell_count (int): number of cell appear in the image
        cell_axes_bounds (dict): specifies the min, max bounds of the principle axes of the ellipsoid
    Returns:
        np.ndarray: [N x 3] array denoting (x, y, z) positions
    """

    cells = []  # list to store cell instances
    emitter_positions = np.empty(shape=(0, 3))
    total_volume = 0  # the total volume of all cells

    # TODO: take out the cell creation portion of this code to a seperate function.
    while len(cells) < cell_count:
        # Generate center position for cells
        cell_pos = np.random.uniform(0, 1, size=(3,))
        cell_pos = cell_pos * np.array(
            [x_dim[1] - x_dim[0], y_dim[1] - y_dim[0], z_dim[1] - z_dim[0]]
        ) + np.array([x_dim[0], y_dim[0], z_dim[0]])

        current_cell = EllipsoidCell(list(cell_pos), cell_axes_bounds)  # Create a cell
        overlap = False  # initial state for overlapping

        for cell in cells:
            if Ellipsoid.check_overlap(current_cell.shape, cell.shape):
                overlap = True  # if the new cells is overlapping with any of the newcell the flag will be set to True
                break

        if overlap != True:  # If the new cell does not overlap with any old cells
            cells.append(current_cell)  # append the new cell to the list of the cells
            total_volume += current_cell.volume()

    total_emitter = 0

    for i, cell in enumerate(cells):
        if i == len(cells) - 1:  # throw all the remaining emitter into the last cell
            emitter_count = num_emitter - total_emitter
        else:
            emitter_count = math.floor((num_emitter / total_volume) * cell.volume())
            total_emitter += emitter_count

        # NOTE: that xyz dim bounds aren't applied here. Could be added in the future
        cell.generate_emitters(
            emitter_count, np.array([x_dim, y_dim, z_dim]), is_nucleus
        )
        emitter_positions = np.vstack((emitter_positions, cell.emitters))

    return emitter_positions, cells



def enforce_min_center_distance(
    emitter_pos: np.ndarray,
    min_dist: float = 6.0,
    max_tries: int = 100000,
    x_dim: tuple = None,
    y_dim: tuple = None,
    z_dim: tuple = None,
) -> np.ndarray:
    """
    Enforce that all emitters are at least `min_dist` pixels apart
    in XY center-to-center distance.

    This uses rejection + resampling.

    Assumes emitter_pos has shape (N, 3) with columns [x, y, z].
    Spacing is checked in the XY plane only.
    """

    if len(emitter_pos) == 0:
        return emitter_pos

    if x_dim is None or y_dim is None or z_dim is None:
        raise ValueError("x_dim, y_dim, and z_dim must be provided for resampling.")

    accepted = []

    def is_valid(candidate, accepted_points):
        cx = int(candidate[0])
        cy = int(candidate[1])

        for pt in accepted_points:
            px = int(pt[0])
            py = int(pt[1])

            dist = ((cx - px) ** 2 + (cy - py) ** 2) ** 0.5
            if dist < min_dist:
                return False

        return True

    def sample_new_point():
        return np.array([
            np.random.uniform(x_dim[0], x_dim[1]),
            np.random.uniform(y_dim[0], y_dim[1]),
            np.random.uniform(z_dim[0], z_dim[1]),
        ])

    for pt in emitter_pos:
        if is_valid(pt, accepted):
            accepted.append(pt)
            continue

        placed = False
        for _ in range(max_tries):
            new_pt = sample_new_point()

            # If you're enforcing spacing after floor() in simulation.py,
            # keep the resampled point on the same integer grid too.
            new_pt[:2] = np.floor(new_pt[:2])

            if is_valid(new_pt, accepted):
                accepted.append(new_pt)
                placed = True
                break

        if not placed:
            raise RuntimeError(
                "Could not place all emitters with the requested minimum center distance. "
                "Reduce emitter density or spacing threshold."
            )

    return np.asarray(accepted)

if __name__ == "__main__":
    pass
