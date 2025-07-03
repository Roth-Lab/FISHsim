import math
import numpy as np
from scipy import interpolate


class SparseMatrix3D:
    """Representation of a 3-dimensional sparse matrix in COO format"""

    def __init__(
        self, data: np.ndarray, indexes: np.ndarray, shape: tuple, subpixel: bool = True
    ):
        """Creates an instance of sparse matrix

        Args:
            data (np.ndarray): Nx1 array representing the data stored in the matrix. N is the number of data points
            indexes (np.ndarray): Nx3 array representing the indexes (row, col, z) for the data points
            shape (tuple): (row, col, z) shape of the matrix. Could be floating point depending on the subpixel flag.
            subpixel (bool): indicates whether the matrix has subpixel indices
        """
        self.shape = shape
        self.data = data
        self.indexes = indexes  # Nx3 array of (row,col,z) indexes
        self.subpixel = subpixel


def sparse_convolve2d(input_matrix: SparseMatrix3D, kernel: np.ndarray) -> np.ndarray:
    """ Takes a sparse matrix and applies 2d convolution 
        Note: this function does not support subpixel 
    """

    N, M = input_matrix.shape
    n, m = kernel.shape
    output = np.zeros(shape=(N, M))

    indexes = input_matrix.indexes

    cnt = 0
    for _, val in np.ndenumerate(input_matrix.data):
        index = indexes[cnt, :]

        # Output matrix coordinates
        rl, ru = index[0] - math.floor(n / 2), index[0] + math.ceil(n / 2)
        rl, ru = rl if 0 <= rl else 0, ru if ru <= N else N
        cl, cu = index[1] - math.floor(m / 2), index[1] + math.ceil(m / 2)
        cl, cu = cl if 0 <= cl else 0, cu if cu <= M else M

        # Kernel matrix coordinates
        rl_k, ru_k = index[0] - math.floor(n / 2), index[0] + n - math.floor(n / 2)
        rl_k, ru_k = rl_k if 0 <= rl_k else 0, ru_k if ru_k <= N else N
        cl_k, cu_k = index[1] - math.floor(m / 2), index[1] + m - math.floor(m / 2)
        cl_k, cu_k = cl_k if 0 <= cl_k else 0, cu_k if cu_k <= M else M

        tl = np.array([index[0] - math.floor(n / 2), index[1] - math.floor(m / 2)])

        output[rl:ru, cl:cu] += (
            kernel[rl_k - tl[0] : ru_k - tl[0], cl_k - tl[1] : cu_k - tl[1]] * val
        )
        cnt += 1

    # Returning the output minus the padding
    return output[
        math.floor(n / 2) : -math.floor(n / 2), math.floor(m / 2) : -math.floor(m / 2)
    ]


def sparse_convolve3d(input_matrix: SparseMatrix3D, kernel: np.ndarray) -> np.ndarray:
    """Convolves a 3d sparse matrix with no padding (mode = valid).
       Currently only supports (matrix, kernel) pairs with the same height in z.
       Indexes need to be unique.

    Args:
        input_matrix (SparseMatrix3D): input sparse matrix(emitter position)
        kernel (np.ndarray): 3D convolution kernel (would be psf)

    Returns:
        np.ndarray: 2D matrix that represents the result of the convolution
    """
    M, N, K = input_matrix.shape
    mid_z_index = math.floor(K / 2)
    # Padding the output to ensure that convolution doesn't overflow.
    output = np.zeros((M + kernel.shape[1] - 1, N + kernel.shape[0] - 1))

    count = 0
    indexes = input_matrix.indexes + np.array(
        [math.floor(kernel.shape[1] / 2), math.floor(kernel.shape[0] / 2), 0]
    )  # row-wise addition. padding the indexes also.

    for val in input_matrix.data:
        # index is [row,col,z] for a single emitter
        index = indexes[count, :]

        if input_matrix.subpixel is True:
            # Z-interpolation (bi-lienar spline)
            section_bottom = (
                kernel[:, :, mid_z_index - (math.floor(index[2]) - mid_z_index)] * val
            )
            section_top = (
                kernel[:, :, mid_z_index - (math.ceil(index[2]) - mid_z_index)] * val
            )
            weight = index[2] % 1 #Get the floating point value
            section = (1 - weight) * section_bottom + weight * section_top #section that you want to interpolate

            # XY-interpolation (bi-cubic spline)
            rows = list(range(section.shape[0]))
            cols = list(range(section.shape[1]))
            #f is a function
            f = interpolate.RectBivariateSpline(rows, cols, section)

            x_shift = index[0] % 1
            y_shift = index[1] % 1
            rows_shifted = np.linspace(
                x_shift, section.shape[0] - 1 + x_shift, section.shape[0]
            )
            cols_shifted = np.linspace(
                y_shift, section.shape[1] - 1 + y_shift, section.shape[1]
            )
            section = f(rows_shifted, cols_shifted)
        else:
            section = kernel[:, :, mid_z_index - (index[2] - mid_z_index)] * val

        # Assigning the output matrix with the corresponding kernel section
        output[
            index[0]
            - math.floor(section.shape[1] / 2) : index[0]
            + math.ceil(section.shape[1] / 2),
            index[1]
            - math.floor(section.shape[0] / 2) : index[1]
            + math.ceil(section.shape[0] / 2),
        ] += section
        count += 1

    # Returning the output minus the padding
    return output[
        kernel.shape[1] - 1 : -kernel.shape[1] + 1,
        kernel.shape[0] - 1 : -kernel.shape[0] + 1,
    ]
