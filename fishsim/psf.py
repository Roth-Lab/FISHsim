import importlib.resources

import numpy as np
import scipy.io

import fishsim.resources


def load_psf():
    resource_dir = importlib.resources.files(fishsim.resources)

    psf_file = resource_dir.joinpath("psf.mat")

    psf = np.array(scipy.io.loadmat(psf_file)["ans"])

    return _process_psf(psf)


def _process_psf(psf):
    psf = psf.astype(np.float64)

    psf = psf / np.sum(psf)

    # Create binary mask for the psf
    radius = 14 / 2

    mask_x = np.arange(-np.ceil(psf.shape[0] / 2 - 1), np.ceil(psf.shape[0] / 2), 1)

    mask_y = np.arange(-np.ceil(psf.shape[1] / 2 - 1), np.ceil(psf.shape[1] / 2), 1)

    xv, yv = np.meshgrid(mask_x, mask_y)

    circle_mask = xv**2 + yv**2 >= radius**2

    # Apply the binary mask to the entire psf, substract average of the edge value
    for i in range(psf.shape[2]):
        zslice = psf[:, :, i]

        z_masked = zslice * circle_mask.astype(int)

        mean_edge = np.sum(z_masked) / np.count_nonzero(circle_mask)

        zslice = zslice - mean_edge

        zslice[zslice < 0] = 0

        psf[:, :, i] = zslice

    return psf
