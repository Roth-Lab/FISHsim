import numpy as np

from fishsim.sparse import SparseMatrix3D, sparse_convolve2d, sparse_convolve3d
from fishsim.utils import TruncatedNormal, make_gaussian_2d


class ImageSimulator(object):
    """Simulate FISH image"""

    def __init__(
        self,
        data_org,
        sim_background,
        sim_camera,
        sim_photon,
    ):
        self.data_org = data_org

        self.sim_background = sim_background

        self.sim_camera = sim_camera

        self.sim_photon = sim_photon

    def get_bit_img(self, bit, probed_fov):
        img_round = self.data_org.get_round(bit)

        img = self.sim_background.get_img(probed_fov.fov, img_round)

        bit_num = self.data_org.get_bit_number(bit)

        img += self.sim_photon.get_img(bit_num, img_round, probed_fov)

        channel = self.data_org.get_channel(bit)

        return self.sim_camera.capture_img(channel, img)


class BackgroundSimulator(object):
    """Simulate image background"""

    def __init__(self, bg_sampling_prob, photon_counts, rng, signal_to_cell_ratio):
        self.bg_sampling_prob = bg_sampling_prob

        self.photon_counts = photon_counts

        self.rng = rng

        self.signal_to_cell_ratio = signal_to_cell_ratio

        self.cell_bg_to_global_bg_ratio = [self.rng.uniform(1.15, 2.5) for _ in photon_counts]

    def get_img(self, fov, img_round):
        cell_background_level = np.floor(self.photon_counts[img_round] / self.signal_to_cell_ratio[img_round])

        img = np.zeros(shape=(fov.y_size, fov.x_size))

        for cell in fov.cells:
            cell_bg_img = cell.background(cell_background_level)

            pos = cell.center[:2].astype(int)

            # Determine grid positions for overlay
            n, m = cell_bg_img.shape

            rl = max(pos[1] - np.floor(n / 2), 0)

            ru = min(pos[1] + n - np.floor(n / 2), fov.y_size)

            cl = max(pos[0] - np.floor(m / 2), 0)

            cu = min(pos[0] + m - np.floor(m / 2), fov.x_size)

            tl = pos[[1, 0]] - np.array([np.floor(n / 2), np.floor(m / 2)])

            rl, ru, cl, cu = int(rl), int(ru), int(cl), int(cu)

            tl = tl.astype(int)

            img[rl:ru, cl:cu] += cell_bg_img[rl - tl[0] : ru - tl[0], cl - tl[1] : cu - tl[1]]

        # Add global background
        img += self._get_global_background_img(cell_background_level, fov, img_round)

        return img

    def _get_global_background_img(self, cell_background_level, fov, img_round):
        # Make the gaussian kernel
        std = [10, 10]

        kernel_shape = (75, 75)  # rvs odd number

        kernel = make_gaussian_2d(std, kernel_shape)

        # Create a sparse matrix that represents the background
        shape = (fov.y_size + kernel_shape[0] - 1, fov.x_size + kernel_shape[1] - 1)

        cnt = self.rng.binomial(shape[0] * shape[1], self.bg_sampling_prob)

        r = self.rng.integers(0, shape[0], size=(cnt, 1))

        c = self.rng.integers(0, shape[1], size=(cnt, 1))

        background_sparse = SparseMatrix3D(np.ones((cnt, 1)), np.hstack((r, c)), shape, subpixel=False)

        # Convolve the kernel with the background
        background = sparse_convolve2d(background_sparse, kernel)

        background_level = np.floor(cell_background_level / self.cell_bg_to_global_bg_ratio[img_round])

        background = background * background_level / (np.mean(background) + 1e-6)  # scale the image

        return background


class PhotonSimulator(object):
    """Simulate emitter photon emission"""

    def __init__(self, photon_counts, psf, well_depth, subpixel=True):
        self.photon_counts = photon_counts

        self.psf = psf

        self.subpixel = subpixel

        self.well_depth = well_depth

        self.light_dists = self._init_light_dists()

    def get_img(self, bit_num, img_round, probed_fov):
        x_dim = probed_fov.fov.x_size + self.psf.shape[1] - 1

        y_dim = probed_fov.fov.y_size + self.psf.shape[0] - 1

        z_dim = self.psf.shape[2]

        idxs = probed_fov.emitter_barcodes[:, bit_num] == 1

        num_on = sum(idxs)

        on_pos = probed_fov.emitter_positions[idxs][:, [1, 0, 2]]

        on_vals = self.light_dists[img_round].rvs(num_on) / np.max(self.psf)

        # Conv with psf
        idxs = np.floor(on_pos) + np.array([np.floor(self.psf.shape[1] / 2), np.floor(self.psf.shape[0] / 2), 0])

        idxs = idxs.astype(int)

        sparse_matrix = SparseMatrix3D(on_vals, idxs, (y_dim, x_dim, z_dim), self.subpixel)

        img = sparse_convolve3d(sparse_matrix, self.psf)

        img = np.clip(img, a_min=0, a_max=None)

        return img

    def _init_light_dists(self):
        light_dists = []

        for c in self.photon_counts:
            v = c / 20 if c >= 20 else 0.5

            light_dists.append(TruncatedNormal(0, self.well_depth, mean=c, var=v))

        return light_dists


class CameraSimulator(object):
    """Simulate camera capture"""

    def __init__(
        self, bias, dark_current, exposure_time, gain, quantum_efficiency, read_noise, rng, threshold_dark_current=True
    ):
        self.bias = bias

        self.dark_current = dark_current

        self.exposure_time = exposure_time

        self.gain = gain

        self.quantum_efficiency = quantum_efficiency

        self.read_noise = read_noise

        self.rng = rng

        self.threshold_dark_current = threshold_dark_current

    def capture_img(self, channel, img):
        mean_dark_electrons = self.dark_current * self.exposure_time[channel]

        if self.threshold_dark_current:
            low_signal_mask = img < 10

            # Add dark current noise as a Poisson process (before QE and shot noise)
            dark_current_noise = self.rng.poisson(mean_dark_electrons, size=img.shape)

            img[low_signal_mask] += dark_current_noise[low_signal_mask]

        else:
            dark_current_noise = self.rng.poisson(mean_dark_electrons, size=img.shape)

            img += dark_current_noise

        # Quantum efficiency (Electron conversion - photodetector)
        img *= self.quantum_efficiency[channel]

        # Add photon shot noise
        img = self.rng.poisson(img).astype(float)

        # Add read noise
        # img += np.clip(np.random.normal(0, self.read_noise, size=img.shape), a_min=0, a_max=None)

        img += self.rng.lognormal(0, self.read_noise, size=img.shape)

        return (img / self.gain) + self.bias
