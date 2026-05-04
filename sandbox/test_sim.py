import matplotlib.pyplot as pp
import numpy as np
import pandas as pd
import scipy.io
import skimage.io

from fishsim.codebook import Codebook
from fishsim.data_organisation import DataOrganisation
from fishsim.fov import FieldOfView, ProbedFieldOfView
from fishsim.psf import _process_psf
from fishsim.simulate import BackgroundSimulator, CameraSimulator, PhotonSimulator

rng = np.random.default_rng(0)

psf_file = "../fishsim/resources/psf.mat"

psf = np.array(scipy.io.loadmat(psf_file)["ans"])

c = psf.shape[2] // 2

# pp.imshow(psf[:, :, c])
#
# pp.show()
#
psf = _process_psf(psf)

# pp.imshow(psf[:, :, c])
#
# pp.show()

boundary_box = {"x": (0, 400), "y": (0, 400), "z": (-7, 7)}

cell_axes_bounds = {"a": (30, 40), "b": (60, 70), "c": (50, 60)}

boundary_box_psf = boundary_box.copy()

boundary_box_psf["z"] = (
    boundary_box_psf["z"][0] + np.floor(psf.shape[2] // 2),
    boundary_box_psf["z"][1] + np.floor(psf.shape[2] // 2),
)

fov = FieldOfView(boundary_box_psf, cell_axes_bounds, 15, rng)

codebook_df = pd.DataFrame(
    [
        {"target": "gene_1", "bit_1": 1, "bit_2": 0, "bit_3": 0, "bit_4": 0},
        {"target": "gene_2", "bit_1": 0, "bit_2": 1, "bit_3": 0, "bit_4": 0},
        {"target": "gene_3", "bit_1": 0, "bit_2": 0, "bit_3": 1, "bit_4": 0},
        {"target": "gene_4", "bit_1": 0, "bit_2": 0, "bit_3": 0, "bit_4": 1},
    ]
)

codebook = Codebook(codebook_df)

probed_fov = ProbedFieldOfView(
    codebook,
    fov,
    5000,
    rng,
    bit_add_prob=0.1,
    bit_drop_prob=0.1,
    sim_nucleus=False,
    subpixel=False,
)

data_org_df = pd.DataFrame(
    [
        {"bit_id": "bit_1", "bit_number": 1, "channel": "650", "imaging_round": 0},
        {"bit_id": "bit_2", "bit_number": 2, "channel": "750", "imaging_round": 0},
        {"bit_id": "bit_3", "bit_number": 3, "channel": "650", "imaging_round": 1},
        {"bit_id": "bit_4", "bit_number": 4, "channel": "750", "imaging_round": 1},
    ]
)

data_org = DataOrganisation(data_org_df)

assert np.all(data_org.bit_ids == codebook.bit_ids)

print(data_org.get_channel("bit_2"))

print(data_org.get_round("bit_3"))

true_barcodes = np.array([codebook.get_barcode(t) for t in probed_fov.emitter_target_ids])

diff_barcodes = np.sum(probed_fov.emitter_barcodes - true_barcodes, axis=1)

df = pd.DataFrame(
    {
        "cell_id": probed_fov.emitter_cell_ids,
        "target": probed_fov.emitter_targets,
        "x": probed_fov.emitter_positions[:, 1],
        "y": probed_fov.emitter_positions[:, 0],
        "z": probed_fov.emitter_positions[:, 2],
        "bit_add": (diff_barcodes == -1).astype(int),
        "bit_drop": (diff_barcodes == 1).astype(int),
    }
)

print(df)

photon_counts = 20000 * np.ones(data_org.num_rounds)

sim_bg = BackgroundSimulator(0.00, photon_counts, rng, 11 * np.ones(data_org.num_rounds))

img = sim_bg.get_img(fov, 0)

pp.imshow(img)

pp.show()


# psf_file = "../fishsim/resources/psf"
#
# with open(str(psf_file), "rb") as fp:
#     psf_data = pickle.load(fp)
# psf = psf_data["psf"]
# psf = psf[:, 0 : psf.shape[1] - 1, 0 : psf.shape[2] - 1]
# psf = np.transpose(psf, (1, 2, 0))

# pp.imshow(psf[:, :, 0])
#
# pp.show()

sim_photon = PhotonSimulator(photon_counts, psf, 80000, subpixel=False)

img += sim_photon.get_img(0, 0, probed_fov)

pp.imshow(img)

pp.show()

sim_camera = CameraSimulator(
    100,
    0,
    {"650": 2.0, "750": 2.0},
    1.33,
    {"650": 0.89, "750": 0.71},
    1,
    rng,
)

img = sim_camera.capture_img(data_org.get_channel("bit_1"), img)


# pp.savefig("/home/andrew/Desktop/sim.png", dpi=1200)

skimage.io.imsave("/home/andrew/Desktop/sim.tif", img.astype(np.uint16))
