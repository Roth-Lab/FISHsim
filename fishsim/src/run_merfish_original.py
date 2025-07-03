import os
import math
import yaml
import time
import re

import shutil
import numpy as np
import scipy
import skimage

from utils import BASE_PROJECT_DIR
from simulation import Simulator
import sparse


# Parse the configuration file
with open(str(BASE_PROJECT_DIR / "resources" / "config.yml"), "r") as f:
    config = yaml.safe_load(f)

# Choose the appropriate camera and optics parameters
simulator = Simulator(
    config["simulation"], config["cameras"]["Prime95B"], config["optics"]["NikonTi2"],
)

# Choose the appropriate psf and codebook file
psf_file = BASE_PROJECT_DIR / config["data"]["psf"]["matlab"]
codebook_file = BASE_PROJECT_DIR / config["data"]["codebook"]["c1e1"]
data_org_file = BASE_PROJECT_DIR / config["data"]["data_organization"]
t = time.strftime("%Y%m%d-%H-%M-%S")
savepath = BASE_PROJECT_DIR / "results" / f"merfish_{t}"


simulator.generate_training_data(
    codebook_file,
    psf_file,
    data_org_file,
    config["simulation"]["tile_count"],
    config["simulation"]["is_matlab"],
    config["simulation"]["is_subpixel"],
    config["simulation"]["is_cell"],
    config["simulation"]["is_nucleus"],
    config["simulation"]["z_bound"],
    savepath,
)

# Set-up the folder structure
experiment_folder = savepath / "1"
experiment_folder.mkdir(parents=True, exist_ok=True)

img_dir_path_473 = experiment_folder / "473nm, Raw"
img_dir_path_473.mkdir(parents=True, exist_ok=True)
img_dir_path_561 = experiment_folder / "561nm, Raw"
img_dir_path_561.mkdir(parents=True, exist_ok=True)
img_dir_path_647 = experiment_folder / "647nm, Raw"
img_dir_path_647.mkdir(parents=True, exist_ok=True)
img_dir_path_750 = experiment_folder / "750nm, Raw"
img_dir_path_750.mkdir(parents=True, exist_ok=True)
groundtruth_path = experiment_folder / "groundtruths"
groundtruth_path.mkdir(parents=True, exist_ok=True)

# Move and rename files appropriately
img_paths = sorted(savepath.glob("generated_data_merfish_*/tile_*/*/merFISH_*"))
groundtruth_paths = sorted(
    savepath.glob("generated_data_merfish_*/tile_*/groundtruth.csv")
)
frequency_paths = sorted(savepath.glob("generated_data_merfish_*/tile_*/frequency.csv"))

img_regex = r".*tile_(\d+)\D(\d{3})\DmerFISH_(\d{2})_(\d+)_(\d{2}).*"
groundtruth_regex = r".*tile_(\d+)\Dgroundtruth.*"
frequency_regex = r".*tile_(\d+)\Dfrequency.*"

for path in img_paths:
    result = re.search(img_regex, str(path))
    round_num = math.ceil(int(result.group(3)) / 2)
    path.rename(
        experiment_folder
        / f"{result.group(2)}nm, Raw"
        / f"merFISH_{round_num:02}_{int(result.group(1)):03}_01.tiff"
    )

for path in groundtruth_paths:
    result = re.search(groundtruth_regex, str(path))
    path.rename(groundtruth_path / f"groundtruth_{result.group(1)}.csv")

for path in frequency_paths:
    result = re.search(frequency_regex, str(path))
    path.rename(groundtruth_path / f"frequency_{result.group(1)}.csv")


# Add 473 and 561 images
mat = scipy.io.loadmat(str(psf_file))
psf = np.array(mat["ans"])
psf = Simulator.process_psf(psf)
max_psf = np.max(psf)

tile_count = config["simulation"]["tile_count"]
img_size = config["simulation"]["image_size"] + psf.shape[0] - 1
photon_count = config["simulation"]["photon_count"]

values = np.array(
    [photon_count / max_psf, photon_count / max_psf, photon_count / max_psf]
).astype(int)
indices = np.array(
    [
        [
            math.ceil(img_size * 3 / 4),
            math.ceil(img_size / 2),
            math.floor(psf.shape[2] / 2),
        ],
        [
            math.ceil(img_size / 2),
            math.ceil(img_size / 2),
            math.floor(psf.shape[2] / 2),
        ],
        [
            math.ceil(img_size * 1 / 4),
            math.ceil(img_size / 2),
            math.floor(psf.shape[2] / 2),
        ],
    ]
)

sparse_matrix = sparse.SparseMatrix3D(
    values, indices, (img_size, img_size, psf.shape[2])
)
bead_im = sparse.sparse_convolve3d(sparse_matrix, psf)
bead_im = bead_im.reshape((bead_im.shape[0], bead_im.shape[1], 1))

for tile_num in range(1, tile_count + 1):
    for round_num in range(1, 12):
        skimage.io.imsave(
            str(img_dir_path_473 / f"merFISH_{round_num:02}_{tile_num:03}_01.tiff"),
            bead_im.astype(np.uint16),
        )
        skimage.io.imsave(
            str(img_dir_path_561 / f"merFISH_{round_num:02}_{tile_num:03}_01.tiff"),
            bead_im.astype(np.uint16),
        )

# Move config file
config_path = sorted(savepath.glob("generate*/config.yml"))[0]
config_path.rename(experiment_folder / "config.yml")

# Clean up
config_path = sorted(savepath.glob("generate*"))[0]
shutil.rmtree(str(config_path))
