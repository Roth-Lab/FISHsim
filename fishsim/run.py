import numpy as np
import pandas as pd
import skimage.io
import yaml

from fishsim.codebook import Codebook
from fishsim.data_organisation import DataOrganisation
from fishsim.fov import FieldOfView, ProbedFieldOfView
from fishsim.psf import load_psf
from fishsim.simulate import BackgroundSimulator, CameraSimulator, ImageSimulator, PhotonSimulator


def simulate(codebook_file, config_file, data_org_file, img_file, seed=None):
    rng = np.random.default_rng(seed)

    with open(config_file, "r") as fh:
        config = yaml.safe_load(fh)

    config_camera = config["camera"]

    config_sim = config["simulation"]

    codebook, data_org = _load_codebook_data_org(codebook_file, data_org_file)

    psf = load_psf()

    boundary_box = _load_boundary_box(config_sim, psf)

    config_cell = config_sim["cells"]

    fov = FieldOfView(
        boundary_box,
        config_cell["axes"],
        config_cell["num_cells"],
        rng,
    )

    probed_fov = ProbedFieldOfView(
        codebook,
        fov,
        config_sim["num_emitters"],
        rng,
        bit_add_prob=config_sim["bit_add_prob"],
        bit_drop_prob=config_sim["bit_drop_prob"],
        sim_nucleus=config_cell["nucleus"],
        subpixel=config_sim["subpixel"],
    )

    sim_background = BackgroundSimulator(
        config_sim["bg_sampling_prob"],
        config_sim["photon_counts"],
        rng,
        config_sim["signal_to_cell_ratio"],
    )

    sim_photon = PhotonSimulator(
        config_sim["photon_counts"],
        psf,
        config_camera["well_depth"],
        config_sim["subpixel"],
    )

    sim_camera = CameraSimulator(
        config_camera["bias"],
        config_camera["dark_current"],
        config_camera["exposure_time"],
        config_camera["gain"],
        config_camera["quantum_efficiency"],
        config_camera["read_noise"],
        rng,
    )

    sim = ImageSimulator(data_org, sim_background, sim_camera, sim_photon)

    imgs = []

    for b in data_org.bit_ids:
        imgs.append(sim.get_bit_img(b, probed_fov))

    imgs = np.array(imgs, dtype=np.uint16)

    skimage.io.imsave(img_file, imgs)


def _load_boundary_box(config_sim, psf):
    boundary_box = config_sim["boundary_box"]

    boundary_box["z"] = (
        boundary_box["z"][0] + np.floor(psf.shape[2] // 2),
        boundary_box["z"][1] + np.floor(psf.shape[2] // 2),
    )

    return boundary_box


def _load_codebook_data_org(codebook_file, data_org_file):
    codebook_df = pd.read_csv(codebook_file, index_col="target", sep="\t")

    data_org_df = pd.read_csv(data_org_file, converters={"bit_number": int}, sep="\t")

    data_org = DataOrganisation(data_org_df)

    codebook_df = codebook_df[data_org.bit_ids]

    codebook = Codebook(codebook_df)

    return codebook, data_org
