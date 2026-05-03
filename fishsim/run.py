import yaml


def simulate(config_file, out_dir):
    with open(config_file, "r") as fh:
        config = yaml.safe_load(fh)

    config_data = config["data"]

    psf = load_psf(config)

    codebook = Codebook.from_file(config_data["codebook"])

    data_org = DataOrganisation.from_file(config_data["data_organisation"])

    probed_fov = _load_fov(codebook, config)

    sim = _load_simulator(config, data_org, psf)

    imgs = []

    for bit in codebook.bits:
        imgs.append(sim.get_bit_img(bit, probed_fov))

    imgs = np.array(imgs)

    skimage.io.imsave(out_file, imgs)
