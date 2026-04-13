import math
import time
import random
import yaml
from pathlib import Path

import numpy as np
import pandas as pd
import pickle

import scipy.io
import scipy.signal
import skimage.io
import matplotlib.pyplot as plt

from .generate_emitters import cell_emitter_position, random_emitter_position, enforce_min_center_distance
from .sparse import SparseMatrix3D, sparse_convolve3d
from .utils import BASE_PROJECT_DIR, glob_background, trunc_norm


class Simulator(object):
    # TODO: make a class for these inputs
    def __init__(
        self, simulation_params: dict, camera_params: dict, optics_params: dict
    ):

        # Simulation parameters
        self.number_emitter = simulation_params["emitter_count"]
        self.image_size = simulation_params["image_size"]
        self.scr = simulation_params["scr"]
        self.photon_count = simulation_params["photon_count"]
        
        # Handle scalar vs list for round-level variation
        self.photon_count_list = (
            self.photon_count if isinstance(self.photon_count, list) else None
            )
        self.scr_list = self.scr if isinstance(self.scr, list) else None
        self.global_background_level_list = None
        self.wavelengths = simulation_params["signal_wavelengths"]
        self.bit_drop = simulation_params["bitdrop_probability"]
        self.bit_add = simulation_params["bitadd_probability"]
        self.bg_sampling_prob = simulation_params["background_sampling_probability"]

        # Cell parameters
        self.cell_count = simulation_params["cells"]["count"]
        self.cell_axes = simulation_params["cells"]["axes"]
        self.cells = []  # list to store the cell instances

        # Camera parameters
        self.QE = camera_params["qe"]
        self.gain = camera_params["gain"]
        self.bias = camera_params["bias"]
        self.dark_current = camera_params["dark_current"]
        self.exposure_times = camera_params["exposure_times"]
        self.read_noise = camera_params["read_noise"]
        self.sensor_size = camera_params["sensor_size"]
        self.well_depth = camera_params["well_depth"]

        # Optical parameters
        self.magnification = optics_params["magnification"]

    def generate_training_data(
        self,
        codebook_filepath: Path,
        psf_filepath: Path,
        data_org_filepath: Path,
        tiles: int = 1,
        is_matlab: bool = True,
        is_subpixel: bool = True,
        is_cell: bool = True,
        is_nucleus: bool = False,
        z_bound: tuple = None,
        savepath: Path = None,
    ) -> None:
        """
        This funtion generate the training data for number of set of merFISH data you want. Called generate_single_traing_data tiles times
        Inputs:  codebook_filename: a csv file containing the merfish library you are using
                 psf_filename: the file containing the psf of your microscope (in resource file if you look at our data)
                 tiles: number of set of data(16 images per set) you want
                 is_matlab: boolean value, true to take psf from .mat file, otherwise use python pickle psf file
                 num_bits: number of barcodes in your library

        Output:  groundtruth.csv: contain the photon count, location, gene and barcode of every beads
                 Tiff images: the tiff image for each bit, total num_bits
        """

        # Load the codebook
        codebook = pd.read_csv(codebook_filepath, skiprows=3)
        data_organization = pd.read_csv(data_org_filepath)
        # TODO: count number of bits in the codebook and remove the argument
        
        # Map each bit to its imaging round (0-based indexing)
        self.bit_to_round = {
            row["bitNumber"] - 1: row["imagingRound"]
            for _, row in data_organization.iterrows()
        }

        # Load psf
        if is_matlab:
            mat = scipy.io.loadmat(str(psf_filepath))
            psf = np.array(mat["ans"])
        else:
            with open(str(psf_filepath), "rb") as fp:
                psf_data = pickle.load(fp)
            psf = psf_data["psf"]
            psf = psf[:, 0 : psf.shape[1] - 1, 0 : psf.shape[2] - 1]
            psf = np.transpose(psf, (1, 2, 0))

        # Normalize the psf and offset the values such that edges are equal to zero
        psf = Simulator.process_psf(psf)
        z_bound = (0, psf.shape[2]) if z_bound is None else z_bound # PSF is 3D, (29,29,1001)

        t = time.strftime("%Y%m%d-%H-%M-%S")
        folder_name = "generated_data_merfish_" + t
        if savepath is None:
            folder_path = BASE_PROJECT_DIR / "results" / folder_name
        else:
            folder_path = savepath / folder_name
        folder_path.mkdir(parents=True, exist_ok=True)

        for i in range(tiles):
            print("starting tile " + str(i + 1))
            self.__generate_single_training_data(
                codebook,
                data_organization,
                psf,
                i,
                is_subpixel,
                is_cell,
                is_nucleus,
                z_bound,
                folder_path,
            )
            
        # Log simulation parameters
        self.print_simulation_params(
            folder_path,
            codebook_filepath,
            psf_filepath,
            data_org_filepath,
            tiles,
            is_matlab,
            is_subpixel,
            is_cell,
            is_nucleus,
        )

    def print_simulation_params(
        self,
        save_path: Path,
        codebook_filepath: Path,
        psf_filepath: Path,
        data_org_filepath: Path,
        tiles: int = 1,
        is_matlab: bool = True,
        is_subpixel: bool = True,
        is_cell: bool = True,
        is_nucleus: bool = False,
    ) -> None:
        # Structuring the parameters in the same format as the config file
        config = {
            "simulation": {
                "tile_count": tiles,
                "emitter_count": self.number_emitter,
                "image_size": self.image_size,
                "signal_wavelengths": self.wavelengths,
                "scr": self.scr,
                "photon_count": self.photon_count,
                "photon_count_list": self.photon_count_list, # returns None if photon_count is a single value
                "scr_list": self.scr_list, # returns None if scr is a single value
                "bitdrop_probability": self.bit_drop,
                "bitadd_probability": self.bit_add,
                "is_cell": is_cell,
                "is_matlab": is_matlab,
                "is_subpixel": is_subpixel,
                "is_nucleus": is_nucleus,
                "background_sampling_probability": self.bg_sampling_prob,
                "global_background_level_list": self.global_background_level_list,
                "cells": {"count": self.cell_count, "axes": self.cell_axes},
            },
            "data": {
                "psf": {"psf": str(psf_filepath)},
                "codebook": {"codebook": str(codebook_filepath)},
                "data_organization": str(data_org_filepath),
            },
            "cameras": {
                "camera": {
                    "qe": self.QE,
                    "gain": self.gain,
                    "bias": self.bias,
                    "dark_current": self.dark_current,
                    "exposure_times": self.exposure_times,
                    "read_noise": self.read_noise,
                    "sensor_size": self.sensor_size,
                }
            },
            "optics": {"microscope": {"magnification": self.magnification}},
        }
        with open(str(save_path / "config.yml"), "w") as outfile:
            yaml.dump(config, outfile, default_flow_style=False)

    def __generate_single_training_data(
        self,
        codebook: pd.DataFrame,
        data_organization: pd.DataFrame,
        psf: np.ndarray,
        tile_num: int,
        is_subpixel: bool,
        is_cell: bool,
        is_nucleus: bool,
        z_bound: tuple,
        folder_path: Path,
    ) -> None:
        """Creates a single set of training data

        Args:
            codebook (pd.DataFrame): codebook dataframe
            data_organization (pd.DataFrame): data organization dataframe
            psf (np.ndarray): psf matrix
            tile_num (int): number of the tile that's being generated
            is_subpixel (bool): if true, emitter locations are floating point (i.e. non-zero shift)
            is_cell (bool): if true, cell structures imposed on the simulation
            is_nucleus (bool): if true, cells have a nucleus
            z_bound (tuple): z range for the emitter/cell center position (centered around zero)
            folder_path (Path): folder where the results are saved
        """
        # Simulate flourescent probes
        genes = codebook["name"]
        gene_ids = codebook["numeric_id"]
        barcodes = codebook["barcode"]
        numofbits = len(barcodes[0].replace(" ", ""))

        # if there is distribution column, draw sample based on that distribution
        if "distribution" in codebook:
            gene_index = list(range(len(genes)))
            distribution = np.array(codebook["distribution"].fillna(0).to_list())
            distribution = [int(item) for item in distribution]

            # r contain the number of lines in the table for the genes
            r = [i for i in gene_index for _ in range(distribution[i])]
            random.shuffle(r)  # shuffle to randomize the order of genes
            num_emitter = len(r)
        
        # if not, draw sample using uniform distribution
        else:
            genes = codebook["name"]
            # r contain the number of lines in the table for the genes
            r = np.random.randint(0, np.size(barcodes), size=self.number_emitter)
            num_emitter = self.number_emitter

        if is_cell:
            emitter_pos, self.cells = cell_emitter_position(
                (0, self.image_size),
                (0, self.image_size),
                (
                    math.floor(psf.shape[2] / 2) + z_bound[0],
                    math.floor(psf.shape[2] / 2) + z_bound[1],
                ),
                num_emitter,
                self.cell_count,
                self.cell_axes,
                is_nucleus,
            )
        else:
            # emitter_pos has shape (num_emitters, 3 columns)
            # col 1: x, col 2: y, col 3: z; see generate_emitters.py
            emitter_pos = random_emitter_position(
                (0, self.image_size),
                (0, self.image_size),
                (
                    math.floor(psf.shape[2] / 2) + z_bound[0],
                    math.floor(psf.shape[2] / 2) + z_bound[1],
                ),
                num_emitter,
            )

        if not is_subpixel:
            emitter_pos = np.floor(emitter_pos)

        # Enforce spacing to prevent transcripts overlapping (on the same X,Y coordinates)
        emitter_pos = enforce_min_center_distance(
            emitter_pos,
            min_dist=6.0,
            x_dim=(0, self.image_size),
            y_dim=(0, self.image_size),
            z_dim=(
                math.floor(psf.shape[2] / 2) + z_bound[0],
                math.floor(psf.shape[2] / 2) + z_bound[1],
            ),
        )

        # Crate a goundtruth data frame
        ground_truth = []

        # emitter_barcodes store only the barcode of emitter
        emitter_barcodes = np.empty((emitter_pos.shape[0], numofbits))

        # Randomly determine the emitters that will have bit drop/add
        is_bit_drop = np.random.rand(emitter_pos.shape[0]) <= self.bit_drop
        is_bit_add = np.random.rand(emitter_pos.shape[0]) <= self.bit_add
        #distribution for the emitter brightness
        
        # Precompute 8 light_dist objects (one per imaging round)
        num_rounds = len(set(self.bit_to_round.values()))  # Should be 8
        light_dists_by_round = {}
        
        # This for-loop will create a light distribution per imaging round
        # For 8 imagingRounds, there will be 8 light distributions.
        for round_idx in range(num_rounds):
            photon_count = (
                self.photon_count_list[round_idx]
                if self.photon_count_list
                else self.photon_count
            )
            light_dists_by_round[round_idx] = trunc_norm(
                0,
                self.well_depth,
                mean=photon_count,
                var=int(photon_count / 20) if photon_count >= 20 else 0.5
            )

        for i in range(emitter_pos.shape[0]):
            gene = genes[r[i]]
            gene_id = gene_ids[r[i]]
            barcode = barcodes[r[i]].replace(" ", "")
            mapped_barcode = Simulator.map_barcode(barcode, data_organization)
            ground_truth.append(
                [   
                    # Generate one photon count for each of the 8 imaging rounds
                    # The original code for this was light_dist.choose(), and it serves no particular purpose.
                    # We follow its original purpose of "generating some light_dist value", here we generate
                    # a list of 8 values assuming 8 imagingRounds.
                    [
                        light_dists_by_round[round_idx].choose()[0]
                        for round_idx in range(8)
                    ],
                    int(
                        math.floor(emitter_pos[i, 1])
                        + self.image_size * math.floor(emitter_pos[i, 0])
                    ),
                    math.floor(emitter_pos[i, 1]),
                    math.floor(emitter_pos[i, 0]),
                    math.floor(emitter_pos[i, 2]),
                    f"{(emitter_pos[i, 1] % 1):5f}",
                    f"{(emitter_pos[i, 0] % 1):5f}",
                    f"{(emitter_pos[i, 2] % 1):5f}",
                    gene,
                    gene_id,
                    f"'{barcode}'",
                    f"'{mapped_barcode}'",
                    is_bit_drop[i],
                    is_bit_add[i],
                ]
            )
            emitter_barcodes[i, :] = list(mapped_barcode)
            if is_bit_drop[i]:
                index_ones = [
                    k for k in range(numofbits) if emitter_barcodes[i, k] == 1
                ]
                dropout_bit = random.choice(index_ones)
                emitter_barcodes[i, dropout_bit] = 0

            if is_bit_add[i]:
                index_zeros = [
                    k for k in range(numofbits) if emitter_barcodes[i, k] == 0
                ]
                add_bit = random.choice(index_zeros)
                emitter_barcodes[i, add_bit] = 1

        # Create a new folder store the image and groundtruth
        tile_folder_name = "tile_" + str(tile_num + 1)
        tile_folder_path = folder_path / tile_folder_name
        tile_folder_path.mkdir(parents=True, exist_ok=True)

        # Create folders for wavelength
        for wavelength in self.wavelengths:
            (tile_folder_path / str(wavelength)).mkdir(parents=True, exist_ok=True)

        # Save the ground truth table
        ground_truth = pd.DataFrame(
            ground_truth,
            columns=[
                "photoncount",
                "pixel_index_num",
                "column",
                "row",
                "z",
                "column_shift",
                "row_shift",
                "z_shift",
                "genes",
                "gene_id",
                "barcode",
                "mapped_barcode",
                "is_bit_drop",
                "is_bit_add",
            ],
        )
        ground_truth.to_csv(tile_folder_path / "groundtruth.csv", index=False)

        # Save the frequency of a gene appers
        frequency = (
            ground_truth.value_counts("genes")
            .rename_axis("genes")
            .to_frame(name="counts")
        )
        frequency.to_csv(tile_folder_path / "frequency.csv")

        # Main loop

        # Grid dimension
        x_dim = self.image_size + psf.shape[1] - 1
        y_dim = self.image_size + psf.shape[0] - 1
        z_dim = psf.shape[2]

        # Creating cellular background image
        bg_lvl_list = np.floor(list(np.array(self.photon_count_list)/np.array(self.scr_list)))
        bg_imgs_by_round = {}
        global_background_level_list = []

        for round_idx in range(num_rounds):
            
            bg_lvl = math.floor(bg_lvl_list[round_idx])
            bg_img = np.zeros(shape=(self.image_size, self.image_size))
            
            # If is_cell: False, then this for-loop will still run, but nothing happens.
            for cell in self.cells:
                cell_bg_img = cell.background(bg_lvl)
                pos = cell.center[:2].astype(int)
    
                # Determine grid positions for overlay
                n, m = cell_bg_img.shape
                rl = pos[0] - math.floor(n / 2)
                rl = rl if rl >= 0 else 0
                ru = pos[0] + n - math.floor(n / 2)
                ru = ru if ru <= self.image_size else self.image_size
                cl = pos[1] - math.floor(m / 2)
                cl = cl if cl >= 0 else 0
                cu = pos[1] + m - math.floor(m / 2)
                cu = cu if cu <= self.image_size else self.image_size
    
                tl = pos - np.array([math.floor(n / 2), math.floor(m / 2)])
                bg_img[rl:ru, cl:cu] += cell_bg_img[
                    rl - tl[0] : ru - tl[0], cl - tl[1] : cu - tl[1]
                ]

            # Add global background
            cell_bg_to_global_bg_ratio = random.uniform(1.15, 2.5)
            glob_bg_lvl = math.floor(bg_lvl / cell_bg_to_global_bg_ratio)
            global_background_level_list.append(glob_bg_lvl)
            bg_img += glob_background(
                shape=(self.image_size, self.image_size),
                sampling_prob=self.bg_sampling_prob,
                bg_lvl=glob_bg_lvl,
            )
            
            # Store the final per-round background image
            bg_imgs_by_round[round_idx] = bg_img
            
        self.global_background_level_list = global_background_level_list

        for i in range(numofbits):
            # Let's fetch the round index, given the bit number i.
            round_idx = self.bit_to_round[i]  # 0 to 7
            light_dist = light_dists_by_round[round_idx]
            bg_img = bg_imgs_by_round[round_idx]

            background_matrix = np.zeros(shape=(x_dim, y_dim, z_dim))
            place_holder = emitter_barcodes[:, i].astype(bool)
            on_emitter_pos = emitter_pos[place_holder.astype(bool), :]
            on_emitter_val = (
                light_dist.choose(np.array(ground_truth.loc[:, "photoncount"][place_holder]).shape[0]) / np.max(psf)   
            ).astype(int)
            wavelength_name = str(self.wavelengths[i % len(self.wavelengths)])
            QE_used = self.QE[
                str(self.wavelengths[i % len(self.wavelengths)])
            ]  # choose the QE and the filename
            # create the emitters
            for j in range(np.shape(on_emitter_pos)[0]):
                background_matrix[
                    math.floor(on_emitter_pos[j, 0]) + math.floor(psf.shape[1] / 2),
                    math.floor(on_emitter_pos[j, 1]) + math.floor(psf.shape[0] / 2),
                    math.floor(on_emitter_pos[j, 2]),
                ] = on_emitter_val[j]

            # Conv with psf
            indexes = np.floor(on_emitter_pos).astype(int) + np.array(
                [math.floor(psf.shape[1] / 2), math.floor(psf.shape[1] / 2), 0]
            )
            #print(np.max(background_matrix))
            sparse_matrix = SparseMatrix3D(
                on_emitter_val, indexes, (x_dim, y_dim, z_dim), is_subpixel
            )
            
            photon_im = sparse_convolve3d(sparse_matrix, psf)
            photon_im = photon_im.reshape((photon_im.shape[0], photon_im.shape[1], 1))
            #print(np.max(photon_im))
            photon_im[:, :, 0] += bg_img
            #print(np.max(photon_im))
            # photon_im = scipy.signal.convolve(background_matrix, psf, mode="valid")
            photon_im[photon_im < 0] = 0
            
            # Simulating dark current noise
            # Given the correct channel_wavelength, we have different exposure time
            exposure_time = self.exposure_times[wavelength_name]
             
            # Dark current has units electrons/pixel/second, multiplied by exposure time to get mean_dark_electrons
            mean_dark_electrons = self.dark_current * exposure_time
            
            # Optionally threshold for low-signal areas (< 10 photon counts)
            self.threshold_dark_current = True
             
            if self.threshold_dark_current:
                low_signal_mask = photon_im < 10
                 
                print(low_signal_mask.sum())
                 
                # Add dark current noise as a Poisson process (before QE and shot noise)
                dark_current_noise = np.random.poisson(mean_dark_electrons, size=photon_im.shape)
                photon_im[low_signal_mask] += dark_current_noise[low_signal_mask]
            else:
                dark_current_noise = np.random.poisson(mean_dark_electrons, size=photon_im.shape)
                photon_im += dark_current_noise

            # Quantum efficiency (Electron conversion - photodetector)
            electron_im = np.array(photon_im) * QE_used

            # Add photon shot noise
            electron_im = [*map(np.random.poisson, electron_im)]
            #print(np.max(electron_im))
            row = np.shape(electron_im)[0]
            col = np.shape(electron_im)[1]
            read_noise_im = np.random.normal(0, self.read_noise, (row, col, 1))

            read_noise_im[read_noise_im < 0] = 0

            electron_im2 = np.array(electron_im) + np.array(read_noise_im)

            # converting electrons to ADUs
            digital_im = electron_im2 / self.gain + self.bias

            # save image
            file_name = (
                "merFISH_" + f"{i+1:02}" + "_" + f"{tile_num+1:03}" + "_01" + ".tiff"
            )
            skimage.io.imsave(
                str(tile_folder_path / wavelength_name / file_name),
                digital_im.astype(np.uint16),
            )
            print("Round %i completed \n" % (i + 1))

    @staticmethod
    def process_psf(psf: np.ndarray) -> np.ndarray:
        psf = psf.astype(np.float64)
        psf = psf / np.sum(psf)  # normalize the psf

        # Create binary mask for the psf
        radius = 14 / 2
        mask_x = np.arange(-np.ceil(psf.shape[0] / 2 - 1), np.ceil(psf.shape[0] / 2), 1)
        mask_y = np.arange(-np.ceil(psf.shape[1] / 2 - 1), np.ceil(psf.shape[1] / 2), 1)
        xv, yv = np.meshgrid(mask_x, mask_y)

        circle_mask = xv ** 2 + yv ** 2 >= radius ** 2

        # Apply the binary mask to the entire psf, substract average of the edge value
        for i in range(psf.shape[2]):
            zslice = psf[:, :, i]
            z_masked = zslice * circle_mask.astype(int)
            mean_edge = np.sum(z_masked) / np.count_nonzero(circle_mask)
            zslice = zslice - mean_edge
            zslice[zslice < 0] = 0
            psf[:, :, i] = zslice

        return psf

    @staticmethod
    def map_barcode(barcode: str, data_organization: pd.DataFrame) -> str:
        bit_number = data_organization["bitNumber"]
        imaging_round = data_organization["imagingRound"]
        color = data_organization["color"]

        mapping = dict()
        for i in range(data_organization.shape[0]):
            mapping[bit_number[i] - 1] = (
                imaging_round[i] * 2 if color[i] == 650 else imaging_round[i] * 2 + 1
            )

        mapped_barcode = ["0"] * len(barcode)
        for i, bit in enumerate(barcode):
            mapped_barcode[mapping[i]] = bit

        return "".join(mapped_barcode)


if __name__ == "__main__":
    pass