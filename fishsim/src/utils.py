import math
import random
from re import I
import numpy as np
import numpy.linalg as LA
import pandas as pd
from scipy.stats import truncnorm
from scipy.ndimage import gaussian_filter
from typing import Callable
from pathlib import Path
import sys
from .sparse import SparseMatrix3D, sparse_convolve2d

BASE_PROJECT_DIR = Path(__file__).resolve().parents[1].absolute()


def add_gene_ids_to_codebook(codebook_filepath, save_filepath):
    # Load the codebook
    codebook = pd.read_csv(codebook_filepath, skiprows=3)

    numeric_id = []
    count = 0
    for name in codebook["name"]:
        if name[:5] == "blank":
            numeric_id.append(-1)
            continue
        numeric_id.append(count)
        count += 1

    codebook.insert(loc=1, column="numeric_id", value=numeric_id)
    codebook.to_csv(save_filepath, index=False)


def read_groundtruth(groundtruth_file: str, savepath: str):
    groundtruth = pd.read_csv(groundtruth_file)

    sorted_groundtruth = groundtruth.sort_values(by="column", ascending=True)
    sorted_groundtruth.to_csv(savepath, index=False)


def generate_barcodes(num_genes, num_bits, num_ones):
    """
    Generates barcodesMap to represent a specified number of genes
        numGenes - the number of genes you want to create barcode
        numBits  - number of bitss
        numOnes  - number of 1s in the barcode
        writecsv - flag to indicate
    """

    # create an empty dataframe to stores the barcode
    code_table = []
    # code_table = pd.DataFrame(code_table, columns=['genes', 'barcode'])

    # for each genes, dp
    for i in range(num_genes):
        barcode_char = generate_single_barcode(num_bits, num_ones)
        code_table.append([i + 1, barcode_char])

    code_table = pd.DataFrame(code_table, columns=["genes", "barcode"])
    table_name = "generated_codebook.csv"
    code_table.to_csv(table_name)


def generate_single_barcode(num_bits, num_ones):
    list_r = list(range(1, num_bits + 1))
    r = random.sample(list_r, num_ones)
    barcode = np.zeros(num_bits)

    for i in range(num_ones):
        barcode[r[i] - 1] = 1

    barcode_char = str(barcode)

    return barcode_char


def make_gaussian_2d(
    std: list,
    shape: tuple,
    rot: float = 0,
    inner_extent: float = None,
    outer_extent: float = None,
    smoothing: float = 0,
) -> np.ndarray:
    """ Creates a 2D gaussian image

        Args:
            std (list): standard deviation in the form [a, b]
            shape (tuple): shape of the returned image in the form (rows, cols)
            rot (float): rotation to be applied to the gaussian in RADs
            inner_extent (float): fraction of the standard deviation
            outer_extent (float): fraction of the standard deviation
            smoothing (float): sigma of the gaussian kernel used for smoothing

        Returns:
            gaussian image: 2D gaussian image (not scaled)

    """

    f = gaussian_2d(std, rot)
    x = np.linspace(-math.floor(shape[0] / 2), math.floor(shape[0] / 2), shape[0])
    y = np.linspace(-math.floor(shape[1] / 2), math.floor(shape[1] / 2), shape[1])
    c, r = np.meshgrid(x, y)  # changed mappings: x -> cols, y -> rows

    R = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]])
    if inner_extent is not None:
        inner_axes = np.array(std) * inner_extent
        inner_V = np.diag(1 / inner_axes ** 2)  # ellipse matrix w/o rotation
        inner_V = R @ inner_V @ R.T  # ellipse matrix w/ rotation applied
        inner_floor = f(R @ np.array([std[0] * inner_extent, 0])) * 0.4

    if outer_extent is not None:
        outer_axes = np.array(std) * outer_extent
        outer_V = np.diag(1 / outer_axes ** 2)  # ellipse matrix w/o rotation
        outer_V = R @ outer_V @ R.T  # ellipse matrix w/ rotation applied

    img = np.zeros(shape)

    for i in range(shape[0]):
        for j in range(shape[1]):
            x = np.array([r[i, j], c[i, j]])
            if inner_extent is not None and x @ inner_V @ x.T <= 1:
                img[i, j] = inner_floor
                continue
            if outer_extent is not None and x @ outer_V @ x.T > 1:
                continue
            img[i, j] = f(x)

    # Smooth the image remove sharp transitions
    if inner_extent is not None or outer_extent is not None:
        img = gaussian_filter(img, sigma=smoothing)

    return img


def gaussian_2d(std: list, rot: float) -> Callable[[np.ndarray], float]:
    """ Creates a 2D gaussian fucntion
        Args: 
            std (list): standard deviations of the gaussian [sigma_1, sigma_2]
            rot (float): rotation angle in radians
        
        Returns:
            2D gaussian function
    """
    V = np.diag(np.array(std) ** 2)  # covariance matrix w/o rotation
    R = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]])
    V = R @ V @ R.T  # cov. matrix w/ rotation applied

    det_V, inv_V = LA.det(V), LA.inv(V)  # precalculating for performance
    coeff = 1 / (np.sqrt((2 * np.pi) ** 2) * det_V)

    return lambda x: coeff * np.exp(-0.5 * (x.T @ inv_V @ x))


def glob_background(shape: tuple, sampling_prob: float, bg_lvl: float) -> np.ndarray:
    """ Creates an image with randomly sampled gaussian noise

        Args:
            shape (tuple): size of the resulting image
            sampling_prob (float): 0 ~ 1 value that indicates the probablity of each pixel contributing to the background
            bg_lvl (float): image will be scaled such that the mean is equal to the bg_lvl

        Returns:
            background_image: image with randomly sampled gaussian noise
    """
    # Make the gaussian kernel
    std = [10, 10]
    kernel_shape = (75, 75)  # choose odd number
    kernel = make_gaussian_2d(std, kernel_shape)

    # Create a sparse matrix that represents the background
    shape = (shape[0] + kernel_shape[0] - 1, shape[1] + kernel_shape[1] - 1)
    cnt = np.random.binomial(shape[0] * shape[1], sampling_prob)
    r = np.random.randint(0, shape[0], size=(cnt, 1))
    c = np.random.randint(0, shape[1], size=(cnt, 1))
    background_sparse = SparseMatrix3D(
        np.ones((cnt, 1)), np.hstack((r, c)), shape, subpixel=False
    )

    # Convolve the kernel with the background
    background = sparse_convolve2d(background_sparse, kernel)
    background = background * bg_lvl / (np.mean(background) + 1e-6)  # scale the image

    return background


class trunc_norm:
    """
        Sytactic sugar to create truncated normal distribution using scipy.stats.truncnorm
        Args:
            a (int): lower bound of distribution
            b (int): upper bound of distribution
            mean (int): mean of normal distribution
            var (int): std deviation of distribution
        Methods:
            choose()
                returns a random variable from the specified distribution
    """
    def __init__(self, a:np.uint16, b:np.uint16, mean:np.uint16=50, var:np.uint16=100):
        self.mean=mean
        self.var=var
        """
            These need to be general normal distributions
            scipy does not use the entire truncated normal distribution formula
            they plug these numbers directly into the bottom half
        """
        self.lower=(a-self.mean)/self.var
        self.upper=(b-self.mean)/self.var
        
        

    def choose(self,num:int=1) -> int:
        return truncnorm.rvs(self.lower,self.upper,loc=self.mean,scale=self.var,size=num)
    @staticmethod
    def choose_from_array(means,a:np.uint16=0,b:np.uint16=80000,var:int=20):
        return np.array([trunc_norm(a,b,mean=mean,var=mean/var).choose() for mean in means])


if __name__ == "__main__":
    pass
