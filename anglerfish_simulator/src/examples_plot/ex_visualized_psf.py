import math
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
from ellipsoid import Ellipsoid
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from cells import EllipsoidCell
import simulation

from pathlib import Path

import scipy.io
import scipy.signal
import skimage.io

import sparse
import utils
from utils import BASE_PROJECT_DIR

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d

    #grab the psf
    mat = scipy.io.loadmat(str("resources/PSF.mat"))
    psf = np.array(mat["ans"])
    psf1 = simulation.Simulator.process_psf(psf)
    print(psf1.shape)


    #generate a matrix with shape
    matrix = np.zeros((100,200,psf1.shape[2]))
    matrix[50, 40,int(np.round(psf1.shape[2]/2) -350)]=1000000
    matrix[50, 70,int(np.round(psf1.shape[2]/2) -300)]=1000000
    matrix[50, 100,int(np.round(psf1.shape[2]/2) -250)]=1000000
    matrix[50, 130,int(np.round(psf1.shape[2]/2))]=1000000
    matrix[50, 160,int(np.round(psf1.shape[2]/2) + 150)]=1000000
    
   
    
    photon_im = scipy.signal.convolve(matrix, psf1, "valid")
    plt.figure(figsize=(5,10))
    plt.imshow(photon_im)
    plt.show()

