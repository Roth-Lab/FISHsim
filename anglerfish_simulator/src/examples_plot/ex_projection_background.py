import math
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
from ellipsoid import Ellipsoid
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from cells import EllipsoidCell

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d

    fig, axes = plt.subplots(2, 2)


    cell_axes_bounds = {"a": [80, 80], "b": [80, 80], "c": [80, 80]}
    cell = EllipsoidCell([100, 0, 100], cell_axes_bounds, [0, 0, 90])
    cell.generate_emitters(1000)
    
    cell2 = EllipsoidCell([0, 0, 100], cell_axes_bounds, [0, 0, 90])
    cell2.generate_emitters(1000)
    gaussian_img = cell2.background(10)
    # print(gaussian_img.shape)
    # print(cell2.center)

    axes[0,0].plot(cell2.emitters[:, 0]+160,cell2.emitters[:, 1]+160,'.')
    axes[0,0].set_title("Projected cell emitters")
    axes[0,0].set_aspect('equal','box')
    axes[0,0].set(xlim=(0, 320), ylim=(0, 320))

    axes[0,1].imshow(gaussian_img, cmap = "gray")
    axes[0,1].set_title("Projected cell background")
    axes[1,0].plot(np.arange(0,320,1),gaussian_img[:,160])
    axes[1,0].set_ylim([0,10])

    axes[1,0].set_title("Gaussian cross-section y")
    axes[1,1].plot(np.arange(0,320,1),gaussian_img[160,:])
    axes[1,1].set_ylim([0,10])
    axes[1,1].set_title("Gaussian cross-section x")
    fig.tight_layout()
    plt.show()