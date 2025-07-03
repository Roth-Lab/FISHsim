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
    # # #generate a single ellipsoid
    e1 = Ellipsoid([4, 0, 0], [2, 2, 2])
    cell_axes_bounds = {"a": [20, 20], "b": [60, 60], "c": [30, 30]}
    cell = EllipsoidCell([100, 0, 0], cell_axes_bounds,[90,30,10])
    cell.generate_emitters(1000)
    # Plot 3D cell structure
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.set_box_aspect([1,1,1])
    ax.set_xlim3d(-100, 100)
    ax.set_ylim3d(-100, 100)
    ax.set_zlim3d(-100, 100)
    ax.scatter3D(cell.emitters[:, 0], cell.emitters[:, 1], cell.emitters[:, 2])
    plt.show()

#Generate two ellpsoid

   


