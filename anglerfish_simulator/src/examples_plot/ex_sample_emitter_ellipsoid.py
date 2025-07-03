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
    
    cell_axes_bounds = {"a": [80, 80], "b": [80, 80], "c": [80, 80]}
    cell = EllipsoidCell([100, 0, 100], cell_axes_bounds, [0, 0, 90])
    cell.generate_emitters(1000)

   # Plot 3D cell structure
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.scatter3D(cell.emitters[:, 0], cell.emitters[:, 1], cell.emitters[:, 2])

    plt.show()