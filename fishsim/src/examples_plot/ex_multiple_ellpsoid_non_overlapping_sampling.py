import math
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
from ellipsoid import Ellipsoid
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from cells import EllipsoidCell
import time
from typing import Tuple
from ellipsoid import Ellipsoid
import generate_emitters

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d
    import time


    fig = plt.figure()
    ax = plt.axes(projection="3d")
    time1 = time.time()
    emitter_positions, cell= generate_emitters.cell_emitter_position(
        (0, 1000),
        (0, 1000),
        (0, 1000),
        5000,
        6,
        {"a": [100, 120], "b": [150, 300], "c": [80, 100]},
        True,
    )
    time2 = time.time()
    print(time2 - time1)
    ax.scatter3D(
        emitter_positions[:, 0], emitter_positions[:, 1], emitter_positions[:, 2], s=3
    )
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, 1000)
    ax.set_zlim(0, 1000)
    plt.show()
