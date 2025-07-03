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
    # cell axis nound
    cell_axes_bounds_1 = {"a": [20, 20], "b": [20, 20], "c": [20, 20]} 
    cell_axes_bounds_2 = {"a": [20, 20], "b": [60, 60], "c": [30, 30]}

    # Create Ellipsoif for checking for overlap
    c1 = Ellipsoid([40,0,0], [20,20,20], [0,0,0])
    c2 = Ellipsoid([0,0,0], [20,60,30], [0,0,90])
    if (Ellipsoid.check_overlap(c1, c2)):
        print("it is overlapping")
    print("it is not overlapping")

    # Reconstruct the EllipsoidCell to visualized
    c11 = EllipsoidCell([40,0,0], cell_axes_bounds_1, [0,0,0])
    c21 = EllipsoidCell([0,0,0], cell_axes_bounds_2, [0,0,90])
    c11.generate_emitters(1000)
    c21.generate_emitters(1000)
    
  
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.set_box_aspect([1,1,1])
    ax.set_xlim3d(-100, 100)
    ax.set_ylim3d(-100, 100)
    ax.set_zlim3d(-100, 100)
    ax.scatter3D(np.append(c11.emitters[:,0], c21.emitters[:,0]), np.append(c11.emitters[:,1], c21.emitters[:,1]), np.append(c11.emitters[:,2], c21.emitters[:,2]))
    plt.show()
   