import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.spatial.transform import Rotation as R

import ellipsoid
import cells

# [33.55934571 65.68648844 50.00579853]
# [ 87.60376451  57.50128419 520.75555253]
# [[ 0.56654561, -0.79989604, -0.1979707 ],
#  [-0.63912717, -0.57819577,  0.50715492],
#  [-0.52013703, -0.16079794, -0.83880957]]


# cell_axes_bounds = {"a": [80, 80], "b": [30, 30], "c": [20, 20]}
# cell = cells.EllipsoidCell(
#     [0, 0, 100], cell_axes_bounds, [50, 10, 60], nucleus_size=0.3
# )
cell_axes_bounds = {"a": [33, 33], "b": [65, 65], "c": [50, 50]}
cell = cells.EllipsoidCell(
    [0, 0, 100], cell_axes_bounds, [180 - 169.1, 31, 180 - 48.4], nucleus_size=0.3
)
# r = R.from_matrix(
#     [
#         [0.56654561, -0.79989604, -0.1979707],
#         [-0.63912717, -0.57819577, 0.50715492],
#         [-0.52013703, -0.16079794, -0.83880957],
#     ]
# )
# a = r.as_euler("xyz", degrees=True)
# print(f"angle: {a}")
# cell.shape.R = np.array(
#     [
#         [0.56654561, -0.79989604, -0.1979707],
#         [-0.63912717, -0.57819577, 0.50715492],
#         [-0.52013703, -0.16079794, -0.83880957],
#     ]
# )
# cell.R = cell.shape.R
cell.generate_emitters(1000)
major, minor, center = cell.shape.project()
print(major)
img = cell.background(20000 / 10)

# plot 3d cell structure
nrow, ncol = 1, 3
fig = plt.figure()
ax1 = fig.add_subplot(nrow, ncol, 1)
ax1.scatter(cell.emitters[:, 0], cell.emitters[:, 1])
ax1.set_aspect(aspect="equal")

# plot cell projection
ax2 = fig.add_subplot(nrow, ncol, 2, sharex=ax1, sharey=ax1)
ax2.plot([major[0], 0, minor[0]], [major[1], 0, minor[1]])
ax2.set_aspect(aspect="equal")
ax3 = fig.add_subplot(nrow, ncol, 3)
ax3.imshow(img)

plt.show()

