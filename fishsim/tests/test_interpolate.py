import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import interpolate
from utils import BASE_PROJECT_DIR
from simulation import Simulator

psf_path = BASE_PROJECT_DIR / "resources/PSF.mat"
mat = scipy.io.loadmat(str(psf_path))
psf = np.array(mat["ans"])
psf = Simulator.process_psf(psf)

print(psf.shape)
# plt.imshow(psf[:, :, 800])
# plt.show()

rows = list(range(psf.shape[0]))
cols = list(range(psf.shape[1]))
f = interpolate.RectBivariateSpline(rows, cols, psf[:, :, 800])
# f = interpolate.interp2d(rows, cols, psf[:, :, 800], kind="cubic") # Not efficient for data in a grid

rows_new = np.arange(0, psf.shape[0] - 1 + 0.5, 0.5)
cols_new = np.arange(0, psf.shape[1] - 1 + 0.5, 0.5)
z_new = f(rows_new, cols_new)
print(z_new.shape)
print(rows_new.shape)
print(psf[13, :, 800].shape)

fig, axs = plt.subplots(3)
axs[0].imshow(psf[:, :, 800])
axs[1].imshow(z_new)
axs[2].plot(rows, psf[15, :, 800], "ro-", rows_new, z_new[30, :], "b-")
plt.show()
