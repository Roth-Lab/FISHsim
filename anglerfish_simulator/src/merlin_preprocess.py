""" Pre-processing pipeline taken from Starfish
"""
from scipy.signal import convolve, fftconvolve
from skimage import img_as_float32

import numpy as np
import skimage.io


def preprocess(in_file, out_file, decon_sigma=2, high_pass_sigma=3, low_pass_sigma=1, num_decon_iters=15):
    imgs = skimage.io.imread(in_file).astype(np.uint16)

    imgs = img_as_float32(imgs)

    kernel_size = int(2 * np.ceil(2 * decon_sigma) + 1)
    
    psf = gaussian_kernel(
        shape=(kernel_size, kernel_size),
        sigma=decon_sigma
    )
    
    for i in range(imgs.shape[0]):
        imgs[i] = levels(
            high_pass_filter(imgs[i], high_pass_sigma)
        ) 

        imgs[i] = richardson_lucy_deconv(imgs[i], num_decon_iters, psf)
        
        imgs[i] = levels(low_pass_filter(imgs[i], low_pass_sigma))
    
    skimage.io.imsave(out_file, imgs, imagej=True)


def high_pass_filter(image, sigma):
    blurred = low_pass_filter(image, sigma)
    blurred = levels(blurred)  # clip negative values to 0.
    filtered = image - blurred
    return filtered


def low_pass_filter(image, sigma):
    return skimage.filters.gaussian(
        image,
        sigma=sigma,
        cval=0,
        preserve_range=True,
        truncate=4.0
    )


def richardson_lucy_deconv(image, iterations, psf):
    direct_time = np.prod(image.shape + psf.shape)

    fft_time = np.sum([n * np.log(n) for n in image.shape + psf.shape])

    time_ratio = 40.032 * fft_time / direct_time

    if time_ratio <= 1 or len(image.shape) > 2:
        convolve_method = fftconvolve

    else:
        convolve_method = convolve

    image = image.astype(float)

    psf = psf.astype(float)

    im_deconv = 0.5 * np.ones(image.shape)

    psf_mirror = psf[::-1,::-1]

    eps = np.finfo(image.dtype).eps

    for _ in range(iterations):
        x = convolve_method(im_deconv, psf, 'same')

        np.place(x, x == 0, eps)

        relative_blur = image / x + eps

        im_deconv *= convolve_method(relative_blur, psf_mirror, 'same')

    if np.all(np.isnan(im_deconv)):
        raise RuntimeError(
            'All-NaN output data detected. Likely cause is that deconvolution has been run for '
            'too many iterations.')

    return np.array(im_deconv)

    
def gaussian_kernel(shape=(3, 3), sigma=0.5):
    m, n = [int((ss - 1.) / 2.) for ss in shape]
    
    y, x = np.ogrid[-m:m + 1, -n:n + 1]
    
    kernel = np.exp(-(x * x + y * y) / (2. * sigma * sigma))
    
    kernel[kernel < np.finfo(kernel.dtype).eps * kernel.max()] = 0
    
    sumh = kernel.sum()
    
    if sumh != 0:
        kernel /= sumh
    
    return kernel


def levels(array, rescale=False, rescale_saturated=False):
    if rescale and rescale_saturated:
        raise ValueError("rescale and rescale_saturated cannot both be set.")

    data = array

    data = data.astype(np.float32, copy=False)

    # we don't want a copy, so we just do it in place.
    aboveone = np.any(data > 1)
    
    do_rescale = bool(rescale or (rescale_saturated and aboveone))
    
    _adjust_image_levels_in_place(data, do_rescale)

    array = data

    return array


def _adjust_image_levels_in_place(array, rescale=False):
    assert array.dtype == np.float32
    
    array[array < 0] = 0
    
    if rescale:
        array /= array.max()
    
    else:
        array[array > 1] = 1
