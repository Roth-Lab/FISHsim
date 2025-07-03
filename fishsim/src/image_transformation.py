import numpy as np
import skimage.io as skio
import matplotlib.pyplot as plt
import seaborn as sns
import os
import cv2
import numba

def plot_X_distribution(X_transformed, transformation_type, output_path):
    # Assuming X_transformed is your transformed 2D array
    
    num_columns = X_transformed.shape[1]  # Number of columns
    
    # Setting up the plot grid
    fig, axs = plt.subplots(num_columns, 1, figsize=(10, 5*num_columns))
    
    # Looping through each column and plot
    for i in range(num_columns):
        column_data = X_transformed[:, i]
        
        # Using Seaborn to plot histogram and KDE
        sns.histplot(column_data, bins=30, kde=True, ax=axs[i])
        
        # Adding title and labels
        axs[i].set_title(f'Distribution of Column {i+1}')
        axs[i].set_xlabel('Value')
        axs[i].set_ylabel('Frequency')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'img_dist_{transformation_type}.png'))


def reshape_data(imgs):
    
    return imgs.reshape(
        (imgs.shape[0], imgs.shape[1] * imgs.shape[2])).T



def apply_ade(imgs):
    
    # Assuming img is float, and not yet in uint16.
    imgs = (imgs * 65535).astype(np.uint16)

    # Assuming img is your 3D image array with shape (16, 1608, 1608)
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    
    # Initialize an empty array to hold the processed images
    imgs_ahe = np.empty_like(imgs)
    
    # Apply AHE to each 2D image in the 3D array
    for i in range(imgs.shape[0]):
        imgs_ahe[i] = clahe.apply(imgs[i])
        
    # Convert uint16 to float32 and normalize to [0.0, 1.0] range
    img_float32 = imgs_ahe.astype(np.float32) / 65535.0

    return img_float32


def apply_gaussian_filter(imgs, sigma=1):

    # Assume 'stack' is your 3D numpy array with shape (16, height, width)
    smoothed_stack = np.empty_like(imgs, dtype=np.float32)  # assuming the dtype of the image is float32
    
    for i in range(imgs.shape[0]):
        smoothed_stack[i] = cv2.GaussianBlur(imgs[i], (5, 5), sigma)  # You can change kernel size and sigma as per your requirement.

    return smoothed_stack


def apply_bilateral_filter(imgs, d=5, sc=2, ss=50):
    
    # Assuming X is your stack of images of shape (16, height, width)
    filtered_images = np.empty_like(imgs)
    
    for i in range(imgs.shape[0]):
        # Apply bilateral filter
        # You can adjust the parameters: d, sigmaColor, and sigmaSpace as needed
        filtered_image = cv2.bilateralFilter(imgs[i], d=d, sigmaColor=sc, sigmaSpace=ss)
        filtered_images[i] = filtered_image

    return filtered_images



def gamma_correction_stack(images, gamma=1.1):
    """Apply gamma correction to a stack of images."""
    inv_gamma = 1.0 / gamma
    corrected_images = np.clip((images ** inv_gamma) * 255, 0, 255).astype(np.uint8)
    return corrected_images / 255.0  # Convert back to float32 between [0, 1]



def reshape_tX(transformed_imgs):
    
    img_axis_size = int(np.sqrt(transformed_imgs.shape[0]))
    ir = transformed_imgs.shape[1]
    
    tX = transformed_imgs.T.reshape((ir,img_axis_size,img_axis_size))
    
    return tX


@numba.njit
def normalizing_data(data):
    """ Normalize data to strictly (0,1) for KDE Mixture Distributions with data boundary requirements """
    # For simplicity and consistency reason, we will normal data to (0,1) for ALL types of Mixtures.
    # Despite (0,1) is only required for Beta Distribution.
    # Define a small epsilon to prevent 0 and 1 values in the normalized data
    epsilon = 1e-5

    # Min-max normalization to (0, 1)
    data_min = data.min()
    
    data_max = data.max()
    
    range = data_max - data_min
    
    normalized_data = (data - data_min) / range

    # Adjust normalization to (epsilon, 1-epsilon)
    normalized_data = (normalized_data * (1 - 2 * epsilon)) + epsilon

    return normalized_data





def compute_transformations(path, input_img_file_name, output_img_file_name, save_boxcox_image=True):
    
    imgs = skio.imread(os.path.join(path, input_img_file_name))

    
    ''' Reshape Data '''
    X = reshape_data(imgs)
    
    X = X.astype(np.float64)
    
    # X with rows = # of pixels and cols = # of imaging rounds
    print(X.shape)
    
    
    ''' Log Transformations '''
    X_log = np.log(X + 1e-3)
                
    X_log = (X_log - X_log.mean(axis=0)) / X_log.std(axis=0)
    
    
    ''' New transformations '''
    from sklearn.preprocessing import PowerTransformer
    
    # Initialize the transformer
    pt = PowerTransformer(method='box-cox', standardize=True)  # 'yeo-johnson' allows zero values, change to 'box-cox' if all values are positive.
    
    # Prepare the data
    X_boxcox = X + 1e-3  # Added 1e-3 to avoid zero
    
    # Fit and Transform
    X_boxcox = pt.fit_transform(X_boxcox)
    
    
    original = reshape_tX(X)
    
    log_transformed_img = reshape_tX(X_log)
    
    boxcox_transformed_img = reshape_tX(X_boxcox)
    
    
    original_normalized_image = np.array([normalizing_data(original[i])for i in range(original.shape[0])])

    log_normalized_image = np.array([normalizing_data(log_transformed_img[i])for i in range(log_transformed_img.shape[0])])
    
    boxcox_transformed_image = np.array([normalizing_data(boxcox_transformed_img[i])for i in range(boxcox_transformed_img.shape[0])])

    
    
    ''' Let's observe the distribution of each imaging round 
    1. Without transformation
    2. With current log + standardization
    3. With box-cox + standardization '''
    
    img_dist_output_path = os.path.join(path, 'img_distribution')
    
    os.makedirs(img_dist_output_path, exist_ok=True)

    plot_X_distribution(reshape_data(original_normalized_image), 'original', img_dist_output_path)
    
    plot_X_distribution(reshape_data(log_normalized_image), 'log', img_dist_output_path)
    
    plot_X_distribution(reshape_data(boxcox_transformed_image), 'boxcox', img_dist_output_path)
    
    # We only save the boxcox transformed image.
    if save_boxcox_image == True:
        
        skio.imsave(os.path.join(path, output_img_file_name), boxcox_transformed_image)


if __name__ == "__main__":

    path = 'C:/Users/jenkints/Documents/GitHub/anglerfish_simulator/anglerfish_simulator/results/merfish_xp8054_simulated_scr1.1_pc500_zdist18/1'
    
    input_img_file_name = 'stacked_image_preprocessed_001.tif'
    
    output_img_file_name = 'stacked_image_preprocessed_boxcox_001.tif'
    
    compute_transformations(path, input_img_file_name, output_img_file_name, save_boxcox_image=False)


