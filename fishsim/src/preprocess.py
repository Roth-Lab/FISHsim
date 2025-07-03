
import os
import numpy as np
import skimage.io as skio
import os
import argparse
import numba

from sklearn.preprocessing import PowerTransformer
from distutils.dir_util import copy_tree

from merlin_preprocess import preprocess
import image_transformation 


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



def savannah_preprocess(output_dir_name, channel_list, channel_1, channel_2, total_bits, total_fov, image_size):
    
    output_dir_path = os.path.join(os.getcwd(), 'results', output_dir_name, '1')
            
    # Check if the directory exists        
    pp_channel_1_path = os.path.join(os.getcwd(), 'results', output_dir_name, '1_pp', f'{channel_list[0]}nm, Raw')
    
    pp_channel_2_path = os.path.join(os.getcwd(), 'results', output_dir_name, '1_pp', f'{channel_list[1]}nm, Raw')
        
    if not os.path.exists(pp_channel_1_path):
        # Create the directory if it does not exist
        os.makedirs(pp_channel_1_path)        

    if not os.path.exists(pp_channel_2_path):
        # Create the directory if it does not exist
        os.makedirs(pp_channel_2_path)
            
    
    for f in range(total_fov):
        
        stacked_image = np.zeros((total_bits,image_size,image_size), dtype=np.int16)
        
        for b in range(total_bits):
            
            if b+1 in channel_1:
        
                channel_num = channel_list[0]
                
                round_num = int(np.ceil((b+1)/2))
                
            else:
                
                channel_num = channel_list[1]
                
                round_num = int(np.ceil((b+1)/2))
        
                
             
            print(channel_num, round_num)
            
            unstacked_img_path = os.path.join(output_dir_path, f'{channel_num}nm, Raw', f'merFISH_{round_num:02}_{f+1:03}_01.tiff')
            
            image = skio.imread(unstacked_img_path)
        
            stacked_image[b, :, :] = image.squeeze()
        
    
        stacked_img_path = os.path.join(output_dir_path, f'stacked_image_{f+1:03}.tif')
        
        skio.imsave(stacked_img_path, stacked_image)
        
        stacked_preprocessed_img_path = os.path.join(output_dir_path, f'stacked_image_preprocessed_{f+1:03}.tif')
        
        preprocess(stacked_img_path, stacked_preprocessed_img_path)


        ''' Now that we have preprocessed them, we will split the preprocessed stacked images and put them back into individual rounds.'''
        stacked_image = skio.imread(stacked_preprocessed_img_path)
        
        for s in range(stacked_image.shape[0]):
            
            if s+1 in channel_1:
                
                round_num = np.argwhere(np.array(channel_1) == s+1)[0][0]
                
                skio.imsave(os.path.join(pp_channel_1_path, f'merFISH_{round_num+1:02}_{f+1:03}_01.tiff'), stacked_image[s,:,:][:,:,None])
                
            elif s+1 in channel_2:
                
                round_num = np.argwhere(np.array(channel_2) == s+1)[0][0]
                
                skio.imsave(os.path.join(pp_channel_2_path, f'merFISH_{round_num+1:02}_{f+1:03}_01.tiff'), stacked_image[s,:,:][:,:,None])
    
    
        ''' We will compute the transformations (original, log and boxcox); and save the intensity distribution plots;
        and also save the boxcox transformed image. (added Dec 2023)
        
        This one takes way too long and visualization of plots is unncessary. Should be separate from this script'''
        #compute_transformations(output_dir_path, f'stacked_image_preprocessed_{f+1:03}.tif', f'stacked_image_preprocessed_boxcox_{f+1:03}.tif')
        
        
        """ BOXCOX Transformation and saving them into a separate folder named '1_pp_boxcox' """
        # Check if the directory exists        
        pp_boxcox_channel_1_path = os.path.join(os.getcwd(), 'results', output_dir_name, '1_pp_boxcox', f'{channel_list[0]}nm, Raw')
        
        pp_boxcox_channel_2_path = os.path.join(os.getcwd(), 'results', output_dir_name, '1_pp_boxcox', f'{channel_list[1]}nm, Raw')
            
        if not os.path.exists(pp_boxcox_channel_1_path):
            # Create the directory if it does not exist
            os.makedirs(pp_boxcox_channel_1_path)        

        if not os.path.exists(pp_boxcox_channel_2_path):
            # Create the directory if it does not exist
            os.makedirs(pp_boxcox_channel_2_path)
        
        stacked_preprocessed_boxcox_normalized_img_path = os.path.join(output_dir_path, f'stacked_image_preprocessed_boxcox_normalized_{f+1:03}.tif')

        
        ''' Read in preprocessed image '''
        imgs = skio.imread(stacked_preprocessed_img_path)

        
        ''' Reshape Data '''
        X = image_transformation.reshape_data(imgs)
        
        X = X.astype(np.float64)
        
        ''' New transformations '''
        print('Performing Boxcox Transformation on preprocessed images...')
        
        # Initialize the transformer
        pt = PowerTransformer(method='box-cox', standardize=True)  # 'yeo-johnson' allows zero values, change to 'box-cox' if all values are positive.
        
        # Prepare the data
        X_boxcox = X + 1e-3  # Added 1e-3 to avoid zero
        
        # Fit and Transform
        X_boxcox = pt.fit_transform(X_boxcox)


        ''' transform back to (ir, y_size, x_size) '''
        stacked_preprocessed_boxcox_image = image_transformation.reshape_tX(X_boxcox)
        
        ''' Normalizing image to (0,1) '''
        stacked_pp_boxcox_normalized_image = np.array([normalizing_data(stacked_preprocessed_boxcox_image[i])for i in range(stacked_preprocessed_boxcox_image.shape[0])])
        
        skio.imsave(stacked_preprocessed_boxcox_normalized_img_path, stacked_pp_boxcox_normalized_image)
        
        ''' Now that we have boxcox transformed them, we will split the images and put them back into individual rounds.'''
        for s in range(stacked_pp_boxcox_normalized_image.shape[0]):
            
            if s+1 in channel_1:
                
                round_num = np.argwhere(np.array(channel_1) == s+1)[0][0]
                
                skio.imsave(os.path.join(pp_boxcox_channel_1_path, f'merFISH_{round_num+1:02}_{f+1:03}_01.tiff'), stacked_pp_boxcox_normalized_image[s,:,:][:,:,None])
                
            elif s+1 in channel_2:
                
                round_num = np.argwhere(np.array(channel_2) == s+1)[0][0]
                
                skio.imsave(os.path.join(pp_boxcox_channel_2_path, f'merFISH_{round_num+1:02}_{f+1:03}_01.tiff'), stacked_pp_boxcox_normalized_image[s,:,:][:,:,None])



        ''' Copy 100nm, Raw; 473nm, Raw and groundtruths folders to 1_pp and 1_pp_boxcox'''
        
        original_groundtruths_directory = os.path.join(output_dir_path, 'groundtruths')
        
        original_100nm_directory = os.path.join(output_dir_path, '100nm, Raw')

        original_473nm_directory = os.path.join(output_dir_path, '473nm, Raw')

        
        destination_directory = os.path.join('results', output_dir_name)
        
        copy_tree(original_groundtruths_directory, os.path.join(destination_directory, '1_pp', 'groundtruths'))
        
        copy_tree(original_100nm_directory, os.path.join(destination_directory, '1_pp', '100nm, Raw'))

        copy_tree(original_473nm_directory, os.path.join(destination_directory, '1_pp', '473nm, Raw'))


        copy_tree(original_groundtruths_directory, os.path.join(destination_directory, '1_pp_boxcox', 'groundtruths'))
        
        copy_tree(original_100nm_directory, os.path.join(destination_directory, '1_pp_boxcox', '100nm, Raw'))

        copy_tree(original_473nm_directory, os.path.join(destination_directory, '1_pp_boxcox', '473nm, Raw'))




def main():
    
    parser = argparse.ArgumentParser(description='Process parameters for Savannah Preprocessing.')
    
    parser.add_argument('--output-dir-name', type=str, required=True, help='Name of the output directory')
    parser.add_argument('--channel-list', nargs='+', type=int, required=True, help='List of channels, example: 650 561')
    parser.add_argument('--channel-1', nargs='+', type=int, required=True, help='Channel 1 list, example: 1 4 5 7 10 11 13 16')
    parser.add_argument('--channel-2', nargs='+', type=int, required=True, help='Channel 2 list, example: 2 3 6 8 9 12 14 15')
    parser.add_argument('--total-bits', type=int, default=16, help='Total number of bits, example: 16')
    parser.add_argument('--total-fov', type=int, default=10, help='Total Field of View, example: 10')
    parser.add_argument('--image-size', type=int, default=1600, help='Image size, example: 1600')

    args = parser.parse_args()

    savannah_preprocess(args.output_dir_name, args.channel_list, args.channel_1, args.channel_2, args.total_bits, args.total_fov, args.image_size)




if __name__ == "__main__":
    
    main()
    
    '''
    channel_list=[650,561]
    
    channel_1=[1,4,5,7,10,11,13,16]
    
    channel_2=[2,3,6,8,9,12,14,15]
    
    total_bits=16
    
    total_fov=10
    
    image_size=1600
    
    output_dir_name = 'merfish_xp8054_simulated_scr2_pc1500_zdist14'
    
    savannah_preprocess(output_dir_name, channel_list, channel_1, channel_2, total_bits, total_fov, image_size)
    '''