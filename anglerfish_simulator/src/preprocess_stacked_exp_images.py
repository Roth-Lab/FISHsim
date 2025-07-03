"""
This is a script to preprocess real images by FOV/Z -- to put them into individual rounds, for subdivide.
"""

import skimage.io as skio
import os
import argparse
import numpy as np


def preprocess_experimental_images(output_dir_name, stacked_preprocessed_img_path, channel_list, channel_1, channel_2, fov_num, z_num, total_bits):

    pp_boxcox_channel_1_path = os.path.join(os.getcwd(), 'results', output_dir_name, '1_pp_boxcox', f'{channel_list[0]}nm, Raw')
        
    pp_boxcox_channel_2_path = os.path.join(os.getcwd(), 'results', output_dir_name, '1_pp_boxcox', f'{channel_list[1]}nm, Raw')
            
    if not os.path.exists(pp_boxcox_channel_1_path):
        # Create the directory if it does not exist
        os.makedirs(pp_boxcox_channel_1_path)        

    if not os.path.exists(pp_boxcox_channel_2_path):
        # Create the directory if it does not exist
        os.makedirs(pp_boxcox_channel_2_path)
            
            
    stacked_image = skio.imread(stacked_preprocessed_img_path)
    
    for s in range(stacked_image.shape[0]):
        
        if s+1 in channel_1:
            
            round_num = np.argwhere(np.array(channel_1) == s+1)[0][0]
            
            skio.imsave(os.path.join(pp_boxcox_channel_1_path, f'merFISH_{round_num+1:02}_{fov_num:03}_{z_num:02}.tiff'), stacked_image[s,:,:][:,:,None])
            
        elif s+1 in channel_2:
            
            round_num = np.argwhere(np.array(channel_2) == s+1)[0][0]
            
            skio.imsave(os.path.join(pp_boxcox_channel_2_path, f'merFISH_{round_num+1:02}_{fov_num:03}_{z_num:02}.tiff'), stacked_image[s,:,:][:,:,None])
            
            


def main():
    
    parser = argparse.ArgumentParser(description='Splitting stacked preprocessed experimental images back into individual rounds, prior to AnglerFISH inference.')
    
    parser.add_argument('--output-dir-name', type=str, required=True, help='Output directory name')
    parser.add_argument('--stacked-preprocessed-img-path', type=str, required=True, help='Input image path')
    parser.add_argument('--channel-list', nargs='+', type=int, required=True, help='List of channels, example: 650 561')
    parser.add_argument('--channel-1', nargs='+', type=int, required=True, help='Channel 1 list, example: 1 4 5 7 10 11 13 16')
    parser.add_argument('--channel-2', nargs='+', type=int, required=True, help='Channel 2 list, example: 2 3 6 8 9 12 14 15')
    parser.add_argument('--fov-num', type=int, required=True, help='fov number of the image you wish to process')
    parser.add_argument('--z-num', type=int, required=True, help='z number of the image you wish to process')
    parser.add_argument('--total-bits', type=int, default=16, help='Total number of bits, example: 16')

    args = parser.parse_args()

    preprocess_experimental_images(args.output_dir_name, args.stacked_preprocessed_img_path, args.channel_list, args.channel_1, args.channel_2, args.fov_num, args.z_num, args.total_bits)




if __name__ == "__main__":
    
    main()
    
    '''     
    python src\preprocess_stacked_exp_images.py --fov-num 18 --z-num 6 --output-dir-name 'merfish_xp8054' --stacked-preprocessed-img-path 'C:/Users/jenkints/Documents/GitHub/anglerfish_simulator/anglerfish_simulator/results/merfish_xp8054/1_pp/xp8054_018_06_boxcox_cropped.tif' --channel-list 650 561 --channel-1 1 4 5 7 10 11 13 16 --channel-2 2 3 6 8 9 12 14 15 --total-bits 16
    '''