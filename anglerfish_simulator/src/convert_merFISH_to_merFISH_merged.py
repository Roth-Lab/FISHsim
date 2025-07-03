# -*- coding: utf-8 -*-
"""
This script converts simulator-generated images that's in two-channel formats to merFISH_merged.
"""

import os
import numpy as np
import skimage.io as skio
import pandas as pd

def convert_to_merFISH_merged(simulated_images_path, data_organization_path, image_size, total_fov):
    
    data_org = pd.read_csv(data_organization_path)
    
    color_channels_list = np.sort(np.unique(data_org[['color']])) # get unique color channels, and also sort so lowest color channel first
    
    # Assuming we have two (in exception case, we have three) color channels
    channel_650_path = os.path.join(simulated_images_path, f'{color_channels_list[0]}nm, Raw')
    channel_750_path = os.path.join(simulated_images_path, f'{color_channels_list[1]}nm, Raw')
    
    fiducial_channel_list = np.unique(data_org[['fiducialColor']])
    
    fiducialchannel_561_path = os.path.join(simulated_images_path, f'{fiducial_channel_list[0]}nm, Raw')
    
    total_imagingRounds = np.unique(data_org[['imagingRound']]).shape[0]
    
    merFISH_merged_folder_path = os.path.join(simulated_images_path, 'merFISH_merged')
    
    os.makedirs(merFISH_merged_folder_path, exist_ok=True)
    
    for f in range(total_fov):
        # There's only 1 z-slice for now. Hence, no need to loop over z.
                
        for ir in range(total_imagingRounds):
            
            stacked_image = np.zeros((8, image_size, image_size), dtype=np.uint16)  # (Z, H, W) format for multi-page TIFF

            stacked_image[0] = skio.imread(os.path.join(channel_650_path, f"merFISH_{ir+1:02}_{f+1:03}_01.tiff")).squeeze()
            
            stacked_image[1] = skio.imread(os.path.join(channel_750_path, f"merFISH_{ir+1:02}_{f+1:03}_01.tiff")).squeeze()
            
            stacked_image[2] = skio.imread(os.path.join(fiducialchannel_561_path, f"merFISH_{ir+1:02}_{f+1:03}_01.tiff")).squeeze()
        
            # Save as a multi-page grayscale TIFF
            skio.imsave(os.path.join(merFISH_merged_folder_path, f"merFISH_merged_{ir+1:02}_{f+1:03}.tiff"), 
                        stacked_image, plugin="tifffile", photometric="minisblack")

if __name__ == "__main__":

    anglerfish_simulator_directory = r"C:\Users\jenkints\Documents\anglerfish_simulator\anglerfish_simulator"

    simulated_images_folder_name = 'merfish_serval_experiment_default'
    
    simulated_images_path = os.path.join(anglerfish_simulator_directory, "results", simulated_images_folder_name, '1')
    
    data_organization_path = os.path.join(anglerfish_simulator_directory, "resources", "data_organizations", "data_organization.csv")
    
    convert_to_merFISH_merged(simulated_images_path, data_organization_path)
            
        