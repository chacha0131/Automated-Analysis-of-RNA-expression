### Detect spots
# Wrapper for bigfish functions.
#
# Author: Heike Schuler
# Date: 25/7/22

# Prep environment
import os, sys
import numpy as np
import bigfish.detection as detection
import bigfish.stack as stack
import bigfish.plot as plot
import matplotlib.pyplot as plt
from torch import matrix_exp


def plot_spot_overlay(image, spots, path_output=None, show=False):
    
    # fig setup
    max_x = image.shape[1]-8
    max_y = image.shape[0]-8
    x = [row[1] for row in spots if 8 < row[0] < max_y and 8 < row[1] < max_x]
    y = [row[0] for row in spots if 8 < row[0] < max_y and 8 < row[1] < max_x]

#plot
    fig = plt.figure()
    px = 1/plt.rcParams['figure.dpi']
    fig.set_size_inches(image.shape[1]*px,image.shape[0]*px)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    fig = ax.imshow(image, cmap='Greys_r')
    plt.scatter(x, y, marker='o', facecolors='none', edgecolors='red', s=75)
    if path_output is not None:
        plt.savefig(path_output)
    if show:
        plt.show()
    else:
        plt.close()
        
# Set up experimental variables
path = "/Users/emilycha/Desktop/Automated_analysis/PFC/cohort1/" # "D:\\PAVish_pilot\\"
channels = [channel for channel in os.listdir(os.path.join(path, 'imgs', 'preproc')) if ((channel.startswith('rna_')) & (channel.endswith('denoised')))]
id_char = 13 #4 WF

# Set up parameters 
voxel_size_nm = 156 #150 WF
spot_radius_nm_list = [624,624,624] #400 WF
alpha_list = [0.5,0.5,0.5] #[0.5,0.5,0.5] WF #[ch1,ch2,ch3] in alphabetical order
beta_list = [1,1,1] #[2,2,1.5] WF
gamma_list = [15,15,15] #[10,10,10] WF
spread_list = [5,5,6] #[5,8,5] WF

channels = [channels[2]]

# Process images
for channel in channels:
    path_in = os.path.join(path, 'imgs', 'preproc', channel)
    files = os.listdir(path_in)
    print('Processing ' + channel[4:] + ' channel')
    if not os.path.isdir(os.path.join(path, 'out', 'rna', channel[4:])):
        os.makedirs(os.path.join(path, 'out', 'rna', channel[4:]))
        os.makedirs(os.path.join(path, 'out', 'rna', channel[4:], 'spots'))
        os.makedirs(os.path.join(path, 'out', 'rna', channel[4:], 'elbow'))
        os.makedirs(os.path.join(path, 'out', 'rna', channel[4:], 'reference'))
        os.makedirs(os.path.join(path, 'out', 'rna', channel[4:], 'overlay'))
    ch_id = channels.index(channel) 
    alpha = alpha_list[ch_id]
    beta = beta_list[ch_id]
    gamma = gamma_list[ch_id]
    spread = spread_list[ch_id]
    spot_radius_nm = spot_radius_nm_list[ch_id]
    log_kernel_size = detection.get_object_radius_pixel(voxel_size_nm = voxel_size_nm, object_radius_nm = spot_radius_nm, ndim=2)
    minimum_distance = (log_kernel_size[0]/2, log_kernel_size[1]/2)
    #files = [files[0]]
    for file in files:
        id = os.path.splitext(os.path.basename(file))[0]
        ### Preprocessing image. 
        print('Preprocessing image ' + id[:id_char])
        img = stack.read_image(os.path.join(path_in, file))
        img_scaled = stack.rescale(img)
        ### Spot detection and decomposition
        print('Detecting ' + channel[4:] + ' signal in image ' + id[:id_char])
        spots = detection.detect_spots(
            images = img_scaled, 
            voxel_size = voxel_size_nm, #WF
            spot_radius = spot_radius_nm, #WF
            log_kernel_size= log_kernel_size,
            minimum_distance= minimum_distance,
            spread = spread
        )
        spots_post_decomposition, _, reference_spot = detection.decompose_dense(
            image = img_scaled, 
            spots = spots, 
            voxel_size = voxel_size_nm, 
            spot_radius = spot_radius_nm, 
            alpha = alpha, beta = beta, gamma = gamma)
        print(str(spots.shape[0]) + ' ' + channel [4:] + ' spots detected - ' +
            'Decomposed into ' + str(spots_post_decomposition.shape[0]) + ' spots')
        ### Save output and quality control plots
        print('Saving ' + channel[4:] + ' channel results for ' + id[:id_char])
        path_data = os.path.join(path, 'out', 'rna', channel[4:], 'spots', id[:id_char])
        np.save(path_data + '.npy', spots_post_decomposition)
        path_reference = os.path.join(path, 'out', 'rna', channel[4:], 'reference', id[:id_char])
        plot.plot_reference_spot(reference_spot, rescale=True, title=id[:id_char], 
            path_output= path_reference, show=False)
        path_elbow = os.path.join(path, 'out', 'rna', channel[4:], 'elbow', id[:id_char])
        plot.plot_elbow(img_scaled, spread, voxel_size = voxel_size_nm, spot_radius = spot_radius_nm, #WF
            title=id[:id_char], path_output=path_elbow, show=False)
        path_overlay = os.path.join(path, 'out', 'rna', channel[4:], 'overlay', id[:id_char])
        plot_spot_overlay(img_scaled, spots_post_decomposition, path_overlay + '.tif')
        del img, img_scaled, spots, spots_post_decomposition, reference_spot


### END ###
        
