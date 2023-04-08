### Baysor Prep
# Prepare files for baysor.
#
# Author: Heike Schuler
# Date: 4/8/22

import os,sys
import numpy as np
import pandas as pd
from PIL import Image

# Set up experimental variables
path = "/Users/emilycha/Desktop/imgs_bagot/PVT"
channels = [channel for channel in os.listdir(os.path.join(path, 'imgs', 'preproc')) if channel.startswith('rna_')]
files = os.listdir(os.path.join(path, 'out', 'rna', channels[0][4:], 'spots'))

for file in files:
    id = os.path.splitext(os.path.basename(file))[0]
    id = id[0:7]
    print(id)
    # Load mask
    mask = Image.open(os.path.join(path, 'mask', 'mask_tif', id + '.tif')).convert('L')
    cells = Image.open(os.path.join(path, 'out', 'nuclei', id + '_nuclei_cp_masks.png')) # '_nuclei_ilcor_cp_masks.png'))
    mask = mask.resize((cells.size[0],cells.size[1])) #Resize to fit images in case masks are scaled (PFC images only(?))
    mask = np.array(mask)
    mask = mask > 0 # mask > 0 (white cutout: >0; white background: <255)
    mask_xy = np.c_[np.where(mask)[1], np.where(mask)[0]]
    mask_xy = pd.DataFrame(mask_xy, columns = ['x','y'])
    # Subset, mask, integrate and save spots files
    all_spots=pd.DataFrame(columns = ['x','y','gene'])
    for channel in channels:
        rna = channel[4:]
        spots = np.load(os.path.join(path, 'out', 'rna', rna, 'spots', id + '.npy'))
        spots = spots[:,(1,0)].astype(int)
        spots = pd.DataFrame(spots, columns = ['x','y'])
        spots_sub = pd.merge(spots, mask_xy, 'inner')
        spots_sub['gene'] = np.array([rna] * spots_sub.shape[0])
        all_spots = all_spots.append(spots_sub)
    all_spots.to_csv(os.path.join(path, 'out', 'rna', 'rna_clean', id + '.csv'), header = True, index = False)
    # Mask and save cell segmentation
    cells = np.array(cells).astype('uint16')
    cells[mask == False] = 0
    cells = Image.fromarray(cells)
    cells.save(os.path.join(path, 'out', 'nuclei_masked', id + '.tif'))

### END ###