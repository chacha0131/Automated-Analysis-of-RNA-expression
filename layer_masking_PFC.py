import os
import numpy as np
import pandas as pd
from PIL import Image

# Set up experimental variables
path = "/Users/emilycha/Desktop/imgs_bagot/redo/PFC/out/quant/" #Change path 
files = os.listdir(path)

for file in files:
    print(file)
    size_ref = Image.open('/Users/emilycha/Desktop/imgs_bagot/redo/PFC/imgs/preproc/rna_Sdk1_denoised/' + file + '_Sdk1_denoised.tif') #Change path
    mask = Image.open('/Users/emilycha/Desktop/imgs_bagot/redo/PFC/mask/mask_tif/' + file + '.tif').convert('L') #Change path
    mask = mask.resize((size_ref.size[0],size_ref.size[1])) #Resize to fit images in case masks are scaled
    mask = np.array(mask)
    cells = pd.read_csv(path + file + '/cells_quant.csv')
    ## Area for each cell
    cells['layer'] = 'none'
    for c in range(len(cells)):
        cell = cells.loc[c]
        area = mask[round(cell.y), round(cell.x)]
        if (area == 64) or (area == 69):      #Change area to fit the different grayscale values of your mask
            cells['layer'].loc[c] = 'L1'
        elif (area == 77) or (area == 82):
            cells['layer'].loc[c] = 'L23'
        elif (area == 89) or (area == 94):
            cells['layer'].loc[c] = 'L5'
        elif (area == 102) or (area == 107):
            cells['layer'].loc[c] = 'L6'
    cells.to_csv(path + file + '/cells_quant_layers.csv')