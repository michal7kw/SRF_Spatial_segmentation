import tifffile
import numpy as np
from pathlib import Path

base_path = Path('.')

# Create max projections for each channel
channels = ['DAPI', 'Pcp4', 'PolyT', 'Prox1', 'Snap25', 'TAG', 'Ttr', 'hSETBP1']
for channel in channels:
    z_stack = []
    for z in range(7):
        img = tifffile.imread(f'{base_path}/mosaic_{channel}_z{z}.tif')
        z_stack.append(img)
    
    max_projection = np.max(z_stack, axis=0)
    tifffile.imwrite(f'{base_path}/mosaic_{channel}_maxproj.tif', max_projection)