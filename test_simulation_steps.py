import nibabel as nib
import os
import numpy as np
from simulation_functions import *

from label import SegmentationLabel
from Volume import Volume

np.set_printoptions(precision=2, suppress=True)


ct_wb_seg = nib.load("post_processing/final_total_seg.nii.gz")

# These variables are instances of the nibabel image
# To work with the data:
# ct_wb_data = ct_wb.get_fdata()
labels = ct_wb_seg.get_fdata()

## Step 1: Creating the dipole kernel
#print(ct_wb_seg.header['pixdim'][1:4])

B0_dir = [0, 0, 1]
voxel_size = ct_wb_seg.header['pixdim'][1:4]
dim = np.array(labels.shape)

D = create_dipole_kernel(B0_dir, voxel_size, dim)

#print("Dipole kernel shape:", D.shape)

# Step 1: complete

print(ct_wb_seg.header)



