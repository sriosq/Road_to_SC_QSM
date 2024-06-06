import nibabel as nib
import os
import numpy as np
import matplotlib.pyplot as plt

from label import SegmentationLabel
from Volume import Volume


np.set_printoptions(precision=2, suppress = True)
 
ct_wb = nib.load("data/ct_wb.nii.gz")
ct_wb_seg = nib.load("data/totalseg.nii.gz")
# This variables are instances of the nibabel image
# To work with the data:
ct_wb_data = ct_wb.get_fdata()
labels = ct_wb_seg.get_fdata()
ct_wb_data.shape


full_body_phantom = Volume(labels)

full_body_phantom.manual_labeling()

full_body_phantom.set_susceptibility()





