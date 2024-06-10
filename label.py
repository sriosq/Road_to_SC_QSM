#Dependencies
import numpy as np


class SegmentationLabel:
    def __init__(self, label_id, name=None, susceptibility=None, ct_number = None):
        self.label_id = label_id
        self.name = name
        self.susceptibility = susceptibility
        self.ct_number = ct_number
        self.M0_val = None
        self.T1_val = None
        self.T2_val = None
        # Key is the name and value is ordered: M0, T1, T2
        self.relax_values = {
            "bone": [None, 1204, 53],  # M0 is often not specified for bone
            "air": [None, None, None],  # Air in lungs doesn't have M0, T1, T2 values?
            "water": [None, 2500, 2500],  # High M0 value
            "CSF": [None, 3200, 2000],  # High M0
            "fat": [None, 380, 108],
            "liver": [None, 809, 34],
            "spleen": [None, 1328, 61],
            "brain_wm" : [None, ],
            "brain_gm" : [None, None, None],
            "sc_wm":[None, None, None],
            "sc_gm":[None, None, None],
            "heart":[1000,1300,55],
            "kidney":[None,1190,56],
            "pancreas":[None,725,43],
            "cartilage":[None,1240,32,],
            "bone_marrow":[None,365,23],
            "spinal_cord":[None,993,78]

        }
    # Literature values from:
    # Jorge Zavala Bojorquez, Stéphanie Bricq, Clement Acquitter, François Brunotte, Paul M. Walker, Alain Lalande, What are normal relaxation times of tissues at 3 T?, Magnetic Resonance Imaging, Volume 35, 2017, Pages 69-80, ISSN 0730-725X, https://doi.org/10.1016/j.mri.2016.08.021.
    # MRI from Picture to Proton
    # Questions and answers in MRI website
    # Stanisz, G.J., Odrobina, E.E., Pun, J., Escaravage, M., Graham, S.J., Bronskill, M.J. and Henkelman, R.M. (2005), T1, T2 relaxation and magnetization transfer in tissue at 3T. Magn. Reson. Med., 54: 507-512. https://doi.org/10.1002/mrm.20605


    def set_name(self, name):
        self.name = name
    
    def set_susceptibility(self, susceptibility):
        self.susceptibility = susceptibility
    
    def __repr__(self):
        return f"SegmentationLabel(label_id={self.label_id}, name={self.name}, susceptibility={self.susceptibility})"