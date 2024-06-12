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
        self.T2star_val = None
        self.PD_val = 0

        # Key is the name and value is ordered: M0, T1, T2, T2*, PD
        self.relax_values = {
            "air": [0, 0, 0, 0, 0],
            "bone": [None, 1204, 53, 3303, 117],  # M0 is often not specified for bone
            "lungs": [None, 1270, None, 0, 0],  # Air in lungs doesn't have M0, T2 values?
            "water": [None, 2500, 2500, 1, 100],  # High M0 value
            "CSF": [None, 3200, 2000, 1, 100],  # High M0
            "sc_csf": [None, 3200, 2000, 1, 100],
            "fat": [None, 380, 108, 0.035, 140], # T2star value : 0.5*70e-3 # Daniel PD=90
            "liver": [None, 809, 34, 34/2, 70],
            "spleen": [None, 1328, 61, 65/2, 80],
            # In this initial segmentation the whole brain will be considered 60% GM and 40% WM
            # Given the values a ponderated estimation is 60.8 ms
            "brain":[None,None,None, 60.8, 90],
            "white_matter": [None, None, None, 53.5/2,0],
            "gray_matter": [None, None, None, 66, 0],
            "sc_wm":[None, None, None, None, 0],
            "sc_gm":[None, None, None, None, 0],
            "heart":[1000 ,1300, 55, 18.5/2, 85],
            "kidney":[None, 1190, 56, 65.4/2, 70],
            "pancreas":[None, 725,43, 37, 75],
            "cartilage":[None, 1240,32, 20, 50], # PD value is a guess
            "bone_marrow":[None, 365, 23, None, 60],  # PD value is a guess
            "SpinalCanal":[None, 993, 78, 60, 100], #
            "esophagus":[None,None, None, 17, 35], #
            "trachea":[None, None, None, 25, 15],
            "organ":[None, 800, 34, 17, 50], # Values similar to those from liver
            "gland":[None, None, None, 50, 100],
            # There are some organs that don't have enough documentation on the literature to complete
            # the required values so an estimation is used for these:
            "extra" : [None, 750, 50, 0.035,120]
        }
    # Literature values from:
    # Jorge Zavala Bojorquez, Stéphanie Bricq, Clement Acquitter, François Brunotte, Paul M. Walker, Alain Lalande, What are normal relaxation times of tissues at 3 T?, Magnetic Resonance Imaging, Volume 35, 2017, Pages 69-80, ISSN 0730-725X, https://doi.org/10.1016/j.mri.2016.08.021.
    #
    # Stanisz, G.J., Odrobina, E.E., Pun, J., Escaravage, M., Graham, S.J., Bronskill, M.J. and Henkelman, R.M. (2005), T1, T2 relaxation and magnetization transfer in tissue at 3T. Magn. Reson. Med., 54: 507-512. https://doi.org/10.1002/mrm.20605
    # Arnold, J., Fidler, F., Wang, T. et al. Imaging lung function using rapid dynamic acquisition of T 1-maps during oxygen enhancement. Magn Reson Mater Phy 16, 246–253 (2004). https://doi.org/10.1007/s10334-004-0034-z
    # Meloni, A., De Marchi, D., Positano, V. et al. Accurate estimate of pancreatic T2* values: how to deal with fat infiltration. Abdom Imaging 40, 3129–3136 (2015). https://doi.org/10.1007/s00261-015-0522-9
    # Hesper, T., Hosalkar, H.S., Bittersohl, D. et al. T2* mapping for articular cartilage assessment: principles, current applications, and future prospects. Skeletal Radiol 43, 1429–1445 (2014). https://doi.org/10.1007/s00256-014-1852-3
    #
    # MRI from Picture to Proton
    # Questions and answers in MRI website
    ########### For some T2star values #############
    # Some T2 star values from the literature are at 1.5T: liver, spleen, kidney, WM and GM, cartilage
    # T2 star values should decrease with higher field strength
    # For the purpose of this code, the values at 1.5T are assumed to half at 3T

    # Cristina Rossi, Andreas Boss, Michael Haap, Petros Martirosian, Claus D. Claussen, Fritz Schick, Whole-body T2⁎ mapping at 1.5 T, Magnetic Resonance Imaging, Volume 27, Issue 4, 2009, Pages 489-496, ISSN 0730-725X, https://doi.org/10.1016/j.mri.2008.08.004.
    # For brain T2star values: Andrew M. Peters, Matthew J. Brookes, Frank G. Hoogenraad, Penny A. Gowland, Susan T. Francis, Peter G. Morris, Richard Bowtell, T2* measurements in human brain at 1.5, 3 and 7 T,
    # Magnetic Resonance Imaging, Volume 25, Issue 6, 2007, Pages 748-753, ISSN 0730-725X, https://doi.org/10.1016/j.mri.2007.02.014.

    # For Proton density:
    # Proton density should be independent of field strenght, we are using a value relative to water being 100
    #


    def set_name(self, name):

        if name in self.relax_values.keys():
            self.name = name
            self.M0 = self.relax_values[name][0]
            self.T1_val = self.relax_values[name][1]
            self.T2_val = self.relax_values[name][2]
            self.T2star_val = self.relax_values[name][3]
            self.PD_val = self.relax_values[name][4]

        else:

            self.name = name
            self.M0 = 0
            self.T1_val = 0
            self.T2_val = 0
            self.T2star = 0
            self.PD_val = 0


    def set_susceptibility(self, susceptibility):
        self.susceptibility = susceptibility

    def set_M0_val(self,M0):
        self.M0_val = M0
    def set_t1_val(self,t1):
        self.T1_val = t1

    def set_t2_val(self,t2):
        self.T2_val = t2
    def set_pd_val(self,pd):
        self.PD_val = pd

    def set_t2star_val(self,t2star):
        self.T2star_val = t2star

    def __repr__(self):
        # Add the latest attributes additioned to the class
        return (f"SegmentationLabel(label_id={self.label_id}, name={self.name}, susceptibility={self.susceptibility},"
                f"M0={self.M0_val}, T1,T2,T2* = {self.T1_val,self.T2_val,self.T2star_val}, PD = {self.PD_val} )")