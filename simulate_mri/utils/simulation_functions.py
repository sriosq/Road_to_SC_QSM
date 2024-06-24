# This python file contains functions needed to simulte GRE acquistion

# Creating a dipole kernel
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib

def create_dipole_kernel(B0_dir, voxel_size,dimensions,buffer=1):

    #B0_dir = np.array(B0_dir)
    #B0_dir = B0_dir / np.linalg.norm(B0_dir)
    dims = dimensions*buffer
    # Extract dimensions for each axis
    kmax = 1/(2*voxel_size)
    interval = 2*kmax/dims
    [kx, ky, kz] = np.meshgrid(np.arange(-kmax[0], kmax[0], interval[0]),
                                np.arange(-kmax[1], kmax[1], interval[1]),
                                np.arange(-kmax[2], kmax[2], interval[2]))
    k_sqrd = kx ** 2 + ky ** 2 + kz ** 2

    with np.errstate(divide='ignore', invalid='ignore'):
        kernel = np.fft.fftshift(1/3 - kz**2/k_sqrd)
        kernel[0,0,0] = 1/3

    return kernel

def show_slices(slices):
   """ Function to display row of image slices """
   fig, axes = plt.subplots(1, len(slices))
   for i, slice in enumerate(slices):
       axes[i].imshow(slice.T, cmap="gray", origin="lower")


def generate_signal(pd, T2star, FA, TE, deltaB0, gamma, handedness):
    # This code is for creating the volume dimensions with new dimension for each TE
    # FA is the flip angle in degrees
    # T2 star value will input in ms
    # TE will be a list of echo times in s
    # B0 is to specify the value of Tesla
    # Gamma in rad*Hz/Tesla

    if handedness == 'left':
        sign = -1
    elif handedness == 'right':
        sign = 1
    else:
        raise ValueError("Invalid handedness value")
    # Some regions of the volume have 0 as T2 star values so:
    #if T2star == 0:
    #    signal = pd * np.sin(np.deg2rad(FA))
    #else:

    # The FA input should be in degrees so we convert to radians
    # The lookup table has values in ms so we need to multiply by 10^-3
    # Because input TE is in seconds
    signal = pd * np.sin(np.deg2rad(FA)) * np.exp(-TE / T2star*1e-3- sign * 1j * gamma * deltaB0 * TE)

    mag = np.abs(signal)
    phase = np.angle(signal)

    return mag,phase

def optimized_signal(pd_vol,T2star_vol, FA, TE, deltaB0_vol, gamma, handedness):
    # This is an optimized version from generate_signal, using numpy array matrices

    decay = np.exp(-TE / (T2star_vol * 1e-3))  # Convert T2* to ms and apply decay
    phase_factor = -1j * gamma * deltaB0_vol * TE * 1e-3 if handedness == 'left' else 1j * gamma * deltaB0_vol * TE * 1e-3
    # Phase factor in radians

    signal = pd_vol * np.sin(FA) * decay * np.exp(phase_factor)

    return np.abs(signal), np.angle(signal) # Abs for the Magnitude whereas angle for Phase