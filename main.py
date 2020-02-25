# %%
from numpy import genfromtxt
import matplotlib.pyplot as plt
from scipy.io import loadmat
import numpy as np
from astropy.io import fits
# %matplotlib
import scipy.stats as st

def create_fits(filename,pixelarray,nan=False):
    if nan:
        pixelarray[pixelarray==0] = np.nan
    fits.writeto(filename+".fits", pixelarray, overwrite=True)


def gauss_kernel(kernlen = 9,nsig = 4):
    interval = (2 * nsig + 1.) / (kernlen)
    x = np.linspace(-nsig - interval / 2., nsig + interval / 2., kernlen + 1)
    kern1d = np.diff(st.norm.cdf(x))
    # kernel_raw = np.sqrt(np.outer(kern1d, kern1d))
    # kernel = kernel_raw/kernel_raw.sum()
    return kern1d
#%%

features_template = genfromtxt("sample_spectra.csv", delimiter=",")

pixelarray = loadmat("./6k_spectra.mat")['P']
pixelarray = np.rot90(pixelarray, 3)
pixelarray = pixelarray.round(6)

wavelengths = features_template[:,0]
wavelengths = wavelengths/10000
wavelengths = wavelengths.round(6)

features = features_template[:,1]

spectra = dict( (wavelengths[i],features[i]) for i in  range(len(features)))

targets = pixelarray.nonzero()
row,col = targets[0], targets[1]

kernel = gauss_kernel()

for i in range(len(row)):
    pixelarray[row[i],col[i]] =  spectra.get(pixelarray[row[i],col[i]],0)
    pixelarray[row[i]-4:row[i]+5, col[i]] = pixelarray[row[i],col[i]]*kernel

create_fits("wavelength_fitted", pixelarray)
