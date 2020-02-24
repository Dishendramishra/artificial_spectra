# %%
from numpy import genfromtxt
import matplotlib.pyplot as plt
from scipy.io import loadmat
import numpy as np
from astropy.io import fits
# %matplotlib

def create_fits(filename,pixelarray,nan=False):
    if nan:
        pixelarray[pixelarray==0] = np.nan
    fits.writeto(filename+".fits", pixelarray, overwrite=True)
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

for i in range(len(row)):
    pixelarray[row[i],col[i]] =  spectra.get(pixelarray[row[i],col[i]],0)

create_fits("wavelength_fitted", pixelarray)