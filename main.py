# %%
from numpy import genfromtxt
import matplotlib.pyplot as plt
from scipy.io import loadmat
import numpy as np
from astropy.io import fits
# %matplotlib

features_template = genfromtxt("sample_spectra.csv", delimiter=" ")
wavelength_x, values_y = features_template[:, 0], features_template[:, 2]

pixelarray = loadmat("./6k_spectra.mat")['P']
pixelarray = np.rot90(pixelarray, 3)

def create_fits(filename,pixelarray,nan=False):
    if nan:
        pixelarray[pixelarray==0] = np.nan
    fits.writeto("./fits/"+filename+".fits", pixelarray, overwrite=True)


## returns starting and ending point of all order
def get_order_locations():
    starting_points = np.nonzero(pixelarray[:, 0])[0]
    ending_points = np.nonzero(pixelarray[:, 6143])[0]
    orders_locations = np.column_stack((starting_points, ending_points))
    return orders_locations


# returns 2 array of wavelength-ranges for orders
def get_wavelength_ranges():
    ol = get_order_locations()
    wavelength_ranges = []

    for i, j in ol:
        wavelength_ranges.append([pixelarray[i, 0], pixelarray[j, 6143]])
    return wavelength_ranges


# %%
from skimage.draw import line_aa
ol = get_order_locations()
ccd_center = (6144//2)-1 # center line of ccd
lower = 4608-1
upper = 1536
ff_upper = 5376 
ff_lower = 768

for i,j in ol:
    # i,j = 5633, 5730 
    width = j-i+1
    sec = np.copy(pixelarray[i:i + width, ])
    sec_center = sec[:, ccd_center].nonzero()[0]
    
    order_center = sec_center[len(sec_center) // 2]
    
    lower_center = sec[:, lower].nonzero()[0]
    lower_center = lower_center[len(lower_center) // 2]
    
    upper_center = sec[:, upper].nonzero()[0][0]
    
    ff_lower_center = sec[:,ff_lower].nonzero()[0][0]
    ff_upper_center = sec[:,ff_upper].nonzero()[0][0]
    
    fill_value = 0
    
    sec[ff_lower_center:,:ff_lower] = fill_value
    sec[:ff_upper_center,ff_upper:] = fill_value
    
    sec[order_center + 1:, :ccd_center] = fill_value
    sec[:order_center, ccd_center:] = fill_value
    sec[lower_center + 1:, ccd_center:lower] = fill_value 
    sec[:upper_center, upper:ccd_center] = fill_value
     
    sec[sec == 0] = np.nan
    
    # plt.imshow(sec)
    # print(i,sec_center, order_center)
    create_fits("sector"+str(i)+"_"+str(j), sec)
    # break
    # if i==3742:
    #     break
#%%
# Now Trace Orders and replace with wavelengths
