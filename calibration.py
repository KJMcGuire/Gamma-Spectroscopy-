import data_arrays
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy.stats import linregress


N = 1011
x = np.arange(N)

###Subtract global background from known sample arrays
Cesium = data_arrays.sp1-data_arrays.sp0
Radium = data_arrays.sp2-data_arrays.sp0
Unknown = data_arrays.sp4-data_arrays.sp0

###Create list of known photopeak energies
energies = [661.69, 609.312, 351.932, 295.224, 241.997, 187.1]

###Create data slices for curve fitting
cesium_s1 = slice(0, 280, 1)  #Slice skips the first item, so this slice will take data lines 1-280
x1 = np.linspace(1, 280, 280)
cesium_s2 = slice(280, 348, 1)
x2 = np.linspace(281, 348, 68)
cesium_s3 = slice(347, 1011, 1)
x3 = np.linspace(348, 1011, 664)
radium_s1 = slice(254, 306, 1)
ra_x1 = np.linspace(255, 306, 52)
radium_s2 = slice(128, 168, 1)
ra_x2 = np.linspace(129, 168, 40)
radium_s3 = slice(99, 134, 1)
ra_x3 = np.linspace(100, 134, 35)
radium_s4 = slice(79, 112, 1)
ra_x4 = np.linspace(80, 112, 33)
radium_s5 = slice(63, 85, 1)
ra_x5 = np.linspace(64, 85, 22)

###Define a gaussian curve plus exponential term
def gaussian(x, sigma, X, A):
    return A*np.exp( (-(x-X)**2)/(2*sigma**2))

###Define an exponential to fit the background
def exponential(x, A, b):
    return A*np.exp(-x*b)

#Subtract local background from known spectra
Radium_clean = Radium - exponential(x, 3500, .0145)
Cesium_clean = Cesium - exponential(x, 3500, .0145)


ce_params1, ce_param_covar1 = curve_fit(gaussian, x2, Cesium[cesium_s2], p0 = [100, 340.,1.])
ra_params1c, ra_param_covar1c = curve_fit(gaussian, ra_x1, Radium_clean[radium_s1], p0 = [100, 340.,1.])
ra_params2c, ra_param_covar2c = curve_fit(gaussian, ra_x2, Radium_clean[radium_s2], p0 = [130., 220.,1.])
ra_params3c, ra_param_covar3c = curve_fit(gaussian, ra_x3, Radium_clean[radium_s3], p0 = [110., 140.,1.])
ra_params4c, ra_param_covar4c = curve_fit(gaussian, ra_x4, Radium_clean[radium_s4], p0 = [90., 115.,1.])
ra_params5c, ra_param_covar5c = curve_fit(gaussian, ra_x5, Radium_clean[radium_s5], p0 = [65., 95.,1.])

plt.figure(figsize=(10,7))
plt.plot(x2, gaussian(x2, *ce_params1), color = '#02ab2e',  linestyle = '-.', label='Fit curves')
plt.plot(ra_x1, gaussian(ra_x1, *ra_params1c), color = '#02ab2e', linestyle = '-.')
plt.plot(ra_x2, gaussian(ra_x2, *ra_params2c), color = '#02ab2e', linestyle = '-.')
plt.plot(ra_x3, gaussian(ra_x3, *ra_params3c), color = '#02ab2e', linestyle = '-.')
plt.plot(ra_x4, gaussian(ra_x4, *ra_params4c), color = '#02ab2e', linestyle = '-.')
plt.plot(ra_x5, gaussian(ra_x5, *ra_params5c), color = '#02ab2e', linestyle = '-.')\

#plt.plot(x, exponential(x, 3500, 0.0145))

plt.annotate('662 keV', xy=(320, 820), xytext=(320, 820))
plt.annotate('609 keV', xy=(260, 820), xytext=(260, 820))
plt.annotate('352 keV', xy=(150, 1600), xytext=(150, 1600))
plt.annotate('295 keV', xy=(110, 1280), xytext=(110, 1280))
plt.annotate('242 keV', xy=(85, 600), xytext=(85, 600))
plt.annotate('197 keV', xy=(52, 530), xytext=(52, 530))
plt.plot(x, Radium_clean, linewidth=2, color = '#2c6fbb', linestyle = '--', label='Ra-226')
plt.plot(x, Cesium_clean, linewidth=1, color = '#ff474c', label='Cs-137')
plt.xlim([0,400])
plt.ylim([0, 3000])
plt.xlabel("Channel No.")
plt.ylabel("Count")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("Gaussian_fit_to_photo_peaks.png", bbox_inches="tight", pad_inches=0)
plt.show()

#print ra_params1c
#print ra_params2c
#print ra_params3c
#print ra_params4c
#print ra_params5c

###Create lists of sigmas and channel numbers: the List Comprehension Way
energy_calibration1 = [t for t in zip(ce_params1, ra_params1c, ra_params2c, ra_params3c, ra_params4c, ra_params5c)[:2]]

#print energy_calibration1

###Perform curve fitting of energy versus channel no.
m, b, r_value, p_value, std_err = linregress(energy_calibration1[1], energies)

#Plot linear fit
x = np.arange(50, 350, 1)
plt.figure(figsize=(10,7))
fig, ax = plt.subplots()
ax.plot((x), m*(x) + b)
ax.scatter(energy_calibration1[1], energies, )
textstr = '\n'.join((
    r'$m=%.2f  (\sigma_{m}=0.02)$' % (m, ),
    r'$b=%.0f   (\sigma_{b}=3.8)$' % (b, )))
props = dict(boxstyle='round', facecolor='gray', alpha=0.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)
#plt.title("Channel No. versus Energy")
B = [312.86, 285.71, 153.46, 125.73, 100.83, 75.85]
for xy in zip(B, energies):
    ax.annotate(' (%s, %s keV)' % xy, xy=xy, textcoords='data')
plt.xlim([0,500])
plt.ylim([0,800])
plt.xlabel("Channel No.")
plt.ylabel("Energy (keV)")
plt.savefig("Linear_fit.png", bbox_inches="tight", pad_inches=0)
#plt.show()


###Create new plot with energy on x axis, using linear fit data
x = np.linspace(44, 2045.78, 1011)  #Define linear space that incorporates linear fit parameters
plt.figure(figsize=(10,7))
plt.plot(x, Cesium, linewidth=1, color = '#ff474c', label='Cs-137')
plt.plot(x, Radium, linewidth=1, color = '#2c6fbb', linestyle = '--', label='Ra-226')
plt.plot(x, Unknown, linewidth=2, color = '#10a674', linestyle = '-.', label='Sample No. 3')
plt.xlim([0,2046])
plt.ylim([0,2500])
#plt.title('Gamma Spectra (Calibrated)')
plt.xlabel("Energy (keV)")
plt.ylabel("Count")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("gamma_spectra_calibrated.png", bbox_inches="tight", pad_inches=0)
#plt.show()
