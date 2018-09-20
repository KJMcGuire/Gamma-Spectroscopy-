import data_arrays
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit


N = 1011
x = np.arange(N)

###Subtract global background from known sample arrays
Cesium = data_arrays.sp1-data_arrays.sp0
Radium = data_arrays.sp2-data_arrays.sp0

###Create list of known photopeak energies
energies = [662, 609, 352, 295, 242, 186]

###Create data slices for curve fitting
cesium_s1 = slice(0, 280, 1)  #Slice skips the first item, so this slice will take data lines 1-280
x1 = np.linspace(1, 280, 280)
cesium_s2 = slice(280, 348, 1)
x2 = np.linspace(281, 348, 68)
cesium_s3 = slice(347, 1011, 1)
x3 = np.linspace(348, 1011, 664)
radium_s1 = slice(257, 314, 1)
ra_x1 = np.linspace(258, 314, 57)
radium_s2 = slice(138, 175, 1)
ra_x2 = np.linspace(139, 175, 37)
radium_s3 = slice(113, 138, 1)
ra_x3 = np.linspace(114, 138, 25)
radium_s4 = slice(90, 112, 1)
ra_x4 = np.linspace(91, 112, 22)
radium_s5 = slice(68, 85, 1)
ra_x5 = np.linspace(69, 85, 17)


###Define a gaussian curve plus exponential term
def gaussian(x, sigma, w, c, a, d):
    return (a * np.exp( (-(x-w)**2)/(2*sigma**2)) + d*np.exp(-c*((x))))


###Define linear function for energy calibration fit
def linear(X, m, b):
    return m*X + b


###Perform curve fitting of photopeaks
init_vals_ce2 = [1, 320, 1, 1, 1]   #Give curve-fitting function initial guesses for parameters
init_vals_ra1 = [210, 340, 1, 210, 1]
init_vals_ra2 = [130, 220, 1, 130, 1]
init_vals_ra3 = [110, 140, 1, 110, 1]
init_vals_ra4 = [90, 115, 1, 90, 1]
init_vals_ra5 = [65, 95, 1, 65, 1]

ce_params2, ce_param_covars2 = curve_fit(gaussian, x2, Cesium[cesium_s2], p0=init_vals_ce2)
ra_params1, ra_param_covar1 = curve_fit(gaussian, ra_x1, Radium[radium_s1], p0 = init_vals_ra1)
ra_params2, ra_param_covar2 = curve_fit(gaussian, ra_x2, Radium[radium_s2], p0 = init_vals_ra2)
ra_params3, ra_param_covar3 = curve_fit(gaussian, ra_x3, Radium[radium_s3], p0 = init_vals_ra3)
ra_params4, ra_param_covar4 = curve_fit(gaussian, ra_x4, Radium[radium_s4], p0 = init_vals_ra4)
ra_params5, ra_param_covar5 = curve_fit(gaussian, ra_x5, Radium[radium_s5], p0 = init_vals_ra5)


###Create lists of sigmas and channel numbers: the List Comprehension Way
energy_calibration1 = [t for t in zip(ce_params2, ra_params1, ra_params2, ra_params3, ra_params4, ra_params5)[:2]]

###The Longer Imperative Way
energy_calibration2 = []
times = 0
for t in zip(ce_params2, ra_params1, ra_params2, ra_params3, ra_params4, ra_params5):
    times+=1
    if times > 2:
        break
    energy_calibration2.append(t)


###Perform curve fitting of energy versus channel no.
init_vals_linear = [10,0]
linear_params, linear_covars = curve_fit(linear, energy_calibration1[1], energies, p0 = init_vals_linear)

print(linear(energy_calibration1[1], *linear_params))
###Plot spectra with fits
plt.figure(figsize=(10,7))
plt.plot(x, Cesium, label='Cs-137')
plt.plot(x, Radium, label='Ra-226')
plt.plot(x2, gaussian(x2, *ce_params2))
plt.plot(ra_x1, gaussian(ra_x1, *ra_params1))
plt.plot(ra_x2, gaussian(ra_x2, *ra_params2))
plt.plot(ra_x3, gaussian(ra_x3, *ra_params3))
plt.plot(ra_x4, gaussian(ra_x4, *ra_params4))
plt.plot(ra_x5, gaussian(ra_x5, *ra_params5))
plt.title("Known Spectra")
plt.xlim([0,600])
plt.ylim([0,3500])
#plt.legend()
plt.show()


plt.figure(figsize=(10,7))
plt.plot(energy_calibration1[1], linear(energy_calibration1[1], *linear_params))
plt.scatter(energy_calibration1[1], energies)
plt.xlim([0,600])
plt.ylim([0,3500])
#plt.legend()
plt.show()
