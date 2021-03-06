import data_arrays_take_2
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy.stats import linregress
#from sympy.interactive import printing
#printing.init_printing(use_latex=True)
#from IPython.display import display, Math, Latex




N = 1011
x = np.arange(N)

###Subtract global background from known sample arrays
Cesium = data_arrays_take_2.sp1-data_arrays_take_2.sp0
Radium = data_arrays_take_2.sp2-data_arrays_take_2.sp0
Unknown = data_arrays_take_2.sp4-data_arrays_take_2.sp0

###Create list of known photopeak energies
energies = [661.69, 609.312, 351.932, 295.224, 241.997, 187.1]

###Create data slices for curve fitting
#cesium_s1 = slice(0, 281, 1)  #Slice skips the first item, so this slice will take data lines 1-280
#x1 = np.linspace(1, 281, 281)
cesium_s2 = slice(306, 378, 1)
x2 = np.linspace(307, 378, 72)
#cesium_s3 = slice(347, 1011, 1)
#x3 = np.linspace(348, 1011, 664)
radium_s1 = slice(276, 338, 1)
ra_x1 = np.linspace(277, 338, 62)
radium_s2 = slice(148, 179, 1)
ra_x2 = np.linspace(149, 179, 31)
radium_s3 = slice(120, 148, 1)
ra_x3 = np.linspace(121, 148, 28)
radium_s4 = slice(97, 120, 1)
ra_x4 = np.linspace(98, 120, 23)
radium_s5 = slice(70, 92, 1)
ra_x5 = np.linspace(71, 92, 22)

#plt.figure(figsize=(10,7))
#plt.plot(x2, Cesium[cesium_s2], linewidth=1)
#plt.plot(ra_x1, Radium[radium_s1], linewidth=1)
#plt.plot(ra_x2, Radium[radium_s2], linewidth=1)
#plt.plot(ra_x3, Radium[radium_s3], linewidth=1)
#plt.plot(ra_x4, Radium[radium_s4], linewidth=1)
#plt.plot(ra_x5, Radium[radium_s5], linewidth=1)
#plt.plot(x, Radium, linewidth=1, label='Ra-226')
#plt.plot(x, Cesium, linewidth=1, label='Cs-137')
#plt.xlim([0,600])
#plt.ylim([0,2500])
#plt.xlabel("Channel No.")
#plt.ylabel("Count")
#plt.show()

###Define a gaussian curve plus exponential term
def gaussian(x, sigma, X, A):
   return (A/(sigma*np.sqrt(2*np.pi))) * np.exp((-(x-X)**2)/(2*sigma**2))

def gaussian2(x, c, X, A):
    return A*np.exp(c*(x-X)**2)

#def exponential(x, A, c):
#    return A*np.exp(-x*c)

###Perform curve fitting of photopeaks
ce_params2, ce_param_covars2 = curve_fit(gaussian, x2, Cesium[cesium_s2], p0= [1., 320., 1.])
ra_params1, ra_param_covar1 = curve_fit(gaussian, ra_x1, Radium[radium_s1], p0 = [210., 340., 1.])
ra_params2, ra_param_covar2 = curve_fit(gaussian, ra_x2, Radium[radium_s2], p0 = [130., 220., 1.])
ra_params3, ra_param_covar3 = curve_fit(gaussian, ra_x3, Radium[radium_s3], p0 = [110., 140., 1.])
ra_params4, ra_param_covar4 = curve_fit(gaussian, ra_x4, Radium[radium_s4], p0 = [90., 115., 1.])
ra_params5, ra_param_covar5 = curve_fit(gaussian, ra_x5, Radium[radium_s5], p0 = [90., 115., 1.])
#print (ce_params2)
#print (ra_params5)
#exp_fit_params, exp_fit_covars = curve_fit(exponential, ra_x3, Radium[radium_s3], p0 = [1., 1.])


###Create lists of sigmas and channel numbers: the List Comprehension Way
energy_calibration1 = [t for t in zip(ce_params2, ra_params1, ra_params2, ra_params3, ra_params4, ra_params5)[:2]]

print ra_params3
#print exp_fit_params
#print energy_calibration1

###The Longer Imperative Way
#energy_calibration2 = []
#times = 0
#for t in zip(ce_params2, ra_params1, ra_params2, ra_params3, ra_params4, ra_params5):
#    times+=1
#    if times > 2:
#        break
#    energy_calibration2.append(t)

###Perform curve fitting of energy versus channel no.
m, b, r_value, p_value, std_err = linregress(energy_calibration1[1], energies)

###Plot spectra with fits
plt.figure(figsize=(10,7))
#plt.plot(ra_x3, exponential(ra_x3,*exp_fit_params))
plt.plot(x, Cesium, linewidth=1, color = '#047495', label='Cs-137')
plt.plot(x, Radium, linewidth=2, color = '#fac205', linestyle = '--', label='Ra-226')
plt.plot(x2, gaussian(x2, *ce_params2), color = '#02ab2e',  linestyle = '-.', label='Fit curves', )
plt.plot(ra_x1, gaussian(ra_x1, *ra_params1), color = '#02ab2e', linestyle = '-.')
plt.annotate('662 keV', xy=(320, 720), xytext=(320, 820))
plt.annotate('609 keV', xy=(260, 770), xytext=(260, 870))
plt.annotate('352 keV', xy=(155, 1780), xytext=(155, 1880))
plt.annotate('295 keV', xy=(110, 1580), xytext=(110, 1680))
plt.annotate('242 keV', xy=(90, 1300), xytext=(90, 1400))
plt.annotate('197 keV', xy=(50, 1430), xytext=(50, 1530))
plt.plot(ra_x2, gaussian(ra_x2, *ra_params2), color = '#02ab2e', linestyle = '-.')
plt.plot(ra_x3, gaussian(ra_x3, *ra_params3), color = '#02ab2e', linestyle = '-.')
plt.plot(ra_x4, gaussian(ra_x4, *ra_params4), color = '#02ab2e', linestyle = '-.')
plt.plot(ra_x5, gaussian(ra_x5, *ra_params5), color = '#02ab2e', linestyle = '-.')
plt.title("Gaussian fit to photopeaks")
plt.xlim([0,600])
plt.ylim([0,2500])
plt.xlabel("Channel No.")
plt.ylabel("Count")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("Gaussian_fit_to_photo_peaks_take_2.png", bbox_inches="tight", pad_inches=0)
plt.show()

###Plot linear fit
x = np.arange(50, 350, 1)
plt.figure(figsize=(10,7))
fig, ax = plt.subplots()
ax.plot((x), m*(x) + b)
ax.scatter(energy_calibration1[1], energies, )
textstr = '\n'.join((
    r'$m=%.2f  (\sigma_{m}=.02)$' % (m, ),
    r'$b=%.0f   (\sigma_{b}=3.5)$' % (b, )))
props = dict(boxstyle='round', facecolor='gray', alpha=0.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)
plt.title("Channel No. versus Energy")
B = [337.825, 307.8098, 164.587, 134.0995, 105.902, 77.895, ]
for xy in zip(B, energies):
    ax.annotate(' (%s, %s keV)' % xy, xy=xy, textcoords='data')
plt.xlim([0,500])
plt.ylim([0,800])
plt.xlabel("Channel No.")
plt.ylabel("Energy (keV)")
plt.savefig("Linear_fit_take_2.png", bbox_inches="tight", pad_inches=0)
#plt.show()

###Create new plot with energy on x axis, using linear fit data
x = np.linspace(49, 1889.02, 1011)  #Define linear space that incorporates linear fit parameters
plt.figure(figsize=(10,7))
plt.plot(x, Cesium, linewidth=1, color = '#047495', label='Cs-137')
plt.plot(x, Radium, linewidth=2, color = '#fac205', linestyle = '--', label='Ra-226')
plt.plot(x, Unknown, linewidth=1, color = '#ae7181', linestyle = '-.', label='Sample No. 3')
plt.xlim([0,1046])
plt.ylim([0,2500])
plt.title('Gamma Spectra (Calibrated)')
plt.xlabel("Energy (keV)")
plt.ylabel("Count")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("gamma_spectra_calibrated_take_2.png", bbox_inches="tight", pad_inches=0)
plt.show()
