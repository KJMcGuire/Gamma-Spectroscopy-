import data_arrays
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy.stats import linregress

plt.rcParams.update({'font.size': 15})

N = 1011
x = np.arange(N)

###Subtract global background from known sample arrays
Cesium = data_arrays.sp1-data_arrays.sp0
Radium = data_arrays.sp2-data_arrays.sp0
Unknown = data_arrays.sp4-data_arrays.sp0

###Create list of known photopeak energies
energies = [1120, 661.69, 609.312, 351.932, 295.224, 241.997, 187.1]

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
unknown_1 = slice(510, 600, 1)
un_x1 = np.linspace(511, 600, 79)
unknown_1 = slice(237, 290, 1)
un_x1 = np.linspace(238, 290, 53)
unknown_2 = slice(495, 575, 1)
un_x2 = np.linspace(496, 575, 80)
radium_s6 = slice(532, 587, 1)
ra_x6 = np.linspace(533, 587, 55)
#radium_x7 = slice(532, 590, 1)
#ra_x7 = np.linspace(533, 590, 58)

###Define a gaussian curve
def gaussian(x, sigma, X, A):
    return A*np.exp( (-(x-X)**2)/(2*sigma**2))

###Define an exponential to fit the background
def exponential(x, A, b, c):
    return A*np.exp(-x*b) + c

#Subtract local background from known spectra
Radium_clean = Radium - exponential(x, 3500, .0145, 45)
Cesium_clean = Cesium - exponential(x, 3500, .0145, 0)
Unknown_clean = Unknown -exponential(x, 3500, 0.0145, 0)

ce_params1c, ce_param_covar1c = curve_fit(gaussian, x2, Cesium_clean[cesium_s2], p0 = [100, 340.,1.])
ra_params1c, ra_param_covar1c = curve_fit(gaussian, ra_x1, Radium_clean[radium_s1], p0 = [100, 340.,1.])
ra_params2c, ra_param_covar2c = curve_fit(gaussian, ra_x2, Radium_clean[radium_s2], p0 = [130., 220.,1.])
ra_params3c, ra_param_covar3c = curve_fit(gaussian, ra_x3, Radium_clean[radium_s3], p0 = [110., 140.,1.])
ra_params4c, ra_param_covar4c = curve_fit(gaussian, ra_x4, Radium_clean[radium_s4], p0 = [90., 115.,1.])
ra_params5c, ra_param_covar5c = curve_fit(gaussian, ra_x5, Radium_clean[radium_s5], p0 = [65., 95.,1.])
ra_params6c, ra_param_covar6c = curve_fit(gaussian, ra_x6, Radium_clean[radium_s6], p0 = [65., 500.,1.])
unknown_params_1, unknown_covar_1 = curve_fit(gaussian, un_x1, Unknown[unknown_1], p0 = [65., 210,1.])
unknown_params_2, unknown_covar_2 = curve_fit(gaussian, un_x2, Unknown[unknown_2], p0 = [65., 550,1.])

print ra_params6c

###Errors for curve fits
err_ce1 = np.sqrt(np.diag(abs(ce_param_covar1c)))   #Gives 1 standard deviation error of fit parameters
err_ra1 = np.sqrt(np.diag(abs(ra_param_covar1c)))
err_ra2 = np.sqrt(np.diag(abs(ra_param_covar2c)))
err_ra3 = np.sqrt(np.diag(abs(ra_param_covar3c)))
err_ra4 = np.sqrt(np.diag(abs(ra_param_covar4c)))
err_ra5 = np.sqrt(np.diag(abs(ra_param_covar5c)))
err_ra6 = np.sqrt(np.diag(abs(ra_param_covar6c)))

#print err_ce1

###Create lists of Gaussian fit errors
Gaussian_fit_errors = [t for t in zip(err_ra6, err_ce1, err_ra1, err_ra2, err_ra3,err_ra4, err_ra5)[:3]]

print Gaussian_fit_errors

#print ra_params6c


#print perr
#print ce_params1c

plt.figure(figsize=(10,7))
plt.plot(x2, gaussian(x2, *ce_params1c), color = '#2c6fbb',  linestyle = '-.', label='Fit curves')
plt.plot(ra_x1, gaussian(ra_x1, *ra_params1c), color = '#2c6fbb', linestyle = '-.')
plt.plot(ra_x2, gaussian(ra_x2, *ra_params2c), color = '#2c6fbb', linestyle = '-.')
plt.plot(ra_x3, gaussian(ra_x3, *ra_params3c), color = '#2c6fbb', linestyle = '-.')
plt.plot(ra_x4, gaussian(ra_x4, *ra_params4c), color = '#2c6fbb', linestyle = '-.')
plt.plot(ra_x5, gaussian(ra_x5, *ra_params5c), color = '#2c6fbb', linestyle = '-.')
plt.plot(ra_x6, gaussian(ra_x6, *ra_params6c), color = '#2c6fbb', linestyle = '-.')

#plt.plot(un_x1, gaussian(un_x1, *unknown_params_1), color = '#02ab2e', linestyle = '-.')
#plt.plot(un_x2, gaussian(un_x2, *unknown_params_2), color = '#02ab2e', linestyle = '-.')
#plt.plot(un_x1, Unknown[unknown_1])
#plt.plot(x, Unknown, linewidth=1)
#plt.plot(x, exponential(x, 3500, 0.0145, 0))
plt.annotate('662 keV', xy=(320, 920), xytext=(320, 920))
plt.annotate('1120 keV', xy=(535, 240), xytext=(535, 240))
plt.annotate('609 keV', xy=(240, 920), xytext=(240, 920))
plt.annotate('352 keV', xy=(120, 1580), xytext=(120, 1580))
plt.annotate('295 keV', xy=(80, 1200), xytext=(80, 1200))
plt.annotate('242 keV', xy=(50, 700), xytext=(50, 700))
plt.annotate('197 keV', xy=(30, 490), xytext=(30, 490))
plt.plot(x, Radium_clean, linewidth=1, color = '#fec615', linestyle = '--', label='Ra-226')
plt.plot(x, Cesium_clean, linewidth=1, color = '#ff474c', label='Cs-137')
plt.xlim([0,600])
plt.ylim([0, 2500])
plt.xlabel("Channel No.")
plt.ylabel("Count")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("Gaussian_fit_to_photo_peaks_w_unknown.png", bbox_inches="tight", pad_inches=0)
plt.show()


###Create lists of sigmas and channel numbers: the List Comprehension Way
energy_calibration1 = [t for t in zip(ra_params6c, ce_params1c, ra_params1c, ra_params2c, ra_params3c, ra_params4c, ra_params5c)[:2]]

#print energy_calibration1[1]

###Create lists of delta_sigmas from gaussian fits
sigma_variances = [t for t in zip(ra_param_covar6c, ce_param_covar1c, ra_param_covar1c, ra_param_covar2c, ra_param_covar3c, ra_param_covar4c, ra_param_covar5c)[:1]]
full_width = [15.75*2, 9.84*2, 9.81*2, 6.20*2, 4.95*2, 3.72*2, 2.67*2]

#print sigma_variances[0][1][0]
#print ra_param_covar1c[0][0]

###Perform curve fitting of energy versus channel no.
m, b, r_value, p_value, std_err = linregress(energy_calibration1[1], energies)

print  "The slope is {m}".format(m=m)
print  "And the y-int is {b}".format(b=b)


#Plot linear fit
x = np.arange(50, 600, 1)
#plt.figure(figsize=(10,7))
fig, ax = plt.subplots()
#ax.scatter(energy_calibration1[1], energies)
ax.errorbar(energy_calibration1[1], energies, xerr=(full_width), capsize=5, fmt='o', color = '#2c6fbb' )
ax.plot((x), m*(x) + b)

###Calculate Chi^2 of linear fit
def chi_2_linear(x_data, y_data, slope, intercept):
    numerator = 0
    denominator = 0
    sum = 0
    for i in range( len(y_data)):
        numerator = ( y_data[i] - (slope*x_data[i]+intercept))**2
        denominator = y_data[i]
        sum+= (numerator/denominator)
    return sum

chi_2_linear_fit = chi_2_linear(energy_calibration1[1], energies, m, b)
reduced_chi_2 = chi_2_linear_fit/(len(energies)-1)

#plt.annotate(' keV', xy=(580, 1000), xytext=(580, 1000))
textstr = '\n'.join((
    r'$m=%.2f \pm 0.02$' % (m, ),
    r'$b=%.0f \pm 6$' % (b, ),
    r'$\tilde{\chi}^2 =%.2f $' % (reduced_chi_2, )))
props = dict(boxstyle='round', facecolor='gray', alpha=0.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)
#plt.title("Channel No. versus Energy")
B = [561.24, 313.09, 285.65, 153.45, 125.73, 100.83, 75.85]
#for xy in zip(B, energies):
#    ax.annotate(' %s, %s keV' % xy, xy=xy, textcoords='data')
plt.xlim([0,650])
plt.ylim([0,1300])
plt.xlabel("Channel No.")
plt.ylabel("Energy (keV)")
plt.savefig("Linear_fit.png", bbox_inches="tight", pad_inches=0)
#plt.show()

Energy1 = unknown_params_1[1]*m+b
Energy2 = unknown_params_2[1]*m+b

#print Energy1, Energy2
#print unknown_params_1[0], unknown_params_2[0]
print  unknown_params_1[1]
print  unknown_params_2[1]

###Create new plot with energy on x axis, using linear fit data
x = np.linspace(42, 2053.89, 1011)  #Define linear space that incorporates linear fit parameters
plt.figure(figsize=(10,7))
plt.plot(x, Unknown_clean, linewidth=2, color = '#ff474c', linestyle = '-.', label='Sample No. 3')
plt.plot(x, Radium_clean, linewidth=1, color = '#2c6fbb', label='Ra-226')
plt.plot(x, Cesium_clean, linewidth=2, color = '#fec615', linestyle = '--', label='Cs-137')
plt.annotate("560$\pm$6.0 keV",  xy=(400, 930))
plt.annotate("1068$\pm$8.7 keV",  xy=(950, 340))
plt.xlim([0,1500])
plt.ylim([0,1500])
#plt.title('Gamma Spectra (Calibrated)')
plt.xlabel("Energy (keV)")
plt.ylabel("Count")
plt.legend(loc=1, prop={'size': 15})
#plt.savefig("gamma_spectra_calibrated.png", bbox_inches="tight", pad_inches=0)
plt.show()
