import data_arrays_take_2
import calibration_take_2
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

N = 1011
x0 = np.arange(N)
x = np.linspace(47, 1887.02, 1011)  #Define Energy space that incorporates linear fit parameters


###Subtract global background plus costant from Ra-226 attenuated and unattenuated spectra
Radium_1mm_lead = data_arrays_take_2.sp3-data_arrays_take_2.sp0 - 100
Radium = data_arrays_take_2.sp2-data_arrays_take_2.sp0 - 100
Radium_7mm_lead = data_arrays_take_2.sp3_1-data_arrays_take_2.sp0 - 100

peak = slice(280, 333, 1)
x1 = np.linspace(281, 333, 53)

###Define an exponential to fit the background
def exponential(x, A, b):
    return A*np.exp(-x*b)

###Define a gaussian curve
def gaussian(x, sigma, X, A):
    return A*np.exp( (-(x-X)**2)/(2*sigma**2))

###Perform fit of 609 keV peak
ra_params1, ra_param_covar1 = curve_fit(gaussian, x1, Radium[peak], p0 = [100, 310.,1.])
ra_params2, ra_param_covar2 = curve_fit(gaussian, x1, Radium_1mm_lead[peak], p0 = [100, 310.,1.])
ra_params3, ra_param_covar3 = curve_fit(gaussian, x1, Radium_7mm_lead[peak], p0 = [100, 310.,1.])

print ra_params1
print ra_params2
print ra_params3

###Solve for counts under the photopeak fits
sum1=0
sum2=0
sum3=0
for i in range(peak.start, peak.stop + 1):
    sum1 += gaussian(x0[i], *ra_params1)   #Unattenuated
    sum2 += gaussian(x0[i], *ra_params2)   #1mm
    sum3 += gaussian(x0[i], *ra_params3)   #7mm

print sum1, sum2, sum3

#Plot 609keV photopeaks with fits and count no.
plt.figure(figsize=(10,7))
plt.plot(x, Radium, label="Unattenuated", color = '#ff474c')
plt.plot(x, gaussian(x0, *ra_params1), color = '#ff474c', linestyle = '-.', label="Fit")
plt.plot(x, Radium_1mm_lead, color = '#2c6fbb', label="1 mm attenuator" )
plt.plot(x, gaussian(x0,  *ra_params2), color = '#2c6fbb', linestyle = '-.', label="Fit")
plt.plot(x, Radium_7mm_lead, color = '#fec615', label="7 mm attenuator")
plt.plot(x, gaussian(x0,  *ra_params3), color = '#fec615', linestyle = '-.', label="Fit")
plt.annotate('14,316', xy=(609, 490), xytext=(609, 490))
plt.annotate('11,945', xy=(609, 390), xytext=(609, 390))
plt.annotate('5,614', xy=(609, 190), xytext=(609, 190))
plt.xlim([549,669])
plt.ylim([0,1000])
plt.xlabel("Energy (keV)")
plt.ylabel("Count")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("Lead_attenuated_radium.png", bbox_inches="tight", pad_inches=0)
plt.show()
