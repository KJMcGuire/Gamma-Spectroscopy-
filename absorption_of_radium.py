import data_arrays_take_2
import calibration_take_2
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

N = 1011
x0 = np.arange(N)
x = np.linspace(47, 1887.02, 1011)  #Define Energy space that incorporates linear fit parameters
plt.rcParams.update({'font.size': 15})


###Define an exponential to fit the background
def exponential(x, A, b):
    return A*np.exp(-x*b)

###Subtract global and local background from Ra-226 attenuated and unattenuated spectra
Radium_1mm_lead = data_arrays_take_2.sp3-data_arrays_take_2.sp0
Radium = data_arrays_take_2.sp2-data_arrays_take_2.sp0
Radium_7mm_lead = data_arrays_take_2.sp3_1-data_arrays_take_2.sp0

Radium_1mm_lead_clean = Radium_1mm_lead - exponential(x0, 3500, .018) -30
Radium_clean = Radium  - exponential(x0, 3500, .018) -20
Radium_7mm_lead_clean = Radium_7mm_lead - exponential(x0, 3500, .018) -30



###Create data slices for curve fitting and error analysis
peak = slice(280, 330, 1)
x1 = np.linspace(280, 330, 50)

Radium_cut = Radium_clean[peak]
Radium_cut_1mm = Radium_1mm_lead_clean[peak]
Radium_cut_7mm = Radium_7mm_lead_clean[peak]

#print Radium_cut_7mm


###Define a gaussian curve
def gaussian(x, sigma, X, A):
    return A*np.exp( (-(x-X)**2)/(2*sigma**2))

###Perform fit of 609 keV peak
ra_params1, ra_param_covar1 = curve_fit(gaussian, x1, Radium_clean[peak], p0 = [90, 310.,1.])
ra_params2, ra_param_covar2 = curve_fit(gaussian, x1, Radium_1mm_lead_clean[peak], p0 = [100, 310.,1.])
ra_params3, ra_param_covar3 = curve_fit(gaussian, x1, Radium_7mm_lead_clean[peak], p0 = [100, 310.,1.])

#print ra_params1
#print ra_params2
#print ra_params3

###Solve for counts under the photopeak fits
sum1=0
sum2=0
sum3=0
for i in range(peak.start, peak.stop + 1):
    sum1 += gaussian(x0[i], *ra_params1)   #Unattenuated
    sum2 += gaussian(x0[i], *ra_params2)   #1mm
    sum3 += gaussian(x0[i], *ra_params3)   #7mm

print sum1, sum2, sum3

###Errors for curve fits
perr_1 = np.sqrt(np.diag(abs(ra_params1)))   #Gives 1 standard deviation error of fit parameters
perr_2 = np.sqrt(np.diag(abs(ra_params2)))
perr_3 = np.sqrt(np.diag(abs(ra_params3)))
print perr_1, perr_2, perr_3
print ra_params1, ra_params2, ra_params3

#print perr

###Calculate Chi^2 of fits

numerator = 0
denominator = 0
sum = 0
for i in range( len(Radium_cut_7mm) ):
    numerator = ( Radium_cut_7mm[i] - gaussian(x1[i], *ra_params3))**2
    denominator = gaussian(x1[i], *ra_params3 )
    sum+= (numerator/denominator)


Chi_2_7mm = sum
reduced_Chi_2_7mm = sum/(len(Radium_cut_7mm)-1)
#print "Chi^2 7mm = {Chi_2_7mm}".format(Chi_2_7mm=Chi_2_7mm)
#print "Reduced Chi^2 7mm = {reduced_Chi_2_7mm}".format(reduced_Chi_2_7mm=reduced_Chi_2_7mm)



numerator = 0
denominator = 0
sum = 0
for i in range( len(Radium_cut) ):
    numerator = ( Radium_cut[i] - gaussian(x1[i], *ra_params1))**2
    denominator = gaussian(x1[i], *ra_params1 )
    sum+= (numerator/denominator)



Chi_2_unattenuated = sum
reduced_Chi_2_unattenuated = sum/(len(Radium_cut)-1)
#print "Chi^2 unattenuated = {Chi_2_unattenuated}".format(Chi_2_unattenuated=Chi_2_unattenuated)
#print "Reduced Chi^2 unattenuated = {reduced_Chi_2_unattenuated}".format(reduced_Chi_2_unattenuated=reduced_Chi_2_unattenuated)

numerator = 0
denominator = 0
sum = 0
for i in range( len(Radium_cut_1mm) ):
    numerator = ( Radium_cut_1mm[i] - gaussian(x1[i], *ra_params2))**2
    denominator = gaussian(x1[i], *ra_params2 )
    sum+= (numerator/denominator)

Chi_2_1mm = sum
reduced_Chi_2_1mm = sum/(len(Radium_cut_1mm)-1)
#print "Chi^2 1mm = {Chi_2_1mm}".format(Chi_2_1mm=Chi_2_1mm)
#print "Reduced Chi^2 1mm = {reduced_Chi_2_1mm}".format(reduced_Chi_2_1mm=reduced_Chi_2_1mm)


###Calculate uncertainty of counts
def count_uncert(sigma, A, dsigma, dA):
    delta_c = np.sqrt(  np.pi*2*(sigma**2)*(dA)**2     +  2*np.pi*(A*dsigma)**2    )
    return delta_c

count_uncertainty_1 = count_uncert( abs(ra_params1[0]), ra_params1[1], perr_1[0,0], perr_1[1,1])
count_uncertainty_2 = count_uncert( abs(ra_params2[0]), ra_params2[1], perr_2[0,0], perr_2[1,1])
count_uncertainty_3 = count_uncert( abs(ra_params3[0]), ra_params3[1], perr_3[0,0], perr_3[1,1])
print count_uncertainty_1, count_uncertainty_2, count_uncertainty_3

#Calculate the uncertainty in lead absorption coefficient

def del_mu(I, I0, x, dI, dI0, dx):
    np.sqrt( (dI0/(x*I0))**2  + (dI/(x*I))**2   + (( np.log(I0/I) * dx)/(x**2))**2    )

delta_mu_1_25mm = del_mu(sum2, sum1, 1.25, 2570, 2578, 0.5)
#print delta_mu_1_25mm


#Plot 609keV photopeaks with fits and count no.
plt.figure(figsize=(10,7))
#plt.plot(x, exponential(x0, 3500, .018) + 40)
plt.plot(x, Radium_clean, label="Unattenuated, \tilde{\chi}^2 = 5.1", color = '#ff474c')
plt.plot(x, gaussian(x0, *ra_params1), color = '#ff474c', linestyle = '-.')
plt.plot(x, Radium_1mm_lead_clean, color = '#2c6fbb', label="1.25 mm attenuator, $\tilde{\chi} = 4.0" )
plt.plot(x, gaussian(x0,  *ra_params2), color = '#2c6fbb', linestyle = '-.')
plt.plot(x, Radium_7mm_lead_clean, color = '#fec615', label="7 mm attenuator, \tilde{\chi}^2 = 2.3")
plt.plot(x, gaussian(x0,  *ra_params3), color = '#fec615', linestyle = '-.')
plt.annotate('17,700 ', xy=(609, 490), xytext=(609, 490))
plt.annotate('14,800', xy=(609, 390), xytext=(609, 390))
plt.annotate('8,400', xy=(609, 190), xytext=(609, 190))
#plt.xlim([0, 1000])
plt.xlim([560,660])
plt.ylim([-10,800])
#plt.title("Lead attenuation of Ra-226 609 keV photopeak")
plt.xlabel("Energy (keV)")
plt.ylabel("Count")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("Lead_attenuated_radium.png", bbox_inches="tight", pad_inches=0)
plt.show()
