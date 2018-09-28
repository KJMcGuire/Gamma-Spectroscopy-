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

Radium_1mm_lead_clean = Radium_1mm_lead +0.35*x0 - 200#- exponential(x0, 3500, .018) -30
Radium_clean = Radium  +0.35*x0 - 200#- exponential(x0, 3500, .018) -20
Radium_7mm_lead_clean = Radium_7mm_lead +0.35*x0 - 175#- exponential(x0, 3500, .018) -30



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
perr_1 = np.sqrt(np.diag(abs(ra_param_covar1)))   #Gives 1 standard deviation error of fit parameters
perr_2 = np.sqrt(np.diag(abs(ra_param_covar2)))
perr_3 = np.sqrt(np.diag(abs(ra_param_covar3)))
print perr_1, perr_2, perr_3
#print ra_params1, ra_params2, ra_params3

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
print "Reduced Chi^2 7mm = {reduced_Chi_2_7mm}".format(reduced_Chi_2_7mm=reduced_Chi_2_7mm)



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
print "Reduced Chi^2 unattenuated = {reduced_Chi_2_unattenuated}".format(reduced_Chi_2_unattenuated=reduced_Chi_2_unattenuated)

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
print "Reduced Chi^2 1mm = {reduced_Chi_2_1mm}".format(reduced_Chi_2_1mm=reduced_Chi_2_1mm)


###Calculate uncertainty of counts
def count_uncert(sigma, A, dsigma, dA):
    delta_c = np.sqrt(  np.pi*2*(sigma**2)*(dA)**2     +  2*np.pi*(A*dsigma)**2    )
    return delta_c

count_uncertainty_1 = count_uncert( abs(ra_params1[0]), ra_params1[2], perr_1[0], perr_1[2])
count_uncertainty_2 = count_uncert( abs(ra_params2[0]), ra_params2[2], perr_2[0], perr_2[2])
count_uncertainty_3 = count_uncert( abs(ra_params3[0]), ra_params3[2], perr_3[0], perr_3[2])
print count_uncertainty_1, count_uncertainty_2, count_uncertainty_3


###Plot absporber thickness versus count and fit to exponential function
count_uncertainties = [count_uncertainty_1, count_uncertainty_2, count_uncertainty_3]
total_counts = [sum1, sum2, sum3]
thicknesses = [0, 1.25, 7]

absorption_fit_coeff, absorption_fit_covars = curve_fit(exponential, thicknesses, total_counts, p0 = [18000.,1.])



###Calculate lead absorption coefficient
def lead_attenuation(I, I0, x):
    mu = (1/x) * np.log(I0/I)
    return mu


mu_7mm = lead_attenuation(sum3, sum1, 7)
print "The attenuation coefficient w/ 7mm absorner is {0:.5}".format(mu_7mm)

###Calculate the uncertainty in lead absorption coefficient

def del_mu(I, I0, x, dI, dI0, dx):
    delta_mu = np.sqrt( (dI0/(x*I0))**2  + (dI/(x*I))**2   + (( np.log(I0/I) * dx)/(x**2))**2    )
    return delta_mu

delta_mu_1_25mm = del_mu(sum2, sum1, 1.25, count_uncertainty_2, count_uncertainty_1, 0.5)
print "The uncertainty in mu from 1.25mm absorber = {delta_mu_1_25mm}".format(delta_mu_1_25mm=delta_mu_1_25mm)

delta_mu_7mm = del_mu(sum3, sum1, 7, count_uncertainty_3, count_uncertainty_1, 0.5)
print "The uncertainty in mu from 7mm absorber = {delta_mu_7mm}".format(delta_mu_7mm=delta_mu_7mm)


###Plot 609keV photopeaks with fits and count no.
plt.figure(figsize=(10,7))
#plt.plot(x0, exponential(x0, 3500, .018) + 40)
#plt.plot(x0, -0.35*x0 + 175)
plt.plot(x, Radium_clean, label="Unattenuated ($\widetilde{\chi^2}$ = 3.2)", color = '#ff474c')
plt.plot(x, gaussian(x0, *ra_params1), color = '#ff474c', linestyle = '-.')
plt.plot(x, Radium_1mm_lead_clean, color = '#2c6fbb', label="1.25 mm attenuator ($\widetilde{\chi^2}$ = 3.2)" )
plt.plot(x, gaussian(x0,  *ra_params2), color = '#2c6fbb', linestyle = '-.')
plt.plot(x, Radium_7mm_lead_clean, color = '#fec615', label="7 mm attenuator ($\widetilde{\chi^2}$ = 1.8)")
plt.plot(x, gaussian(x0,  *ra_params3), color = '#fec615', linestyle = '-.')
plt.annotate('14,800$\pm$300 ', xy=(609, 490), xytext=(609, 490), fontweight='bold', fontsize = 15)
plt.annotate('12,400$\pm$300', xy=(609, 390), xytext=(609, 390), fontweight='bold', fontsize = 15)
plt.annotate('7,200$\pm$200', xy=(609, 190), xytext=(609, 190), fontweight='bold', fontsize = 15)
#plt.xlim([0, 1000])
plt.xlim([540,680])
plt.ylim([-10,900])
#plt.title("Lead attenuation of Ra-226 609 keV photopeak")
plt.xlabel("Energy (keV)")
plt.ylabel("Count")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("Lead_attenuated_radium.png", bbox_inches="tight", pad_inches=0)
#plt.show()



###Chi^2 of exponential fit
Chi_2_exp_fit = (14819-exponential(0,  *absorption_fit_coeff) )**2/exponential(0,  *absorption_fit_coeff)
+(12436-exponential(1.25,  *absorption_fit_coeff) )**2/exponential(1.25,  *absorption_fit_coeff)  + (7228-exponential(7,  *absorption_fit_coeff) )**2/exponential(7,  *absorption_fit_coeff)

reduced_Chi_2_exp_fit = Chi_2_exp_fit/2
#print "Chi^2 7mm = {Chi_2_7mm}".format(Chi_2_7mm=Chi_2_7mm)
print "Reduced Chi^2 of exponential fit = {reduced_Chi_2_exp_fit}".format(reduced_Chi_2_exp_fit=reduced_Chi_2_exp_fit)

full_width = [2*300, 2*300, 2*200]



#plt.figure(figsize=(10,7))
fig, ax = plt.subplots()
plt.plot(x0, exponential(x0, *absorption_fit_coeff), label = "Fitted curve ($\widetilde{\chi^2}$ = 2.4)")
ax.errorbar(thicknesses, total_counts, yerr=(full_width), capsize=5, fmt='o', color = '#2c6fbb')
plt.xlabel("Absorber Thickness (mm)")
plt.ylabel("Count")
plt.xlim([-1,10])
plt.ylim([0,18000])
plt.legend(loc=1, prop={'size': 15})
plt.savefig("Lead_attenuation_versus_thickness.png", bbox_inches="tight", pad_inches=0)
plt.show()
