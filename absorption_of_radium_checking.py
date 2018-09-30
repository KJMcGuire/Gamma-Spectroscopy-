import data_arrays_take_2
import calibration_take_2
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

N = 1011
x0 = np.arange(N)
x = np.linspace(47, 1897.13, 1011)  #Define Energy space that incorporates linear fit parameters
plt.rcParams.update({'font.size': 20})

###Subtract global and local background from Ra-226 attenuated and unattenuated spectra
Radium_1mm_lead = data_arrays_take_2.sp3-data_arrays_take_2.sp0
Radium = data_arrays_take_2.sp2-data_arrays_take_2.sp0
Radium_7mm_lead = data_arrays_take_2.sp3_1-data_arrays_take_2.sp0

Radium_1mm_lead_clean = Radium_1mm_lead + 0.35*x0 - 180
Radium_clean = Radium  + 0.35*x0 - 180
Radium_7mm_lead_clean = Radium_7mm_lead + 0.35*x0 - 180


###Define an exponential to fit the background
def exponential(x, A, b):
    return A*np.exp(-x*b)

###Create data slices for curve fitting and error analysis
peak = slice(280, 331, 1)  ###THIS AND NEXT LINE IS THE CORRECT WAY TO DO THESE SLICES!!!
x1 = np.linspace(280, 330, 51)

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

print "Counts under the curves are: {sum1}, {sum2}, {sum3}".format(sum1=sum1, sum2=sum2, sum3=sum3)




###Errors for curve fits
perr_1 = np.sqrt(np.diag(abs(ra_param_covar1)))   #Gives 1 standard deviation error of fit parameters
perr_2 = np.sqrt(np.diag(abs(ra_param_covar2)))
perr_3 = np.sqrt(np.diag(abs(ra_param_covar3)))
#print perr_1, perr_2, perr_3
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
print "Count uncertainties are: {count_uncertainty_1}, {count_uncertainty_2}, {count_uncertainty_3}".format(count_uncertainty_1=count_uncertainty_1, count_uncertainty_2=count_uncertainty_2, count_uncertainty_3=count_uncertainty_3)


###Plot absporber thickness versus count and fit to exponential function
count_uncertainties = [count_uncertainty_1, count_uncertainty_2, count_uncertainty_3]
total_counts = [sum1, sum2, sum3]
thicknesses = [0, 1, 7]

absorption_fit_coeff, absorption_fit_covars = curve_fit(exponential, thicknesses, total_counts, p0 = [18000.,1.])



###Calculate mass attenuation coefficient
def mass_attenuation(I, I0, x, rho):
    mu_rho = (1/(x*rho)) * np.log(I0/I)
    return mu_rho


mu_7 = mass_attenuation(sum3, sum1, .7, 1)
mu_1 = mass_attenuation(sum2, sum1, .125, 1)
mu_rho_7mm = mass_attenuation(sum3, sum1, .7, 11.34)
mu_rho_1mm = mass_attenuation(sum2, sum1, .125, 11.34)
print "The attenuation coefficient of lead w/ 7mm absorber is {0:.5}".format(mu_7)
print "The attenuation coefficient of lead w/ 1mm absorber is {0:.5}".format(mu_1)
print "The mass attenuation coefficient of lead w/ 7mm absorber is {0:.5}".format(mu_rho_7mm)
print "The mass attenuation coefficient of lead w/ 1mm absorber is {0:.5}".format(mu_rho_1mm)


###Calculate the uncertainty in mass attenuation coefficient

def del_mu_rho(I, I0, x, dI, dI0, dx, rho):
    delta_mu_rho = np.sqrt( (dI0/(x*I0*rho))**2  + (dI/(x*I*rho))**2   + (( np.log(I0/I) * dx)/(rho*x**2))**2    )
    return delta_mu_rho

delta_mu_7mm = del_mu_rho(sum3, sum1, .7, count_uncertainty_3, count_uncertainty_1, 0.05, 11.35)
print "The uncertainty in mass attenuation coeff. from 7mm absorber = {delta_mu_7mm}".format(delta_mu_7mm=delta_mu_7mm)

delta_mu_1mm = del_mu_rho(sum2, sum1, .125, count_uncertainty_2, count_uncertainty_1, 0.05, 11.35)
print "The uncertainty in mass attenuation coeff. from 1mm absorber = {delta_mu_1mm}".format(delta_mu_1mm=delta_mu_1mm)


###Plot 609keV photopeaks with fits and count no.
plt.figure(figsize=(10,7))
#plt.plot(x0, exponential(x0, 3500, .018) + 40)
#plt.plot(x, -0.35*x0 + 175)
plt.plot(x, Radium_clean, label="Unattenuated ($\\tilde{\chi}^{2}$ = 3.9)", color = '#ff474c')
plt.plot(x, gaussian(x0, *ra_params1), color = '#ff474c', linestyle = '-.')
plt.plot(x, Radium_1mm_lead_clean, color = '#2c6fbb', label="1.25 mm attenuator ($\\tilde{\chi}^{2}$ = 3.4)" )
plt.plot(x, gaussian(x0,  *ra_params2), color = '#2c6fbb', linestyle = '-.')
plt.plot(x, Radium_7mm_lead_clean, color = '#fec615', label="7 mm attenuator ($\\tilde{\chi}^{2}$ = 1.8)")
plt.plot(x, gaussian(x0,  *ra_params3), color = '#fec615', linestyle = '-.')
plt.annotate('15,600$\pm$300 ', xy=(609, 490), xytext=(609, 490), fontweight='bold', fontsize = 15)
plt.annotate('13,300$\pm$300', xy=(609, 390), xytext=(609, 390), fontweight='bold', fontsize = 15)
plt.annotate('6,900$\pm$200', xy=(609, 190), xytext=(609, 190), fontweight='bold', fontsize = 15)
#plt.xlim([0, 1000])
plt.xlim([550,670])
plt.ylim([-10,800])
#plt.title("Lead attenuation of Ra-226 609 keV photopeak")
plt.xlabel("Energy (keV)")
plt.ylabel("Count")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("Lead_attenuated_radium.png", bbox_inches="tight", pad_inches=0)
plt.show()

###Fit the NIST data for mass attenuation coefficient to an exponential

NIST_energies = [200, 300, 400, 500, 600, 800, 1000, 1250]
NIST_coeff  = [.9985, .4031, .2323, .1614, .1248, .08870, .07102, .05876]

NIST_params, NIST_covars = curve_fit(exponential, NIST_energies, NIST_coeff, p0 = [.1, .00000001])

print NIST_params

N = 5000

plt.figure(figsize=(10,7))
plt.plot(x0, exponential(x0, *NIST_params), label = "Fit curve ($\widetilde{\chi^2}$ = )")
plt.scatter(NIST_energies, NIST_coeff, color = '#2c6fbb')
plt.xlabel("$\mu$/$\\rho$ (cm$^{2}$/g)")
plt.ylabel("Energy (keV)")
plt.ylim([0,1.2])
plt.xlim([0, 1400])
#plt.legend(loc=1, prop={'size': 15})
plt.savefig("Mass_attenuation_lead_nist_data.png", bbox_inches="tight", pad_inches=0)
plt.show()


###Chi^2 of exponential fit
numerator = 0
denominator = 0
sum = 0
for i in range( len(NIST_coeff) ):
    numerator = ( NIST_coeff[i] - exponential(NIST_energies[i], *NIST_params))**2
    denominator = exponential(x1[i], *NIST_params )
    sum+= (numerator/denominator)

Chi_2_NIST = sum
reduced_Chi_2_NIST = Chi_2_NIST/ (len(NIST_coeff) - 1)
print "Chi^2 for NIST fit = {Chi_2_NIST}".format(Chi_2_NIST=Chi_2_NIST)
print "Reduced Chi^2 for NIST fit = {reduced_Chi_2_NIST}".format(reduced_Chi_2_NIST=reduced_Chi_2_NIST)
