import calibration
import calibration_take_2
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp




#Calculate errors for Channel No. vs. Energy linear fit
def Sigma_y(A, B, x, y):
    sum = 0.0
    for i in range( len(y) ):
            sum += (y[i] - A - B*x[i])**2
    return np.sqrt( (1/ (len(y)-2.) )  *  sum )

sigma_y = Sigma_y(calibration.b, calibration.m, calibration.energy_calibration1[1], calibration.energies)

#Calculate chi^2 uncertainty

#Calculate uncertainty in parameters m and b in y = mx + b
def Delta(x):
    sum1 = 0.0
    sum2 = 0.0
    for i in range( len(x) ):
        sum1 += x[i]**2
        sum2 += (x[i])
    return len(x)*sum1-(sum2)**2

def Sigma_b(sig_y, x, delta):
    sum = 0.0
    for i in range (  len(x)):
        sum += x[i]**2
    return sig_y* np.sqrt( sum / delta )

def Sigma_m(sig_y, x, delta):
    return sig_y* np.sqrt( len(x) / delta )

delta = Delta(calibration.energy_calibration1[1])
sigma_b = Sigma_b(sigma_y, calibration.energy_calibration1[1], delta)
sigma_m = Sigma_m(sigma_y, calibration.energy_calibration1[1], delta)

print sigma_b
print sigma_m
#print calibration.m


###Calculate relative resolution of the detector for each photopeak and compute their average
def res(sigma, E):
    res = 0.0
    sum = 0.0
    for i in range ( len(calibration.energy_calibration1[0])):
        res = (2.*abs(sigma[i])*100.*calibration.m)/E[i]
#        print res
#        resolutions.fill(res())
        sum += res
    print (sum/ [len(calibration.energy_calibration1[0])])


#print res(calibration.energy_calibration1[0], calibration.energies)

###Calculate error on detector resolution
def res_error(sigma_vars, sig_m, E):
    sum = 0.0
    for i in range( len(calibration.energies) ):
        sum += np.sqrt( ((200*calibration.m* np.sqrt(sigma_vars[0][i][0]) )/E[i])**2 + ( (200*calibration.energy_calibration1[0][i]* sigma_m)/E[i])**2 )
#        print sum
    return (sum/ len(calibration.energies) )

resolution_error =  res_error(calibration.sigma_variances, sigma_m, calibration.energies)

print resolution_error
