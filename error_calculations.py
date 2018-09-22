import calibration
#import calibration_take_2
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

#Calculate errors for Channel No. vs. Energy linear fit
def sigma_y(A, B, x, y):
    sum = 0.0
    for i in range( len(y) ):
            sum += (y[i] - A - B*x[i])**2
    return np.sqrt( (1/ (len(y)-2.) )  *  sum )

sigma_y = sigma_y(calibration.b, calibration.m, calibration.energy_calibration1[1], calibration.energies)

#Calculate chi^2 uncertainty

#Calculate uncertainty in parameters m and b in y = mx + b
def delta(x):
    sum1 = 0.0
    sum2 = 0.0
    for i in range( len(x) ):
        sum1 += x[i]**2
        sum2 += (x[i])
    return len(x)*sum1-(sum2)**2

def sigma_b(sig_y, x, delta):
    sum = 0.0
    for i in range (  len(x)):
        sum += x[i]**2
    return sig_y* np.sqrt( sum / delta )

def sigma_m(sig_y, x, delta):
    return sig_y* np.sqrt( len(x) / delta )

delta = delta(calibration.energy_calibration1[1])
sigma_b = sigma_b(sigma_y, calibration.energy_calibration1[1], delta)
sigma_m = sigma_m(sigma_y, calibration.energy_calibration1[1], delta)

print sigma_b
print sigma_m
#print calibration.m




#Calculate relative resolution of the detector for each photopeak and compute their average/ with uncertainty

#resolutions = np.array( [[len(calibration.energy_calibration1[1])], [len(calibration.energy_calibration1[1])]], np.int32)

def res(sigma, E):
    res = 0.0
    sum = 0.0
    for i in range ( len(calibration.energy_calibration1[0])):
        res = (2.*abs(sigma[i])*100.*calibration.m)/E[i]
        print res
#        resolutions.fill(res())
        sum += res
    print (sum/ [len(calibration.energy_calibration1[0])])



#resolution = res(calibration.energy_calibration1[0], calibration.energies)
#print resolution
