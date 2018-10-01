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
    return (sum/ [len(calibration.energy_calibration1[0])])


resolution = res(calibration.energy_calibration1[0], calibration.energies)

print "Detector resolution = {resolution}".format(resolution=resolution)

###Calculate error on detector resolution
def res_error(sigma_vars, sig_m, E):
    sum = 0.0
    for i in range( len(calibration.energies) ):
        sum += np.sqrt( ((200*calibration.m* np.sqrt(sigma_vars[0][i][0]) )/E[i])**2 + ( (200*calibration.energy_calibration1[0][i]* sigma_m)/E[i])**2 )
#        print sum
    return (sum/ len(calibration.energies) )
resolution_error =  res_error(calibration.sigma_variances, sigma_m, calibration.energies)
print "Detector resolution error = {resolution_error}".format(resolution_error=resolution_error)


###Find the average of the channel No. values
def avg_x(x):
    sum = 0.0
    for i in range ( len(calibration.energy_calibration1[1]) ):
        sum+= x[i]
    return( sum / len(calibration.energy_calibration1[1]))


channel_avg = avg_x(calibration.energy_calibration1[1])
#print "The channel avg is {channel_avg}".format(channel_avg=channel_avg)




###Use fitting results to find photopeak error plus uncertainty

def Energy(A, B, x):
    return(B*(x) + A)

def delE(B, x, delA, delB, delx, x0):
    return ( np.sqrt(delA**2 + ((x-x0)*(delB))**2 + (B*delx)**2 ))

energy_1 = Energy(calibration.b, calibration.m, calibration.unknown_params_1[1])
energy_2 = Energy(calibration.b, calibration.m, calibration.unknown_params_2[1])
del_E_1 = delE(calibration.m, calibration.unknown_params_1[1], sigma_b, sigma_m, calibration.unknown_params_1[0], channel_avg)
del_E_2 = delE(calibration.m, calibration.unknown_params_2[1], sigma_b, sigma_m, calibration.unknown_params_2[0], channel_avg)
print "Energy of unknown photopeak 1 is {energy_1}".format(energy_1=energy_1)
print "Energy of unknown photopeak 2 is {energy_2}".format(energy_2=energy_2)
print "The error for unknown 1 is {del_E_1}".format(del_E_1=del_E_1)
print "The error for unknown 2 is {del_E_2}".format(del_E_2=del_E_2)
