import data_arrays
import matplotlib.pyplot as plt
import numpy as np


N = 1011
x = np.arange(N)


plt.figure(figsize=(10,7))
plt.plot(x, data_arrays.sp0, label="Background", linestyle = '--', color = '#02ab2e')
plt.plot(x, data_arrays.sp1, label="Cs 137", linestyle = '-.', color = '#fac205')
plt.plot(x, data_arrays.sp2, label="Radium 226", color = '#047495')
plt.xlim([0,1011])
plt.ylim([0,3500])
plt.xlabel("Channel No.")
plt.ylabel("Count")
#plt.title("Known Spectra w/ Background")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("known_spectra_w_background_DRAFT.png", bbox_inches="tight", pad_inches=0)
plt.show()

plt.figure(figsize=(10,7))
plt.plot(x, data_arrays.sp2, label="Radium 226", color = '#fac205')
plt.plot(x, data_arrays.sp4_1, label="Ore B", linestyle = '--', color = '#02ab2e')
plt.xlim([0,1011])
plt.ylim([0,4000])
plt.xlabel("Channel No.")
plt.ylabel("Count")
#plt.title("Visual comparison of Radium 226 and 'Ore B'")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("radium_ore_b.png", bbox_inches="tight", pad_inches=0)
plt.show()
