import data_arrays_take_2
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import numpy as np

N = 1011
x = np.arange(N)

#cm.register_cmap(name = "none", cmap = 'Dark2', lut=5)

plt.figure(figsize=(10,7))
plt.plot(x, data_arrays_take_2.sp3, label="1 mm", color = '#fac205')
plt.plot(x, data_arrays_take_2.sp3_1, label="7 mm", linestyle = '--', color = '#02ab2e')
plt.plot(x, data_arrays_take_2.sp2, label="Unattenuated", color = '#047495')
plt.xlim([0,1011])
plt.ylim([0,4000])
plt.xlabel("Channel No.")
plt.ylabel("Count")
#plt.title("Visual comparison of Radium 226 and 'Ore B'")
plt.legend(loc=1, prop={'size': 15})
plt.savefig("Lead_attenuated_radium.png", bbox_inches="tight", pad_inches=0)
plt.show()
