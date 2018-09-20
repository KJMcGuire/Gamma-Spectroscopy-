import numpy as np
import matplotlib.pyplot as plt

N = 1011
x = np.arange(N)

#Create data arrays from spectra
sp0 = np.genfromtxt('spec0_background.Spe',
                skip_header=12,
                max_rows=N)
sp1 = np.genfromtxt('spec1_cs137.Spe',
                skip_header=12,
                max_rows=N)
sp2 = np.genfromtxt('spec2_ra226.Spe',
                skip_header=12,
                max_rows=N)
sp3 = np.genfromtxt('spec3_ra226_attenuated.Spe',
                skip_header=12,
                max_rows=N)
sp3_1 = np.genfromtxt('spec3_ra226_attenuated2.Spe',
                skip_header=12,
                max_rows=N)
sp4 = np.genfromtxt('spec_4_sample3.Spe',
                skip_header=12,
                max_rows=N)
sp5 = np.genfromtxt('background_300s.Spe',
                skip_header=12,
                max_rows=N)
