import numpy as np
import matplotlib.pyplot as plt

N = 1011
x = np.arange(N)

#Create data arrays from spectra
sp0 = np.genfromtxt('Spectrum_0_background_200s.Spe',
                skip_header=12,
                max_rows=N)
sp1 = np.genfromtxt('Spectrum_1_cs137_200s.Spe',
                skip_header=12,
                max_rows=N)
sp2 = np.genfromtxt('Spectrum_2_radium_227_200s.Spe',
                skip_header=12,
                max_rows=N)
sp3 = np.genfromtxt('Spectrum_3_lead_attenuated_cs137_200s.Spe',
                skip_header=12,
                max_rows=N)
sp4 = np.genfromtxt('Spectrum_4_1_unknown_3.Spe',
                skip_header=12,
                max_rows=N)
sp4_1 = np.genfromtxt('Spectrum_4_ore_b.Spe',
                skip_header=12,
                max_rows=N)
sp5 = np.genfromtxt('Spectrum_5_background_300s.Spe',
                skip_header=12,
                max_rows=N)
