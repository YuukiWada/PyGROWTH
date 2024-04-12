#!/usr/bin/env python
import astropy.io.fits as fitsio
import matplotlib.pyplot as plt
import numpy as np
import math
import signal
import sys
signal.signal(signal.SIGINT, signal.SIG_DFL)

# reading argiment parameters
if len(sys.argv)<2:
    print("Usage: makeQLspectrum.py <input file> <adc channel> [bin width (MeV)]")
    exit()

input_file = sys.argv[1]
adc_channel = int(sys.argv[2])
if len(sys.argv)>3:
    bin_width = float(sys.argv[3])
    if (bin_width>10) or (bin_width<+0):
        print("The bin width should be more than 0 MeV and less than 10 MeV")
        bin_width = 0.05
else:
    bin_width = 0.05

bin_num = int(math.ceil(25/bin_width))
bin_max = bin_width*bin_num

# reading fits file
fits_file = fitsio.open(input_file)
event = fits_file[1].data
event = event[~(event["phaMax"]>10000)]
duration = event["unixTime"][-1]-event["unixTime"][0]
event = event[event["boardIndexAndChannel"]==adc_channel]
print("Duration: {} sec".format(round(duration,1)))

energy = (event["energy"]+(np.random.default_rng().random(event["energy"].shape[0])-0.5)*fits_file[1].header["BINW_CH{}".format(adc_channel)])/1000

# plot histogram
hist, bins = np.histogram(energy,bins=bin_num,range=(0,bin_max))
ehist = hist**0.5/(bin_width*duration)
hist = hist/(bin_width*duration)
ebins = (bins[1:]-bins[0:-1])/2
bins = (bins[0:-1]+bins[1:])/2

plt.errorbar(bins, hist, xerr=ebins, yerr=ehist, capsize=0, linestyle="None", ecolor="black")
plt.xlabel("Energy (MeV)")
plt.ylabel("Spectrum (counts/s/MeV)")
plt.show()
