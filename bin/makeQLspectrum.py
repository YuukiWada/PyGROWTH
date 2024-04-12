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
    print("Usage: makeQLspectrum.py <input file> <adc channel> [binning]")
    exit()

input_file = sys.argv[1]
adc_channel = int(sys.argv[2])
if len(sys.argv)>3:
    rebin = int(sys.argv[3])
    if (rebin>2048) or (rebin<1):
        print("The number of binning should be more than 0 and less than 2048")
        rebin = 1
    elif not math.log2(rebin).is_integer():
        print("The number of binning should be 2 to the power of n")
        rebin = 1
else:
    rebin = 1

bin_num = int(2048/rebin)

# reading fits file
fits_file = fitsio.open(input_file)
event = fits_file[1].data
event = event[~(event["phaMax"]>10000)]

clock = 1.0e8
start_count = event["timeTag"][0]
end_count = event["timeTag"][-1]
if start_count>end_count:
    end_count+=2**40
duration = float(end_count-start_count)/clock
event = event[event["boardIndexAndChannel"]==adc_channel]
print("Duration: {} sec".format(round(duration,1)))

# plot histogram
hist, bins = np.histogram(event["phaMax"]-event["phaMin"],bins=bin_num,range=(-0.5,2047.5))
ehist = hist**0.5/(rebin*duration)
hist = hist/(rebin*duration)
ebins = (bins[1:]-bins[0:-1])/2
bins = (bins[0:-1]+bins[1:])/2

plt.errorbar(bins, hist, xerr=ebins, yerr=ehist, capsize=0, linestyle="None", ecolor="black")
plt.yscale("log")
plt.xlabel("Channel")
plt.ylabel("Spectrum (counts/ch/s)")
plt.show()
