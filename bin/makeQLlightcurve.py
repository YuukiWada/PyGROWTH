#!/usr/bin/env python
import astropy.io.fits as fitsio
import matplotlib.pyplot as plt
import numpy as np
import math
import signal
import sys
signal.signal(signal.SIGINT, signal.SIG_DFL)

# reading argiment parameters
if len(sys.argv)<6:
    print("Usage: makeQLspectrum.py <input file> <adc channel> <bin width> <lower channel> <upper channel>")
    exit()

input_file = sys.argv[1]
adc_channel = int(sys.argv[2])
bin_width = float(sys.argv[3])
threshold = [float(sys.argv[4]), float(sys.argv[5])]

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
bin_num = int(math.ceil(duration/bin_width))
bin_max = bin_num*bin_width
event = event[event["boardIndexAndChannel"]==adc_channel]
mask = (event["phaMax"]-event["phaMin"]>=threshold[0])&(event["phaMax"]-event["phaMin"]<=threshold[1])
event = event[mask]
time_tag = event["timeTag"]
time_tag[time_tag-time_tag[0]<-1*2**20]+=2**40
seconds = (time_tag-time_tag[0])/clock
print("Duration: {} sec".format(round(duration,1)))

# plot histogram
hist, bins = np.histogram(seconds,bins=bin_num,range=(0,bin_max))
ehist = hist**0.5/bin_width
hist = hist/bin_width
ebins = (bins[1:]-bins[0:-1])/2
bins = (bins[0:-1]+bins[1:])/2

plt.errorbar(bins, hist, xerr=ebins, yerr=ehist, capsize=0, linestyle="None", ecolor="black")
plt.xlabel("Seconds")
plt.ylabel("Count Rate (counts/s)")
plt.show()
