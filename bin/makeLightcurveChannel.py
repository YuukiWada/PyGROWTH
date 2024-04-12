#!/usr/bin/env python
from matplotlib import dates as mdates
import astropy.io.fits as fitsio
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import math
import signal
import sys
signal.signal(signal.SIGINT, signal.SIG_DFL)

# reading argiment parameters
if len(sys.argv)<2:
    print("Usage: makeQLspectrum.py <input file> <adc channel> <bin width (sec)> <lower channel> <upper channel>")
    exit()

input_file = sys.argv[1]
adc_channel = int(sys.argv[2])
bin_width = float(sys.argv[3])
threshold = [float(sys.argv[4]),float(sys.argv[5])]

# reading fits file
fits_file = fitsio.open(input_file)
event = fits_file[1].data
event = event[~(event["phaMax"]>10000)]
duration = event["unixTime"][-1]-event["unixTime"][0]
event = event[event["boardIndexAndChannel"]==adc_channel]
print("Duration: {} sec".format(round(duration,1)))

bin_min = math.floor(event["unixTime"][0]/bin_width)*bin_width
bin_max = math.ceil(event["unixTime"][-1]/bin_width)*bin_width
bin_num = int((bin_max-bin_min)/bin_width)
channel = event["phaMax"]-event["phaMin"]
unix_time = event["unixTime"][(channel>=threshold[0])&(channel<=threshold[1])]

# plot histogram
hist, bins = np.histogram(unix_time,bins=bin_num,range=(bin_min,bin_max))
ehist = hist**0.5/(bin_width)
hist = hist/(bin_width)
ebins = (bins[1:]-bins[0:-1])/2
bins = (bins[0:-1]+bins[1:])/2

time = []
etime = []
for i in range(bins.shape[0]):
    time.append(dt.datetime.fromtimestamp(bins[i]))
    etime.append(dt.timedelta(seconds=ebins[i]))
plt.errorbar(time, hist, xerr=etime, yerr=ehist, capsize=0, linestyle="None", ecolor="black")
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
plt.setp(plt.gca().get_xticklabels(), rotation=10, fontsize=10)
plt.xlabel("Time")
plt.ylabel("Count Rate (counts/s)")
plt.show()
