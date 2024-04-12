#!/usr/bin/env python
import astropy.io.fits as fitsio
import matplotlib.pyplot as plt
import numpy as np
import math
import signal
import sys
signal.signal(signal.SIGINT, signal.SIG_DFL)

# reading argiment parameters
if len(sys.argv)<5:
    print("Usage: makeQLscatter.py <input file> <adc channel> <center time (sec)> <time width>")
    exit()

input_file = sys.argv[1]
adc_channel = int(sys.argv[2])
center_sec = float(sys.argv[3])
width = float(sys.argv[4])

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
time_tag = event["timeTag"]
time_tag[time_tag-time_tag[0]<-1*2**20]+=2**40
seconds = (time_tag-time_tag[0])/clock
print("Duration: {} sec".format(round(duration,1)))

mask = (seconds>=center_sec-width/2)&(seconds<=center_sec+width/2)
event = event[mask]
seconds = seconds[mask]

# plot histogram
plt.scatter(seconds, event["phaMax"], color="black", s=5)
plt.scatter(seconds, event["phaMin"], color="red", s=5)
plt.xlabel("Seconds")
plt.ylabel("ADC channel")
plt.xlim(center_sec-width/2,center_sec+width/2)
plt.show()
