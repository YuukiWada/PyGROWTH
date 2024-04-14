#!/usr/bin/env python
import astropy.io.fits as fitsio
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import subprocess
import glob
import math
import sys
import os

if len(sys.argv)<7:
    print("Usage: makeQLlightcurveDays.py <input directory> <output directory> <start date> <end date> <adc channel> <threshold>")
    exit()

input_dir = sys.argv[1]
output_dir = sys.argv[2]
date = [sys.argv[3], sys.argv[4]]
adc_channel = int(sys.argv[5])
threshold = int(sys.argv[6])
clock = 1.0e8

date_start = dt.datetime.strptime(date[0], "%Y%m%d")
date_end = dt.datetime.strptime(date[1], "%Y%m%d")
num_days = (date_end-date_start).days+1

# create output directories
#bin_width = [0.01, 1.0, 10.0]
bin_width = [1.0, 10.0]
#dir_name = ["10ms", "1s", "10s"]
dir_name = ["1s", "10s"]
low_th = [threshold, 30]

for day in range(num_days):
    date_obj = date_start+dt.timedelta(days=day)
    date_calc = date_obj.strftime("%Y%m%d")
    print(date_calc)
    input_file_names = sorted(glob.glob("{}/**/{}*.fits.gz".format(input_dir,date_calc), recursive=True))
    
    if len(input_file_names)!=0:
        duration = 0
        for i, input_file in enumerate(input_file_names):
            #print(input_file)
            try:
                fits_file = fitsio.open(input_file)
                event = fits_file[1].data
                if event["timeTag"][0]>event["timeTag"][-1]:
                    obs_time = (event["timeTag"][-1]-event["timeTag"][0]+2**40)/clock
                else:
                    obs_time = (event["timeTag"][-1]-event["timeTag"][0])/clock
                event_high = event[(event["boardIndexAndChannel"]==adc_channel)&(event["phaMax"]-event["phaMin"]>low_th[0])]
                event_all = event[(event["boardIndexAndChannel"]==adc_channel)&(event["phaMax"]-event["phaMin"]>low_th[1])]
                time_tag = event_all["timeTag"]
                time_tag[time_tag[0]-time_tag>2**12]+=2**40
                time_tag_high = event_high["timeTag"]
                time_tag_high[time_tag_high[0]-time_tag_high>2**12]+=2**40

                time_calc_all = duration+(time_tag-time_tag[0])/(clock*3600)
                time_calc_high = duration+(time_tag_high-time_tag_high[0])/(clock*3600)
                if i==0:
                    time_all = time_calc_all
                    time_high = time_calc_high
                else:
                    time_all = np.concatenate([time_all,time_calc_all])
                    time_high = np.concatenate([time_high,time_calc_high])
                duration+=obs_time/3600.0
            except:
                print("The FITS file is corrupted.")

        for i in range(2):
            os.makedirs("./{}/high_{}".format(output_dir,dir_name[i]), exist_ok=True)  
            bin_max = 25.0
            bin_num = round(bin_max*3600.0/bin_width[i])
            output_file = "./{}/high_{}/{}.png".format(output_dir,dir_name[i],date_calc)
            hist, bins = np.histogram(time_high,bins=bin_num,range=(0,bin_max))
            ehist = hist**0.5/(bin_width[i])
            hist = hist/(bin_width[i])
            ebins = (bins[1:]-bins[0:-1])/2
            bins = (bins[0:-1]+bins[1:])/2
            #plt.errorbar(bins, hist, xerr=ebins, yerr=ehist, capsize=0, linestyle="None", ecolor="black")
            plt.plot(bins, hist, color="black")
            plt.xlabel("Hours")
            plt.ylabel("Count Rate (counts/sec)")
            plt.savefig(output_file)
            plt.clf()

            os.makedirs("./{}/all_{}".format(output_dir,dir_name[i]), exist_ok=True)  
            bin_max = 25.0
            bin_num = round(bin_max*3600.0/bin_width[i])
            output_file = "./{}/all_{}/{}.png".format(output_dir,dir_name[i],date_calc)
            hist, bins = np.histogram(time_all,bins=bin_num,range=(0,bin_max))
            ehist = hist**0.5/(bin_width[i])
            hist = hist/(bin_width[i])
            ebins = (bins[1:]-bins[0:-1])/2
            bins = (bins[0:-1]+bins[1:])/2
            #plt.errorbar(bins, hist, xerr=ebins, yerr=ehist, capsize=0, linestyle="None", ecolor="black")
            plt.plot(bins, hist, color="black")
            plt.xlabel("Hours")
            plt.ylabel("Count Rate (counts/sec)")
            plt.savefig(output_file)
            plt.clf()
