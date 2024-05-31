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
dir_name = ["1s", "10s"]
low_th = [threshold, 30]

for day in range(num_days):
    date_obj = date_start+dt.timedelta(days=day)
    date_calc = date_obj.strftime("%Y%m%d")
    unixtime_range = dt.datetime.timestamp(date_obj)
    print(date_calc)
    date_calc_before = (date_obj-dt.timedelta(days=1)).strftime("%Y%m%d")
    date_calc_after = (date_obj+dt.timedelta(days=1)).strftime("%Y%m%d")

    input_file_names = []
    input_files = sorted(glob.glob("{}/**/{}*.fits.gz".format(input_dir,date_calc_before), recursive=True))
    if len(input_files)>0:
        input_file_names += [input_files[-1]]
    input_files = sorted(glob.glob("{}/**/{}*.fits.gz".format(input_dir,date_calc), recursive=True))
    if len(input_files)>0:
        input_file_names += input_files
    input_files = sorted(glob.glob("{}/**/{}*.fits.gz".format(input_dir,date_calc_after), recursive=True))
    if len(input_files)>0:
        input_file_names += [input_files[0]]
    
    if len(input_file_names)!=0:
        duration = 0
        for i, input_file in enumerate(input_file_names):
            try:
                fits_file = fitsio.open(input_file)
                event = fits_file[1].data
                if event["timeTag"][0]>event["timeTag"][-1]:
                    obs_time = (event["timeTag"][-1]-event["timeTag"][0]+2**40)/clock
                else:
                    obs_time = (event["timeTag"][-1]-event["timeTag"][0])/clock
                if obs_time>2000:
                    print("Observation duration is incorrect: {}".format(input_file))
                    print(fits_file)
                if obs_time<0:
                    print("Observation duration is incorrect: {}".format(input_file))
                    print(fits_file)
                event_high = event[(event["boardIndexAndChannel"]==adc_channel)&(event["phaMax"]-event["phaMin"]>low_th[0])]
                event_all = event[(event["boardIndexAndChannel"]==adc_channel)&(event["phaMax"]-event["phaMin"]>low_th[1])]
                time_tag = event_all["timeTag"]
                time_tag[time_tag[0]-time_tag>2**12]+=2**40
                time_tag_high = event_high["timeTag"]
                time_tag_high[time_tag_high[0]-time_tag_high>2**12]+=2**40

                if i==0:
                    date_file = dt.datetime.strptime(input_file.split("/")[-1].split(".")[0], "%Y%m%d_%H%M%S")
                    unixtime_all = dt.datetime.timestamp(date_file)+(time_tag.astype(np.float128)-event["timeTag"][0])/clock
                    unixtime_high = dt.datetime.timestamp(date_file)+(time_tag_high.astype(np.float128)-event["timeTag"][0])/clock
                    time_tag_last = event["timeTag"][-1]
                    unixtime_last = unixtime_all[-1] 
                else:
                    date_file = dt.datetime.strptime(input_file.split("/")[-1].split(".")[0], "%Y%m%d_%H%M%S")
                    if abs(unixtime_last-dt.datetime.timestamp(date_file))>=60:
                        unixtime_all = np.concatenate([unixtime_all,dt.datetime.timestamp(date_file)+(time_tag.astype(np.float128)-event["timeTag"][0])/clock])
                        unixtime_high = np.concatenate([unixtime_high,dt.datetime.timestamp(date_file)+(time_tag_high.astype(np.float128)-event["timeTag"][0])/clock])
                        time_tag_last = event["timeTag"][-1]
                        unixtime_last = unixtime_all[-1]
                    elif abs(time_tag_last-event["timeTag"][0])<60*clock:
                        unixtime_all = np.concatenate([unixtime_all,unixtime_last+(time_tag.astype(np.float128)-time_tag_last)/clock])
                        unixtime_high = np.concatenate([unixtime_high,unixtime_last+(time_tag_high.astype(np.float128)-time_tag_last)/clock])
                        time_tag_last = event["timeTag"][-1]
                        unixtime_last = unixtime_all[-1]
                    elif abs(time_tag_last-(event["timeTag"][0]+2**40))<60*clock:
                        unixtime_all = np.concatenate([unixtime_all,unixtime_last+(time_tag.astype(np.float128)+2**40-time_tag_last)/clock])
                        unixtime_high = np.concatenate([unixtime_high,unixtime_last+(time_tag_high.astype(np.float128)+2**40-time_tag_last)/clock])
                        time_tag_last = event["timeTag"][-1]
                        unixtime_last = unixtime_all[-1]
                    else:
                        unixtime_all = np.concatenate([unixtime_all,dt.datetime.timestamp(date_file)+(time_tag.astype(np.float128)-event["timeTag"][0])/clock])
                        unixtime_high = np.concatenate([unixtime_high,dt.datetime.timestamp(date_file)+(time_tag_high.astype(np.float128)-event["timeTag"][0])/clock])
                        time_tag_last = event["timeTag"][-1]
                        unixtime_last = unixtime_all[-1]
                        
                time_all  = (unixtime_all-unixtime_range)/3600
                time_high = (unixtime_high-unixtime_range)/3600
                itme_all  = time_all[(time_all>=0)&(time_all<=24)]
                itme_high = time_high[(time_high>=0)&(time_high<=24)]
            except Exception as e:
                print("The FITS file {} is corrupted.".format(input_file))
                print(e)

        for i in range(2):
            os.makedirs("./{}/high_{}".format(output_dir,dir_name[i]), exist_ok=True)  
            bin_max = 24.0
            bin_num = round(bin_max*3600.0/bin_width[i])
            output_file = "./{}/high_{}/{}.png".format(output_dir,dir_name[i],date_calc)
            hist, bins = np.histogram(time_high,bins=bin_num,range=(0,bin_max))
            ehist = hist**0.5/(bin_width[i])
            hist = hist/(bin_width[i])
            hist[hist==0] = np.nan
            ebins = (bins[1:]-bins[0:-1])/2
            bins = (bins[0:-1]+bins[1:])/2
            #plt.errorbar(bins, hist, xerr=ebins, yerr=ehist, capsize=0, linestyle="None", ecolor="black")
            plt.plot(bins, hist, color="black")
            plt.grid(axis="x")
            plt.xlabel("Hours")
            plt.ylabel("Count Rate (counts/sec)")
            plt.xlim(0,bin_max)
            plt.xticks(np.arange(0,bin_max+1,1))
            
            plt.savefig(output_file)
            plt.clf()

            os.makedirs("./{}/all_{}".format(output_dir,dir_name[i]), exist_ok=True)  
            bin_max = 24.0
            bin_num = round(bin_max*3600.0/bin_width[i])
            output_file = "./{}/all_{}/{}.png".format(output_dir,dir_name[i],date_calc)
            hist, bins = np.histogram(time_all,bins=bin_num,range=(0,bin_max))
            ehist = hist**0.5/(bin_width[i])
            hist = hist/(bin_width[i])
            hist[hist==0] = np.nan
            ebins = (bins[1:]-bins[0:-1])/2
            bins = (bins[0:-1]+bins[1:])/2
            #plt.errorbar(bins, hist, xerr=ebins, yerr=ehist, capsize=0, linestyle="None", ecolor="black")
            plt.plot(bins, hist, color="black")
            plt.grid(axis="x")
            plt.xlabel("Hours")
            plt.ylabel("Count Rate (counts/sec)")
            plt.xlim(0,bin_max)
            plt.xticks(np.arange(0,bin_max+1,1))
            plt.savefig(output_file)
            plt.clf()

