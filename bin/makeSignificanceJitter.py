#!/usr/bin/env python
import astropy.io.fits as fitsio
from astropy.table import vstack
import numpy as np
import math
import sys
import os

# reading argiment parameters
if len(sys.argv)<6:
    print("Usage: calcSignificance.py <file list> <output file> <adc channel> <bin width (sec)> <threshold (sigma)>")
    exit()

output_file = sys.argv[2]
adc_channel = int(sys.argv[3])
bin_width = int(sys.argv[4])
scan_threshold = float(sys.argv[5])
energy_threshold = 3.0 # MeV
scan_num = 10
scan_int = bin_width/scan_num

# reading file files
with open(output_file, mode="w") as output:
    with open("./error_message.log", mode="a") as e_output:
        with open(sys.argv[1], "r") as f:
            file_list = f.read().split("\n")
        for input_file in file_list:
            if os.path.exists(input_file):
                try:
                    file_name = os.path.basename(input_file).split('.', 1)[0]
                    fits_file = fitsio.open(input_file)
                    event = fits_file[1].data
                    duration = event["unixTime"][-1]-event["unixTime"][0]
                    event = event[event["boardIndexAndChannel"]==adc_channel]
                    energy = (event["energy"]+(np.random.default_rng().random(event["energy"].shape[0])-0.5)*fits_file[1].header["BINW_CH{}".format(adc_channel)])/1000
                    event = event[energy>energy_threshold]
                    time = event["unixTime"]-event["unixTime"][0]

                    max_significance = 0.0
                    bin_num = int(math.ceil(duration/bin_width))

                    for i in range(scan_num):
                        bin_min = i*scan_int
                        bin_max = bin_num*bin_width+bin_min
                        lc, bins = np.histogram(time, bins=bin_num, range=(bin_min, bin_max))
                        std = np.std(lc)
                        ave = np.mean(lc)
                        significance = (lc-ave)/std

                        if np.amax(significance)>3.0:
                            th_select = ave+3.0*std
                            lc_select = lc[(lc<th_select)]
                            std_sel = np.std(lc_select)
                            ave_sel = np.mean(lc_select)
                            sig_sel = (lc-ave_sel)/std_sel
                            if np.amax(sig_sel)>max_significance:
                                max_significance = np.amax(sig_sel)
                                
                    output.write("{}, {}\n".format(file_name,round(max_significance,3)))
                    if max_significance>scan_threshold:
                        print("{}, {}".format(file_name,round(max_significance,3)))
                except:
                    print("Error: {}".format(input_file))
                    e_output.write("Error: {}\n".format(input_file))
