#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyGROWTH import pipeline
from matplotlib import dates as mdates
import matplotlib.pyplot as plt
import astropy.io.fits as fitsio
import datetime as dt
import numpy as np
import subprocess
import shutil
import json
import yaml
import glob
import sys
import os

current_dir = os.getcwd()
process_name = current_dir.split("/")[-1]
fiscal_year = current_dir.split("/")[-2]
det_id = current_dir.split("/")[-3]

process_type = sys.argv[1]

if process_type=="ql":
    if len(sys.argv)==3:
        f = open(sys.argv[2], "r")
        file_list = f.read().split("\n")
        f.close()
    elif len(sys.argv)==2:
        data_dir = "{}_fits_lv0".format(process_name)
        file_list = sorted(glob.glob("{}/*.fits*".format(data_dir)))
    else:
        print("Usage: growth-pipeline.py ql [list file]")
        print("Usage: growth-pipeline.py process <configuration file> [list file]")
        exit()
    work_dir = "{}_work".format(process_name)
    os.makedirs(work_dir, exist_ok=True)
    products_dir = "{}_work/products/".format(process_name)
    os.makedirs(products_dir, exist_ok=True) 
    fits = pipeline(file_list[0])
    fits.ql_spec(products_dir)
    fits.close()
        
elif process_type=="process":
    detector_par = ["UNDEFINED","UNDEFINED","UNDEFINED","UNDEFINED","UNDEFINED","UNDEFINED","UNDEFINED"]
    caldb_dir = subprocess.run("echo $GROWTHCALDB",capture_output=True,text=True,shell=True).stdout.strip()
    caldb_file = "{}/{}_winter/{}.json".format(caldb_dir,fiscal_year.upper(),det_id)
    if os.path.exists(caldb_file):
        with open(caldb_file, "r") as json_file:
            json_load = json.load(json_file)
            detector_par[0] = json_load["detectorInfo"][1]["OBS_SITE"]
            detector_par[1] = json_load["detectorInfo"][2]["DET_CH0"]
            detector_par[2] = json_load["detectorInfo"][3]["DET_CH1"]
            detector_par[3] = json_load["detectorInfo"][4]["DET_CH2"]
            detector_par[4] = json_load["detectorInfo"][5]["DET_CH3"]
            detector_par[5] = json_load["detectorInfo"][6]["INS_DATE"]
            detector_par[6] = json_load["detectorInfo"][7]["RMV_DATE"]
    if len(sys.argv)==4:
        f = open(sys.argv[3], "r")
        file_list = f.read().split("\n")
        f.close()
    elif len(sys.argv)==3:
        data_dir = "{}_fits_lv0".format(process_name)
        file_list = sorted(glob.glob("{}/*.fits*".format(data_dir)))
    else:
        print("Usage: growth-pipeline.py ql [list file]")
        print("Usage: growth-pipeline.py process <configuration file> [list file]")
        exit()
    lv2_dir = "{}_fits_lv2".format(process_name)
    os.makedirs(lv2_dir, exist_ok=True)

    yaml_file = open(sys.argv[2], "r")
    config = yaml.safe_load(yaml_file)

    peak_energy = [1460.8, 2614.5]
    prefit = np.zeros((4,3,2))
    for i in range(4):
        if config["ch{}".format(i)][0]=="track":
            prefit[i,0,0] = float(config["ch{}".format(i)][1])*peak_energy[0]/(float(config["ch{}".format(i)][2])*1000)
            prefit[i,0,1] = float(config["ch{}".format(i)][1])*peak_energy[1]/(float(config["ch{}".format(i)][2])*1000)
            if config["ch{}".format(i)][3]=="bgo":
                resolution = 0.15
                prefit[i,1,0] = prefit[i,0,0]*resolution/2.35 
                prefit[i,1,1] = prefit[i,0,1]*resolution/2.35 
            else:
                resolution = 0.1
                prefit[i,1,0] = prefit[i,0,0]*resolution/2.35 
                prefit[i,1,1] = prefit[i,0,1]*resolution/2.35 
            prefit[i,2,0] = 1000
            prefit[i,2,1] = 1000
        elif config["ch{}".format(i)][0]=="const":
            prefit[i,0,0] = float(config["ch{}".format(i)][1])*peak_energy[0]/(float(config["ch{}".format(i)][2])*1000)
            prefit[i,0,1] = float(config["ch{}".format(i)][1])*peak_energy[1]/(float(config["ch{}".format(i)][2])*1000)

    result_peak_k = [[],[],[],[]]
    result_peak_tl = [[],[],[],[]]
    result_time = [[],[],[],[]]
    for input_file in file_list:
        if not os.path.exists(input_file):
            print("{} does not exist.".format(input_file))
            continue
        fits = pipeline(input_file,credit=False)
        if fits.error:
            error_dir = "{}_fits_error".format(process_name)
            os.makedirs(error_dir, exist_ok=True)
            shutil.copy2(input_file,error_dir)
            continue
        if det_id!=fits.det_id:
            print("{}: Detector ID has inconsistency.".format(input_file))
            continue

        # Energy calibration
        switch = [False,False,False,False]
        p0 = [0,0,0,0]
        p1 = [0,0,0,0]
        for ch in range(4):
            if config["ch{}".format(ch)][0]=="track":
                result, prefit_flag, center_time, peak, sigma, norm = fits.peak_fit(ch,prefit[ch,0,:],prefit[ch,1,:],prefit[ch,2,:],plot=False,plot_dir="{}_work/fit_results".format(process_name))
                if result:
                    p1[ch] = (peak_energy[1]-peak_energy[0])/(peak[1]-peak[0])
                    p0[ch] = peak_energy[0]-peak[0]*p1[ch]
                    result_peak_k[ch].append(peak[0])
                    result_peak_tl[ch].append(peak[1])
                    result_time[ch].append(center_time)
                    switch[ch] = True
                if prefit_flag:
                    prefit[ch,0,:] = np.array(peak)
                    prefit[ch,1,:] = np.array(sigma)
                    prefit[ch,2,:] = np.array(norm)
            elif config["ch{}".format(ch)][0]=="const":
                    p1[ch] = (peak_energy[1]-peak_energy[0])/(prefit[ch,0,1]-prefit[ch,0,0])
                    p0[ch] = peak_energy[0]-prefit[ch,0,0]*p1[ch]
                    switch[ch] = True
        fits.energy_cal(p0,p1)

        # Timing calibration
        fits.timing_cal()

        # Write fits file
        output_file = "{}/{}".format(lv2_dir,os.path.basename(input_file))
        primary_hdu = fits.f[0]
        ecol1 = fitsio.Column(name="boadIndexAndChannel", format="B", array=fits.adc_channel)
        ecol2 = fitsio.Column(name="timeTag", format="K", array=fits.time_tag)
        ecol3 = fitsio.Column(name="triggerCount", format="J", array=fits.trigger_count)
        ecol4 = fitsio.Column(name="phaMax", format="I", array=fits.pha_max)
        ecol5 = fitsio.Column(name="phaMaxTime", format="I", array=fits.pha_max_time)
        ecol6 = fitsio.Column(name="phaMin", format="I", array=fits.pha_min)
        ecol7 = fitsio.Column(name="phaFirst", format="I", array=fits.pha_first)
        ecol8 = fitsio.Column(name="phaLast", format="I", array=fits.pha_last)
        ecol9 = fitsio.Column(name="maxDerivative", format="I", array=fits.max_derivative)
        ecol10 = fitsio.Column(name="unixTime", format="D", unit="sec", array=fits.unixtime)
        ecol11 = fitsio.Column(name="preciseTime", format="E", unit="sec", array=fits.precise_time)
        ecol12 = fitsio.Column(name="energy", format="E", unit="keV", array=fits.energy)
        ecol  = fitsio.ColDefs([ecol1, ecol2, ecol3, ecol4, ecol5, ecol6, ecol7, ecol8, ecol9, ecol10, ecol11, ecol12])
        event_hdu = fitsio.BinTableHDU.from_columns(ecol,header=fits.f[1].header,name="EVENTS")
        event_hdu.header.append(("PIPELINE", "level-2", "present pipeline process level"))
        event_hdu.header.append(("PL1_DATE", dt.datetime.now().strftime("%Y%m%d_%H%M%S"), "pipeline level-1 processing date"))
        event_hdu.header.append(("PL1_VER" , fits.pipeline_ver, "pipeline level-1 version"))
        event_hdu.header.append(("OBS_SITE", detector_par[0], "observation site"))
        event_hdu.header.append(("DET_CH0" , detector_par[1], "scintillator ID of Channel 0"))
        event_hdu.header.append(("DET_CH1" , detector_par[2], "scintillator ID of Channel 1"))
        event_hdu.header.append(("DET_CH2" , detector_par[3], "scintillator ID of Channel 2"))
        event_hdu.header.append(("DET_CH3" , detector_par[4], "scintillator ID of Channel 3"))
        event_hdu.header.append(("INS_DATE", detector_par[5], "installation date"))
        event_hdu.header.append(("RMV_DATE", detector_par[6], "removal date"))
        event_hdu.header.append(("PL2_DATE", dt.datetime.now().strftime("%Y%m%d_%H%M%S"), "pipeline level-2 processing date"))
        event_hdu.header.append(("PL2_VER" , fits.pipeline_ver, "pipeline level-2 version"))
        if fits.gps_flag:
            event_hdu.header.append(("TIME_CAL", "2", "method of time calibration"))
        else:
            event_hdu.header.append(("TIME_CAL", "1", "method of time calibration"))
        event_hdu.header.add_comment("TIME_CAL: method of time calibration 0:from file name, 1:from UNIXTIME, 2:from GPS time")
        for ch in range(4):
            if switch[ch]:
                event_hdu.header.append(("CAL_CH{}".format(ch) , "YES", "Channel {} are calibrated or not.".format(ch)))
                event_hdu.header.append(("LNK_CH{}".format(ch) , round(result_peak_k[ch][-1],2), "40K peak mean in Channel {} (adc channel)".format(ch)))
                event_hdu.header.append(("LNTL_CH{}".format(ch), round(result_peak_tl[ch][-1],2), "208Tl peak mean in Channel {} (adc channel)".format(ch)))
                event_hdu.header.append(("BINW_CH{}".format(ch), round(p1[ch],4), "Bin width of Channel {} in Energy (keV)".format(ch)))
            else:
                event_hdu.header.append(("CAL_CH{}".format(ch) , "NO", "Channel {} are calibrated or not.".format(ch)))
                event_hdu.header.append(("LNK_CH{}".format(ch) , "NONE", "40K peak mean in Channel {} (adc channel)".format(ch)))
                event_hdu.header.append(("LNTL_CH{}".format(ch), "NONE", "208Tl peak mean in Channel {} (adc channel)".format(ch)))
                event_hdu.header.append(("BINW_CH{}".format(ch), "NONE", "Bin width of Channel {} in Energy (keV)".format(ch)))
        gcol1 = fitsio.Column(name="fpgaTimeTag", format="K", array=fits.gps_timetag)
        gcol2 = fitsio.Column(name="unixTime", format="J", array=fits.gps_unixtime)
        gcol3 = fitsio.Column(name="gpsTime", format="14A", array=fits.gps_string)
        gcol  = fitsio.ColDefs([gcol1, gcol2, gcol3])
        gps_hdu = fitsio.BinTableHDU.from_columns(gcol,header=fits.f[2].header,name="GPS")
        hdu_list = fitsio.HDUList(hdus=[primary_hdu,event_hdu,gps_hdu])
        hdu_list.writeto(output_file,overwrite=True)

    # Plot peak variation
    for ch in range(4):
        if config["ch{}".format(ch)][0]=="track":
            output_file = "{}_work/peak_transition_ch{}.pdf".format(process_name,ch)
            plt.plot(result_time[ch],result_peak_k[ch],color="black",label="40K")
            plt.plot(result_time[ch],result_peak_tl[ch],color="red",label="208Tl")
            plt.legend()
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
            plt.setp(plt.gca().get_xticklabels(), rotation=10, fontsize=10)
            plt.xlabel("Date")
            plt.ylabel("Channel")
            plt.savefig(output_file)
            plt.clf()

else:
    print("Error: Mode should be 'ql' or 'process'.")
    exit()
