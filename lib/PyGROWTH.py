#!/usr/bin/env python
# -*- coding: utf-8 -*-
from scipy.optimize import curve_fit
import astropy.io.fits as fitsio
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import warnings
import yaml
warnings.simplefilter("ignore", UserWarning)

class pipeline:
    def __init__(self,input_file,credit=True):
        try:
            if credit:
                self.credit()
            self.f = fitsio.open(input_file)
            self.error = False
            if len(self.f)<2:
                self.error = True
                print("Error: {} is broken.".format(input_file))
            event = self.f[1].data
            self.det_id = self.f[1].header["DET_ID"]
            self.file_date = self.f[1].header["FILEDATE"]
            self.initial_selection(event)
            self.calc_duration()
            self.pipeline_ver = "PyGROWTH Ver.1.0"
        except:
            self.error = True
            print("Error: {} is broken or not a FITS file.".format(input_file))
            
    def close(self):
        self.f.close()

    def info(self):
        self.f.info()
        
    def credit(self):
        print("")
        print("  ------------------------------------------------------")
        print("    PyGROWTH: Pipeline and analysis tool for GROWTH")
        print("    Version 1.0 (April 11, 2024")
        print("    Developed by Yuuki Wada (Osaka University)")
        print("    https://github.com/YuukiWada/PyGROWTH")
        print("")
        print("    The author is not responsible for") 
        print("    any disadvantages and damages caused by this tool.")
        print("  ------------------------------------------------------")
        print("")

    def initial_selection(self,event):
        adc_channel = np.array(event["boardIndexAndChannel"])
        time_tag = np.array(event["timeTag"])
        trigger_count = np.array(event["triggerCount"],dtype=np.uint16)
        pha_max = np.array(event["phaMax"])
        pha_max_time = np.array(event["phaMaxTime"])
        pha_min = np.array(event["phaMin"])
        pha_first = np.array(event["phaFirst"])
        pha_last = np.array(event["phaLast"])
        max_derivative = np.array(event["maxDerivative"])
        mask = ~(pha_max>10000)
        self.adc_channel = adc_channel[mask]
        self.time_tag = time_tag[mask]
        self.trigger_count = trigger_count[mask]
        self.pha_max = pha_max[mask]
        self.pha_max_time = pha_max_time[mask]
        self.pha_min = pha_min[mask]
        self.pha_first = pha_first[mask]
        self.pha_last = pha_last[mask]
        self.max_derivative = max_derivative[mask]

    def calc_duration(self):
        time_tag = self.time_tag
        time_tag[time_tag-time_tag[0]<-1*2**20]+=2**40
        self.duration = (time_tag[-1]-time_tag[0])*1e-8
                
    def ql_spec(self,output_dir="."):
        for ch in range(4):
            output_file = "{}/{}_ch{}_line.pdf".format(output_dir,self.file_date,ch)
            pha = self.pha_max[self.adc_channel==ch] - self.pha_min[self.adc_channel==ch]
            hist, bins = np.histogram(pha,bins=2048,range=(-0.5,2047.5))
            bins = (bins[0:-1]+bins[1:])/2
            plt.step(bins, hist, where="mid")
            plt.xlabel("Channel")
            plt.ylabel("Counts/bin")
            plt.yscale("log")
            plt.xlim(-0.5, 500.5)
            plt.savefig(output_file)
            plt.clf()

    def energy_cal(self,p0,p1):
        energy = np.zeros(self.pha_max.shape)
        for ch in range(4):
            energy[self.adc_channel==ch] = p0[ch]+(self.pha_max[self.adc_channel==ch] - self.pha_min[self.adc_channel==ch])*p1[ch]
        self.energy = energy

        
    def peak_fit(self,ch,pre_peak,pre_sigma,pre_norm,plot=False,plot_dir="."):
        def gaussian(x,a,b,c):
            return a*np.exp(-1*(x-b)**2/(2*c**2))
        center_time = dt.datetime.strptime(self.file_date,"%Y%m%d_%H%M%S")+dt.timedelta(seconds=self.duration/2)
        fit_range_const = [0.6, 1.4]
        hard_limit = [0.85, 1.15] 
        if (self.duration>600) and (self.duration<2400):
            pha = self.pha_max[self.adc_channel==ch] - self.pha_min[self.adc_channel==ch]
            hist, bins = np.histogram(pha,bins=2048,range=(-0.5,2047.5))
            bins = (bins[0:-1]+bins[1:])/2
            peak = [0,0]
            sigma = [0,0]
            norm = [0,0]
            if plot:
                plt.step(bins, hist, where="mid")
            
            for i in range(2):
                # First trial
                fit_range = [pre_peak[i]-pre_sigma[i]*fit_range_const[0],pre_peak[i]+pre_sigma[i]*fit_range_const[1]]
                hist_sel = hist[(bins>=fit_range[0])&(bins<=fit_range[1])]
                bins_sel = bins[(bins>=fit_range[0])&(bins<=fit_range[1])]
                popt,pocv = curve_fit(gaussian, bins_sel, hist_sel, p0=[pre_norm[i],pre_peak[i],pre_sigma[i]], maxfev=50000)
                if plot:
                    plt.plot(bins_sel,gaussian(bins_sel,popt[0],popt[1],popt[2]))

                # Second trial
                fit_range = [popt[1]-popt[2]*fit_range_const[0],popt[1]+popt[2]*fit_range_const[1]]
                hist_sel = hist[(bins>=fit_range[0])&(bins<=fit_range[1])]
                bins_sel = bins[(bins>=fit_range[0])&(bins<=fit_range[1])]
                popt,pocv = curve_fit(gaussian, bins_sel, hist_sel, p0=popt, maxfev=50000)
                if plot:
                    plt.plot(bins_sel,gaussian(bins_sel,popt[0],popt[1],popt[2]))
                norm[i] = popt[0]
                peak[i] = popt[1]
                sigma[i] = popt[2]

            if plot:
                os.makedirs(plot_dir, exist_ok=True)
                output_file = "{}/{}_ch{}_fitting.pdf".format(plot_dir,self.file_date,ch)
                plt.xlim(-0.5, 500.5)
                plt.yscale("log")
                plt.savefig(output_file)

            if (peak[0]<pre_peak[0]*hard_limit[0]) or (peak[0]>pre_peak[0]*hard_limit[1]) or (peak[1]<pre_peak[1]*hard_limit[0]) or (peak[1]>pre_peak[1]*hard_limit[1]):
                prefit = False
            else:
                prefit = True
            result = True
        else:
            result = False
            prefit = False
            peak = [0,0]
            sigma = [0,0]
            norm = [0,0]
        return result, prefit, center_time, peak, sigma, norm

    def timing_cal(self):
        self.gps_verification()
        if self.gps_flag:
            print("GPS available: {}".format(self.file_date))
            gps_timetag = self.gps_timetag
            time_tag = self.time_tag
            unixtime = np.zeros(time_tag.shape,np.float128)
            precise_time = np.zeros(time_tag.shape)
            time_tag[time_tag-time_tag[0]<-1*2**20]+=2**40
            gps_unixtime_start = self.gps_unixtime_start            
            gps_timetag = gps_timetag&0xFFFFFFFFFF
            gps_timetag[gps_timetag<gps_timetag[0]]+=2**40
            gps_unixtime = np.round((gps_timetag-gps_timetag[0])*1e-8)+gps_unixtime_start
            if gps_timetag.shape[0]==3:
                gradient = (gps_unixtime[1]-gps_unixtime[0])/(gps_timetag[1]-gps_timetag[0])
                unixtime[time_tag<=gps_timetag[1]] = gps_unixtime[0] + gradient*(time_tag[time_tag<=gps_timetag[1]]-gps_timetag[0])
                gradient = (gps_unixtime[2]-gps_unixtime[1])/(gps_timetag[2]-gps_timetag[1])
                unixtime[time_tag>gps_timetag[1]] = gps_unixtime[1] + gradient*(time_tag[time_tag>gps_timetag[1]]-gps_timetag[1])
            else:
                gradient = (gps_unixtime[1]-gps_unixtime[0])/(gps_timetag[1]-gps_timetag[0])
                mask = (time_tag<=gps_timetag[1])
                unixtime[mask] = gps_unixtime[0] + gradient*(time_tag[mask]-gps_timetag[0])
                for i in range(gps_timetag.shape[0]-3):
                    gradient = (gps_unixtime[i+2]-gps_unixtime[i+1])/(gps_timetag[i+2]-gps_timetag[i+1])
                    mask = (time_tag>gps_timetag[i+1])&(time_tag<=gps_timetag[i+2])
                    unixtime[mask] = gps_unixtime[i+1]+ gradient*(time_tag[mask]-gps_timetag[i+1])
                gradient = (gps_unixtime[-1]-gps_unixtime[-2])/(gps_timetag[-1]-gps_timetag[-2])
                mask = (time_tag>gps_timetag[-2])
                unixtime[mask] = gps_unixtime[-2] + gradient*(time_tag[mask]-gps_timetag[-2])
            self.unixtime = unixtime.astype(np.float64)
            self.precise_time = (unixtime-np.floor(unixtime)).astype(np.float64)
        else:
            print("GPS unavailable: {}".format(self.file_date))
            time_tag = self.time_tag
            time_tag[time_tag-time_tag[0]<-1*2**20]+=2**40
            self.unixtime = dt.datetime.strptime(self.file_date,"%Y%m%d_%H%M%S").timestamp()+(time_tag-time_tag[0])*1e-8
            self.precise_time = (time_tag-time_tag[0])*1e-8-np.floor((time_tag-time_tag[0])*1e-8)

    def gps_verification(self):
        self.gps_flag = False
        if len(self.f)!=3:
            return 0
        gps = self.f[2].data
        gps_timetag = np.array(gps["fpgaTimeTag"],dtype=np.uint64)
        gps_unixtime = np.array(gps["unixTime"],dtype=np.uint64)
        gps_string = gps["gpsTime"]
        gps_string_select = []
        gps_select = np.ones(gps_timetag.shape, dtype=np.bool)*True
        gps_string_select.append(gps_string[0])
        if (gps_string[0].split("\s")[0]=="GP") or (gps_string[0].split("\s")[0]==""):
            return 0
        for i in range(gps_timetag.shape[0]-1):
            if gps_string[i+1].split("\s")[0]=="":
                gps_select[i+1] = False
            elif gps_string[i+1]==gps_string[i]:
                gps_select[i+1] = False
            else:
                gps_string_select.append(gps_string[i+1])
        gps_timetag = gps_timetag[gps_select]
        gps_unixtime = gps_unixtime[gps_select]
        if gps_timetag.shape[0]<3:
            return 0
        time_string = gps_string_select[0][8:14]
        start_time = dt.datetime.fromtimestamp(gps_unixtime[0])
        for i in range(1200):
            time_calc = start_time+dt.timedelta(seconds=(i-600))-dt.timedelta(hours=9)
            if time_calc.strftime("%H%M%S")==time_string:
                self.gps_flag = True
                self.gps_unixtime_start = time_calc.timestamp()
                break
        self.gps_timetag = gps_timetag
        self.gps_unixtime = gps_unixtime
        self.gps_string = gps_string_select

