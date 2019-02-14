# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 17:53:53 2018

@author: mitadm
"""
import numpy as np
import scipy.optimize as opt

from mpl_toolkits.mplot3d import Axes3D
from numpy import ndarray
from datetime import datetime


import matplotlib.pyplot as pl

import scipy.optimize as optimization

from scipy.optimize import curve_fit

import scipy.signal as sgl
import pylab
from scipy.signal import butter, lfilter, freqz
import os
import re
import math
import scipy
import pickle  
import sklearn as sk

font =11.5
#print pl.rcParams.keys()
pl.rcParams.update({'font.size': font})
#pl.rc('text', usetex=False)
pl.rc('figure' , figsize=(10,5))
pl.figure(figsize=(10,5))
pl.rcParams['axes.linewidth'] = 2
pl.rcParams['axes.grid'] = False
pl.rcParams['axes.grid.which'] = 'both'
pl.rcParams['xtick.minor.visible'] = True
pl.rcParams['xtick.minor.size'] = 4
pl.rcParams['xtick.minor.width'] = 1
pl.rcParams['xtick.major.size'] = 6
pl.rcParams['xtick.major.width'] = 2
pl.rcParams['ytick.minor.visible'] = True
pl.rcParams['ytick.minor.size'] = 4
pl.rcParams['ytick.minor.width'] = 1
pl.rcParams['ytick.major.size'] = 6
pl.rcParams['ytick.major.width'] = 2
pl.rcParams['agg.path.chunksize'] = 10000
#pl.rcParams['axes.formatter.use_locale'] = 'right'
pl.rcParams['legend.numpoints'] = 1
pl.rcParams['legend.frameon'] = False

fe=1000.0
def getWindowedFFT(x, scaling="density", window = "boxcar", oneSided = True):
    '''
    set scaling to "density" for comparison with LPSD
    returns frequencies, FFT
    '''
    win_vals = sgl.get_window(window, len(x))
    windowed_x = win_vals*x
    FTfreq=np.arange(0,len(x))*fe/len(x)
    FT=np.fft.fft(windowed_x)
    if scaling == "density":
        S1 = np.sum(win_vals)
        S2 = np.sum(win_vals**2)
        ENBW = fe*S2/(float(S1**2))

        FT = FT/(S1*np.sqrt(ENBW))
    if oneSided:
        FTfreq=FTfreq[0:int(1+len(FTfreq)/2)]
        FT=FT[0:int(1+len(FT)/2)]
        for j in range(1, 1+len(FT)/2):
            FT[j]=2*FT[j]
            
    return FTfreq, FT

def highPassFilter(timeseries, filtersize=2,  inverse=False):
        b, a =sgl.butter(1,2*filtersize/fe,btype="highpass")
        if inverse==False:
            z=lfilter(b,a,timeseries)
        
        if inverse==True:
            z=lfilter(a,b,timeseries)
        
        return z

def removeDC(timeseries):
    fourier = np.fft.fft(timeseries)
    fourier[0]=0
    return np.real(np.fft.ifft(fourier))

def analyzeFitResults(list_ofPulse_Dict_Lists):
    chis=[]
    peaks=[]
    obs_pulses=[]
    templates=[]
    slow_amps=[]
    fast_amps=[]
    times=[]
    amps=[]
    delays=[]
    errors=[]
    for run in list_ofPulse_Dict_Lists:
        for result in run: #result = a pulse
            chis.append(result.get("chi")    )
            peaks.append(result.get("peak")    )
            obs_pulses.append(result.get("obs pulse") )
            templates.append(result.get("template") )
            slow_amps.append(result.get("slow amplitude") )
            fast_amps.append(result.get("fast amplitude") )
            amps.append(result.get("amplitude") )
            times.append( result.get("time")   )
            delays.append(result.get("delay") )
            errors.append(result.get("amp error") )
     
    total=0
    for i in list_ofPulse_Dict_Lists:
        print "Number of Pulse Fit is ", len(i)
        total += len(i)
    
    if total==0:
        print "NO PULSES FOUND"
        return
    
    DoF=len(getWindowedFFT(templates[0])[1])

    amps=np.array(amps)        
    chis_red=np.array(chis)/DoF
        

        

    print "Number of Total Pulses", total
    one_second = np.linspace(0,1,1000)
    mean_pulse = np.mean(obs_pulses, axis=0) 
    error_on_mean = np.std(obs_pulses, axis=0)/np.sqrt(len(obs_pulses))
    for i in np.random.randint(0, len(obs_pulses)-1, 2):
        pl.plot(one_second, obs_pulses[i], label='Pulse ' + str(i))
    pl.errorbar(one_second, mean_pulse, error_on_mean, label='Average Pulse')
    pl.legend()
    pl.show()

    pl.hist(delays, range=(-0.1,0.37), bins=100)
    pl.xlabel('Time Delay')
    pl.ylabel('Events')
    pl.show() 

#    pl.hist(errors, bins=50)
#    pl.xlabel('Error on Amplitudes')
#    pl.ylabel('Events')
#    pl.show() 
#    
#    pl.scatter(amps, errors)
#    pl.ylabel('Error on Amplitudes')
#    pl.xlabel('Amplitude [A.D.U.]')
#    pl.xlim([0,200])
#    pl.show() 
#
#    pl.scatter(delays, errors)
#    pl.ylabel('Error on Amplitudes')
#    pl.xlabel('Delay of Pulse [s]')
#    pl.show() 


    amp_hist = pl.hist(amps, range=(0,200), bins = 100, color='black')
    pl.xlabel('Amplitude [A.D.U.]')
    pl.ylabel('Events')
    pl.xlim([0,200])
    pl.show()    
    
    amp_hist = pl.hist(amps, range=(0,200), bins = 200, color='black')
    pl.xlabel('Amplitude [A.D.U.]')
    pl.ylabel('Events')
    pl.xlim([0,200])
    pl.show()   
    
    amp_hist2 = pl.hist(amps, range=(0,200), bins = 100, color='black', log=True)
    pl.xlabel('Amplitude [A.D.U.]')
    pl.ylabel('Events')
    pl.xlim([0,200])
    pl.show()    
    
    bin_values_list = list(amp_hist[0])
    peak_1_location = bin_values_list.index(max(bin_values_list))
    print "1st peak location", amp_hist[1][peak_1_location]
    
#    pl.hist(amps, range=(0,200), bins = 100)
#    pl.xlabel('Amplitude [A.D.U.]')
#    pl.ylabel('Events')
#    pl.yscale('log', nonposy='clip')
#    pl.xlim([0,200])
#    pl.show()  
    
#    pl.hist(amps, range=(0,200), bins = 100)
#    pl.xlabel('Amplitude [A.D.U.]')
#    pl.ylabel('Events')
#    pl.yscale('log', nonposy='clip')
#    pl.xscale('log', nonposx='clip')
#    pl.xlim([0,200])
#    pl.show()   
#    
    
    pl.scatter(np.asarray(times)[:,1], chis_red)
    pl.xlabel("time of pulse")
    pl.ylabel("Reduced X-square")
    pl.xlim([0, 3600])
    pl.ylim([0,5])
    pl.show()
#    pl.hist(peaks)
#    pl.xlabel('Peak Height')
#    pl.ylabel('Events')
#    pl.show()
        

    fig = pl.figure()
    ax = pl.gca()
    ax.plot(amps , chis, '.', c='black', markeredgecolor='none')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Chi Sq.')
    ax.minorticks_on()
    x_plot = np.logspace(np.log10(min(amps)), np.log10(max(amps)), 250 )
    if source=="Neutron":
        ax.plot(x_plot, loglog_line(x_plot, n0_PU_points_below), color='blue')
        ax.plot(x_plot, loglog_line(x_plot, n0_PU_points_above), color='green')
        ax.plot(x_plot, loglog_line(x_plot, n0_spike_points_above), color='orange')
        ax.plot(x_plot, cut_parab(x_plot, (a1, c1)), color='red')
    elif source=="Barium":
        ax.plot(x_plot, loglog_line(x_plot, ba_PU_points_above), color='blue')
        ax.plot(x_plot, loglog_line(x_plot, ba_spike_points_below), color='orange')
        ax.plot(x_plot, cut_parab(x_plot, (a1, c1)), color='red')
    elif source=="Fe":
        ax.plot(x_plot, loglog_line(x_plot, fe_PU_points_above), color='blue')
        ax.plot(x_plot, loglog_line(x_plot, fe_spike_points_below), color='orange')
        ax.plot(x_plot, cut_parab(x_plot, (a2, c2)), color='red')
    ax.set_ylim([10**1, 10*np.max(chis)])
    ax.set_xlim([10**0, 3*10**3])

    pl.show()
    
    if slow_amps[0]!=None:
        pl.hist(slow_amps)
        pl.xlabel('Amplitude of Slow Response')
        pl.ylabel('Events')
        pl.show()
        
        pl.hist(fast_amps)
        pl.xlabel('Amplitude of Fast Response')
        pl.ylabel('Events')
        pl.show()
        
        slow_amps=np.asarray(slow_amps)
        fast_amps=np.asarray(fast_amps)
        ratio= fast_amps/slow_amps
        pl.hist(ratio)
        pl.xlabel('Ratio of Fast Amplitude to Slow Amplitude')
        pl.ylabel('Events')
        pl.show() 
    
        pl.hist(1.0/ratio)
        pl.xlabel('Ratio of Slow Amplitude to Fast Amplitude')
        pl.ylabel('Events')
        pl.show() 
        
        pl.scatter(slow_amps, chis)
        pl.xlabel('Slow Amplitude')
        pl.ylabel('Chi Sq')
        pl.show()
        
        pl.scatter(fast_amps, chis)
        pl.xlabel('Fast Amplitude')
        pl.ylabel('Chi Sq')
        pl.show()
        
        pl.scatter(slow_amps, fast_amps, c="black", marker=".")
        pl.xlabel('Slow Amplitude')
        pl.ylabel('Fast Amplitude')
        pl.xlim([0,200])
        pl.ylim([0,500])
        pl.show()
        
        pl.scatter(ratio, chis)
        pl.xlabel('Ratio of Fast Amplitude to Slow Amplitude')
        pl.ylabel('Chi Sq')
        pl.show()
        
        fig = pl.figure()
        ax = Axes3D(fig)
        ax.scatter(slow_amps, fast_amps, chis)
        ax.set_xlabel('Slow Amplitude')
        ax.set_ylabel('Fast Amplitude')
        ax.set_zlabel('Chi Sq')
        pl.show()
    
    
    return
    
def getPulseSubset(pulse_dict_list, cut_list):
    #provide list of cuts where each element is of the form (var, lower_bound, upper_bound)
    subset = []   
    for result in pulse_dict_list:
        veto = False
        for cut in cut_list:
            var = cut[0]
            lower_bound = cut[1]
            upper_bound = cut[2]
            if result[var]>=upper_bound or result[var] <= lower_bound:
                veto = True
                break #skip event if any of cuts not met
        if not veto:
            subset.append(result)
    return subset
 
def cut_parab(amp, (a, c) ):
      return  a*(amp)**2 + c

def line(amp, (m, b)):
    return m*amp + b
  
def step(x):
    return 1 * (x > 0) 
    
def expo_temp1(t, (Tr,Tf)):
    return (-np.exp(-t/Tr)+np.exp(-t/Tf)) 

def A_theta_exp(t, a, t_0, Tr, Tf):
    offset=[]
    offset.append(t_0)
    z=len(t)*offset
    z=np.array(z)
    m=expo_temp1(t-z, (Tr,Tf))
    return step(t-z)*m*a

def AB_theta_exp(t, t_0, Tr, Tf1, Tf2, A, B):
    shifted_t=np.asarray(t)-t_0
    m=expo_temp1(shifted_t, (Tr, Tf1))
    n=expo_temp1(shifted_t, (Tr, Tf2))
    
    H=step(shifted_t)
    
    return H*(A*m+B*n)

def theta_exp(t, (t_0,Tr,Tf)):
    
    offset=[]
    offset.append(t_0)
    z=len(t)*offset
    z=np.array(z)
    m=expo_temp1(t-z, (Tr,Tf))
        
    return step(t-z)*m

def template_1amp(t, (t_0,)):
    exp=theta_exp(t, (t_0, Tr_fit, Tf_fit))
    return exp

def loglog_line(t, (x1,y1,x2,y2)):
    '''
    given two points on a line in a log-log plot (x1,y1) and (x2,y2)
    returns the power law function through the 2 points
    '''
    m = np.log(y2/float(y1))/( np.log(x2/float(x1)) )
    return y1*(t/float(x1))**m
   
def generalizedCut(pulse_dict_list, x_var, y_var, function, parameters, veto_greater_than_func):
    '''
    str: x_var typically "amplitude"
    str: y_var typically "chi" 
    function should be function of x_var with parameters of the form parameters
    veto_greater_than_func is true if you want to veto points above the function
        and false if you want to veto points below the curve  
    '''    
    subset = []
    for result in pulse_dict_list:
        if veto_greater_than_func:
            if result[y_var] <= function(result[x_var], parameters) :
                subset.append(result)
        else:
            if result[y_var] > function(result[x_var], parameters) :
                subset.append(result)
    return subset
            
            
            
            
def compareFitResults(list_of_pulse_dict_list1, list_of_pulse_dict_list2, list_of_pulse_dict_list3 = [] ):
    chis1=[]
    chis2=[]
    slow_amps=[]
    fast_amps=[]
    amps=[]
    assert len(list_of_pulse_dict_list1)==len(list_of_pulse_dict_list2)
    total_pulses=0
    for chunk in range(len(list_of_pulse_dict_list1)):
        print "Total pulses", total_pulses
        assert len(list_of_pulse_dict_list1[chunk])==len(list_of_pulse_dict_list2[chunk])
        total_pulses += len(list_of_pulse_dict_list1[chunk])
        for k in range(len(list_of_pulse_dict_list1[chunk]) ):
            chis1.append(list_of_pulse_dict_list1[chunk][k].get("chi")    )
            chis2.append(list_of_pulse_dict_list2[chunk][k].get("chi")    )
            slow_amps.append(list_of_pulse_dict_list2[chunk][k].get("slow amplitude") )
            fast_amps.append(list_of_pulse_dict_list2[chunk][k].get("fast amplitude") )
            amps.append(list_of_pulse_dict_list1[chunk][k].get("amplitude") )
        
    print "Number of Pulse Fit is ", total_pulses
    chis1 = np.asarray(chis1)
    chis2 = np.asarray(chis2)
    deltaChi = chis1 - chis2
    pl.hist(deltaChi, range=(-1000,1000),bins=100)
    pl.xlabel('Delta Chi Sq')
    pl.ylabel('Events')
    pl.show()

    pl.scatter(amps, deltaChi)
    pl.xlabel('Amplitude')
    pl.xlim([0,200])
    pl.ylim([-2000,2000])
    pl.ylabel('Delta Chi Sq')
    pl.show()
    
    fig = pl.figure()
    ax = pl.gca()
    ax.plot(amps , deltaChi, '.', c='black', markeredgecolor='none')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Delta Chi Sq.')
#    ax.set_ylim([10**2, 10**8])
    pl.show()
    
    fig = pl.figure()
    ax = pl.gca()
    ax.plot(chis1 , chis2, '.', c='black', markeredgecolor='none')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Chi Sq. of One Amp. Fit')
    ax.set_ylabel('Chi Sq. of Two Amp. Fit')
#    ax.set_ylim([10**2, 10**8])
    pl.show()
    
    if len(list_of_pulse_dict_list3) > 0:    
        signal_chis_oneAmp = []
        signal_amp_oneAmp = []
        signal_chis_twoAmp = []
        signal_slowAmp_twoAmp = []
        signal_fastAmp_twoAmp = []
        signal_chis_spike = []
        signal_amp_spike = []
        
        spike_chis_oneAmp = []
        spike_amp_oneAmp = []
        spike_chis_twoAmp = []
        spike_slowAmp_twoAmp = []
        spike_fastAmp_twoAmp = []
        spike_chis_spike = []
        spike_amp_spike = []
        
        PU_chis_oneAmp = []
        PU_amp_oneAmp = []
        PU_chis_twoAmp = []
        PU_slowAmp_twoAmp = []
        PU_fastAmp_twoAmp = []
        PU_chis_spike = []
        PU_amp_spike = []
        
        misc_chis_oneAmp = []
        misc_amp_oneAmp = []
        misc_chis_twoAmp = []
        misc_slowAmp_twoAmp = []
        misc_fastAmp_twoAmp = []
        misc_chis_spike = []
        misc_amp_spike = []
        
        signalBand_twoamp_pulse_dict_list = []
        
        assert len(list_of_pulse_dict_list3) == len(list_of_pulse_dict_list2)
        for k in range(len(list_of_pulse_dict_list1)): #for each run
            for p in range(len(list_of_pulse_dict_list1[k])): #for each pulse
                pulse_dict_oneAmp = list_of_pulse_dict_list1[k][p]
                pulse_dict_twoAmp = list_of_pulse_dict_list2[k][p]
                pulse_dict_spike = list_of_pulse_dict_list3[k][p]
                
                if source == "Neutron":
                    if (pulse_dict_oneAmp["chi"] < cut_parab(pulse_dict_oneAmp["amplitude"], (a1,c1)) 
                        and pulse_dict_oneAmp["amplitude"] > 0):
                    #signal band
                        signal_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        signal_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        signal_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        signal_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        signal_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        signal_chis_spike.append(pulse_dict_spike["chi"])
                        signal_amp_spike.append(pulse_dict_spike["amplitude"])
                        
                        signalBand_twoamp_pulse_dict_list.append(pulse_dict_twoAmp)
                        
                    elif (pulse_dict_oneAmp["chi"] > 2000 
                          and pulse_dict_oneAmp["amplitude"] > 0 
                          and pulse_dict_oneAmp["chi"] > loglog_line(pulse_dict_oneAmp["amplitude"], n0_spike_points_below ) 
                          and pulse_dict_oneAmp["chi"] < loglog_line(pulse_dict_oneAmp["amplitude"], n0_spike_points_above ) ):   
                    #spike band
                        spike_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        spike_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        spike_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        spike_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        spike_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        spike_chis_spike.append(pulse_dict_spike["chi"])
                        spike_amp_spike.append(pulse_dict_spike["amplitude"])
                    
                    elif (pulse_dict_oneAmp["chi"] > 2000 
                          and pulse_dict_oneAmp["amplitude"] > 0 
                          and pulse_dict_oneAmp["chi"] > loglog_line(pulse_dict_oneAmp["amplitude"], n0_PU_points_below )
                          and pulse_dict_oneAmp["chi"] < loglog_line(pulse_dict_oneAmp["amplitude"], n0_PU_points_above ) ):
                   #PU
                        PU_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        PU_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        PU_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        PU_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        PU_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        PU_chis_spike.append(pulse_dict_spike["chi"])
                        PU_amp_spike.append(pulse_dict_spike["amplitude"])
                    elif pulse_dict_oneAmp["amplitude"] > 0:
                        #misc
                        misc_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        misc_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        misc_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        misc_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        misc_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        misc_chis_spike.append(pulse_dict_spike["chi"])
                        misc_amp_spike.append(pulse_dict_spike["amplitude"])
                elif source=="Barium":
                    if (pulse_dict_oneAmp["chi"] < cut_parab(pulse_dict_oneAmp["amplitude"], (a1,c1)) 
                        and pulse_dict_oneAmp["amplitude"] > 0):
                    #signal band
                        signal_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        signal_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        signal_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        signal_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        signal_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        signal_chis_spike.append(pulse_dict_spike["chi"])
                        signal_amp_spike.append(pulse_dict_spike["amplitude"])
                        
                        signalBand_twoamp_pulse_dict_list.append(pulse_dict_twoAmp)

                        
                    elif (pulse_dict_oneAmp["chi"] > 2000 
                          and pulse_dict_oneAmp["amplitude"] > 0 
                          and pulse_dict_oneAmp["chi"] > loglog_line(pulse_dict_oneAmp["amplitude"], ba_spike_points_below  )):
                    #spike band
                        spike_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        spike_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        spike_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        spike_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        spike_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        spike_chis_spike.append(pulse_dict_spike["chi"])
                        spike_amp_spike.append(pulse_dict_spike["amplitude"])
                    
                    elif (pulse_dict_oneAmp["chi"] > 2000 
                          and pulse_dict_oneAmp["amplitude"] > 0 
                          and pulse_dict_oneAmp["chi"] > cut_parab(pulse_dict_oneAmp["amplitude"], (a1,c1))
                          and pulse_dict_oneAmp["chi"] < loglog_line(pulse_dict_oneAmp["amplitude"], ba_PU_points_above ) ):
                   #PU
                        PU_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        PU_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        PU_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        PU_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        PU_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        PU_chis_spike.append(pulse_dict_spike["chi"])
                        PU_amp_spike.append(pulse_dict_spike["amplitude"])
                    elif pulse_dict_oneAmp["amplitude"] > 0:
                        #misc
                        misc_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        misc_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        misc_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        misc_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        misc_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        misc_chis_spike.append(pulse_dict_spike["chi"])
                        misc_amp_spike.append(pulse_dict_spike["amplitude"])
                elif source=="Fe":
                    if (pulse_dict_oneAmp["chi"] < cut_parab(pulse_dict_oneAmp["amplitude"], (a2,c2)) 
                        and pulse_dict_oneAmp["amplitude"] > 0):
                    #signal band
                        signal_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        signal_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        signal_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        signal_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        signal_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        signal_chis_spike.append(pulse_dict_spike["chi"])
                        signal_amp_spike.append(pulse_dict_spike["amplitude"])
                        
                        signalBand_twoamp_pulse_dict_list.append(pulse_dict_twoAmp)

                        
                    elif (pulse_dict_oneAmp["chi"] > 10**4 
                          and pulse_dict_oneAmp["amplitude"] > 0 
                          and pulse_dict_oneAmp["chi"] > loglog_line(pulse_dict_oneAmp["amplitude"], fe_spike_points_below  )):
                    #spike band
                        spike_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        spike_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        spike_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        spike_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        spike_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        spike_chis_spike.append(pulse_dict_spike["chi"])
                        spike_amp_spike.append(pulse_dict_spike["amplitude"])
                    
                    elif (pulse_dict_oneAmp["chi"] > 2000 
                          and pulse_dict_oneAmp["amplitude"] > 0 
                          and pulse_dict_oneAmp["chi"] > cut_parab(pulse_dict_oneAmp["amplitude"], (a2,c2))
                          and pulse_dict_oneAmp["chi"] < loglog_line(pulse_dict_oneAmp["amplitude"], fe_PU_points_above ) ):
                   #PU
                        PU_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        PU_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        PU_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        PU_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        PU_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        PU_chis_spike.append(pulse_dict_spike["chi"])
                        PU_amp_spike.append(pulse_dict_spike["amplitude"])
                    elif pulse_dict_oneAmp["amplitude"] > 0:
                        #misc
                        misc_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        misc_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        misc_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
                        misc_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
                        misc_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        misc_chis_spike.append(pulse_dict_spike["chi"])
                        misc_amp_spike.append(pulse_dict_spike["amplitude"])
                
                        
        fig = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , signal_chis_oneAmp, '.', c='red', markeredgecolor='none')
        ax.plot(spike_amp_oneAmp, spike_chis_oneAmp, '.', c='orange', markeredgecolor='none')
        ax.plot(PU_amp_oneAmp , PU_chis_oneAmp, '.', c='blue', markeredgecolor='none')
        ax.plot(misc_amp_oneAmp , misc_chis_oneAmp, '.', c='black', markeredgecolor='none')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Amplitude')
        ax.set_ylabel('Chi Sq.')
        ax.minorticks_on()
#        x_plot = np.logspace(10**1, 1*10**3, 250 )
#        ax.plot(x_plot, cut_parab(x_plot, (a1, c1)), color='red')
#        ax.plot(x_plot, loglog_line(x_plot, n0_PU_points_below), color='blue')
#        ax.plot(x_plot, loglog_line(x_plot, n0_PU_points_above), color='green')
#        ax.plot(x_plot, loglog_line(x_plot, n0_spike_points_above), color='orange')
    #    ax.set_ylim([10**1, 10*np.max(chis)])
        pl.show()
        
        fig1 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , np.asarray(signal_chis_twoAmp) - np.asarray(signal_chis_oneAmp)  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_amp_oneAmp, np.asarray(spike_chis_twoAmp) - np.asarray(spike_chis_oneAmp), '.', c='orange', markeredgecolor='none')
        ax.plot(PU_amp_oneAmp , np.asarray(PU_chis_twoAmp) - np.asarray(PU_chis_oneAmp), '.', c='blue', markeredgecolor='none')
        ax.plot(misc_amp_oneAmp , np.asarray(misc_chis_twoAmp) - np.asarray(misc_chis_oneAmp), '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
        ax.set_xlabel('Amplitude of One Amplitude Fit')
        ax.set_ylabel('Delta Chi Sq. (Two Amp - One Amp)')
        ax.minorticks_on()
#        ax.set_ylim([-500*2, 500*2])
        ax.set_xlim([0, 200])
        pl.show()
        
        fig3 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , np.asarray(signal_chis_twoAmp) - np.asarray(signal_chis_oneAmp)  , '.', c='red', markeredgecolor='none')
        ax.set_xlabel('Amplitude of One Amplitude Fit')
        ax.set_ylabel('Delta Chi Sq. (Two Amp - One Amp)')
        ax.minorticks_on()
#        ax.set_ylim([-500*2, 500*2])
        ax.set_xlim([0, 200])
        pl.show()
        
        fig2 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , np.asarray(signal_chis_spike) - np.asarray(signal_chis_oneAmp)  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_amp_oneAmp, np.asarray(spike_chis_spike) - np.asarray(spike_chis_oneAmp), '.', c='orange', markeredgecolor='none')
        ax.plot(PU_amp_oneAmp , np.asarray(PU_chis_spike) - np.asarray(PU_chis_oneAmp), '.', c='blue', markeredgecolor='none')
        ax.plot(misc_amp_oneAmp , np.asarray(misc_chis_spike) - np.asarray(misc_chis_oneAmp), '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
        ax.set_xlabel('Amplitude of One Amplitude Fit')
        ax.set_ylabel('Delta Chi Sq. (Spike - One Amp)')
        ax.minorticks_on()
        ax.set_ylim([-3*10**2, 2*10**2])
        ax.set_xlim([0, 20])
        pl.show()
        
        fig5 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , np.asarray(signal_chis_spike) - np.asarray(signal_chis_oneAmp)  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_amp_oneAmp, np.asarray(spike_chis_spike) - np.asarray(spike_chis_oneAmp), '.', c='orange', markeredgecolor='none')
        ax.plot(PU_amp_oneAmp , np.asarray(PU_chis_spike) - np.asarray(PU_chis_oneAmp), '.', c='blue', markeredgecolor='none')
        ax.plot(misc_amp_oneAmp , np.asarray(misc_chis_spike) - np.asarray(misc_chis_oneAmp), '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
        ax.set_xlabel('Amplitude of One Amplitude Fit')
        ax.set_ylabel('Delta Chi Sq. (Spike - One Amp)')
        ax.minorticks_on()
        ax.set_ylim([-2*10**3, 2*10**3])
        ax.set_xlim([0, 200])
        pl.show()
        
        fig6 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_slowAmp_twoAmp , signal_fastAmp_twoAmp  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_slowAmp_twoAmp, spike_fastAmp_twoAmp, '.', c='orange', markeredgecolor='none')
        ax.plot(PU_slowAmp_twoAmp ,  PU_fastAmp_twoAmp, '.', c='blue', markeredgecolor='none')
        ax.plot(misc_slowAmp_twoAmp , misc_fastAmp_twoAmp, '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
        ax.set_ylabel('Slow Amplitude in Two Template Fit')
        ax.set_xlabel('Spike Amplitude in Two Template Fit')
        ax.minorticks_on()
        ax.set_ylim([0, 500])
        ax.set_xlim([0, 200])
        pl.show()
        
        fig6 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_slowAmp_twoAmp , signal_fastAmp_twoAmp  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_slowAmp_twoAmp, spike_fastAmp_twoAmp, '.', c='orange', markeredgecolor='none')
        ax.plot(PU_slowAmp_twoAmp ,  PU_fastAmp_twoAmp, '.', c='blue', markeredgecolor='none')
        ax.plot(misc_slowAmp_twoAmp , misc_fastAmp_twoAmp, '.', c='black', markeredgecolor='none')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylabel('Slow Amplitude in Two Template Fit')
        ax.set_xlabel('Spike Amplitude in Two Template Fit')
        ax.minorticks_on()
#        ax.set_ylim([0, 500])
#        ax.set_xlim([0, 200])
        pl.show()
        
        fig6 = pl.figure()
        ax = pl.gca()
#        ax.plot(signal_slowAmp_twoAmp , signal_fastAmp_twoAmp  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_slowAmp_twoAmp, spike_fastAmp_twoAmp, '.', c='orange', markeredgecolor='none')
#        ax.plot(PU_slowAmp_twoAmp ,  PU_fastAmp_twoAmp, '.', c='blue', markeredgecolor='none')
#        ax.plot(misc_slowAmp_twoAmp , misc_fastAmp_twoAmp, '.', c='black', markeredgecolor='none')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylabel('Slow Amplitude in Two Template Fit')
        ax.set_xlabel('Spike Amplitude in Two Template Fit')
        ax.minorticks_on()
#        ax.set_ylim([0, 500])
#        ax.set_xlim([0, 200])
        pl.show()
        
        
        fig7 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_slowAmp_twoAmp , signal_fastAmp_twoAmp  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_slowAmp_twoAmp, spike_fastAmp_twoAmp, '.', c='orange', markeredgecolor='none')
        ax.plot(PU_slowAmp_twoAmp ,  PU_fastAmp_twoAmp, '.', c='blue', markeredgecolor='none')
        ax.plot(misc_slowAmp_twoAmp , misc_fastAmp_twoAmp, '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
        ax.set_ylabel('Slow Amplitude in Two Template Fit')
        ax.set_xlabel('Spike Amplitude in Two Template Fit')
        ax.minorticks_on()
        ax.set_ylim([0, 500])
        ax.set_xlim([0, 5*10**2])
        pl.show()
        
        signal_slowAmp_twoAmp = np.asarray(signal_slowAmp_twoAmp)
        signal_fastAmp_twoAmp = np.asarray(signal_fastAmp_twoAmp)
        spike_slowAmp_twoAmp = np.asarray(spike_slowAmp_twoAmp)
        spike_fastAmp_twoAmp = np.asarray(spike_fastAmp_twoAmp)
        PU_slowAmp_twoAmp = np.asarray(PU_slowAmp_twoAmp)
        PU_fastAmp_twoAmp = np.asarray(PU_fastAmp_twoAmp)
        misc_slowAmp_twoAmp = np.asarray(misc_slowAmp_twoAmp)
        misc_fastAmp_twoAmp = np.asarray(misc_fastAmp_twoAmp)
        
        fig8 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_slowAmp_twoAmp , signal_fastAmp_twoAmp/signal_slowAmp_twoAmp  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_slowAmp_twoAmp, spike_fastAmp_twoAmp/spike_slowAmp_twoAmp, '.', c='orange', markeredgecolor='none')
        ax.plot(PU_slowAmp_twoAmp ,  PU_fastAmp_twoAmp/PU_slowAmp_twoAmp, '.', c='blue', markeredgecolor='none')
        ax.plot(misc_slowAmp_twoAmp , misc_fastAmp_twoAmp/misc_slowAmp_twoAmp, '.', c='black', markeredgecolor='none')
        ax.set_yscale('log')
#        ax.set_xscale('log')
        ax.set_ylabel('Slow Amp/Spike Amp in Two Template Fit')
        ax.set_xlabel('Spike Amplitude in Two Template Fit')
        ax.minorticks_on()
        ax.set_ylim([10**-5, 10**3])
        ax.set_xlim([0, 5*10**2])
        pl.show()
        

        fig9 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_slowAmp_twoAmp/signal_fastAmp_twoAmp , signal_fastAmp_twoAmp*1.092 + signal_slowAmp_twoAmp*0.00133  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_slowAmp_twoAmp/spike_fastAmp_twoAmp, spike_fastAmp_twoAmp*1.092 + spike_slowAmp_twoAmp*0.00133 , '.', c='orange', markeredgecolor='none')
        ax.plot(PU_slowAmp_twoAmp/PU_fastAmp_twoAmp ,  PU_fastAmp_twoAmp*1.092 + PU_slowAmp_twoAmp*0.00133 , '.', c='blue', markeredgecolor='none')
        ax.plot(misc_slowAmp_twoAmp/misc_fastAmp_twoAmp , misc_fastAmp_twoAmp*1.092 + misc_slowAmp_twoAmp*0.00133 , '.', c='black', markeredgecolor='none')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Spike Amp/Slow Amp in Two Template Fit')
        ax.set_ylabel('Integral of Two Template Fit')
        ax.minorticks_on()
        ax.set_ylim([10**-1, 10**3])
        ax.set_xlim([10**-2, 10**4])
        pl.show()
        
        fig89 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_fastAmp_twoAmp*1.092 + signal_slowAmp_twoAmp*0.00133, signal_chis_twoAmp  , '.', c='red', markeredgecolor='none')
        ax.plot( spike_fastAmp_twoAmp*1.092 + spike_slowAmp_twoAmp*0.00133, spike_chis_twoAmp , '.', c='orange', markeredgecolor='none')
        ax.plot(  PU_fastAmp_twoAmp*1.092 + PU_slowAmp_twoAmp*0.00133 , PU_chis_twoAmp, '.', c='blue', markeredgecolor='none')
        ax.plot( misc_fastAmp_twoAmp*1.092 + misc_slowAmp_twoAmp*0.00133, misc_chis_twoAmp , '.', c='black', markeredgecolor='none')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylabel('Chi Sq. of Two Template Fit')
        ax.set_xlabel('Integral of Two Template Fit')
        ax.minorticks_on()
#        ax.set_ylim([10**-1, 10**3])
#        ax.set_xlim([10**-2, 10**4])
        pl.show()
        
        
        print "now spike study"
        amp_hist = pl.hist(spike_amp_spike, range=(0,500), bins = 50, color='black')
        pl.xlabel('Amplitude [A.D.U.]')
        pl.ylabel('Events')
#        pl.xlim([0,200])
        pl.show()    
        
        amp_hist = pl.hist(spike_amp_spike, range=(0,10**3), bins = 100, color='black')
        pl.xlabel('Amplitude [A.D.U.]')
        pl.ylabel('Events')
#        pl.xlim([0,200])
        pl.show()   
        
        amp_hist2 = pl.hist(spike_amp_spike, range=(0,10**3), bins = 50, color='black')
        pl.xlabel('Amplitude [A.D.U.]')
        pl.ylabel('Events')
#        pl.xlim([0,200])
        pl.show()    
        
        fig17 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_spike , signal_chis_spike, '.', c='red', markeredgecolor='none')
        ax.plot(spike_amp_spike, spike_chis_spike, '.', c='orange', markeredgecolor='none')
        ax.plot(PU_amp_spike , PU_chis_spike, '.', c='blue', markeredgecolor='none')
        ax.plot(misc_amp_spike , misc_chis_spike, '.', c='black', markeredgecolor='none')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Amplitude')
        ax.set_ylabel('Chi Sq. of Spike Fit')
        ax.minorticks_on()
        ax.set_xlim([10**0, 10**4])
        ax.set_ylim([10**2, 10**7])
        pl.show()
        
        fig18 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_spike , signal_chis_spike, '.', c='red', markeredgecolor='none')
        ax.plot(spike_amp_spike, spike_chis_spike, '.', c='orange', markeredgecolor='none')
        ax.plot(PU_amp_spike , PU_chis_spike, '.', c='blue', markeredgecolor='none')
        ax.plot(misc_amp_spike , misc_chis_spike, '.', c='black', markeredgecolor='none')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Amplitude')
        ax.set_ylabel('Chi Sq. of Spike Fit')
        ax.minorticks_on()
#        ax.set_xlim([10**0, 10**4])
#        ax.set_ylim([10**2, 10**7])
        pl.show()
                
        print "Total pulses", (len(signal_amp_oneAmp) + len(spike_amp_oneAmp) + len(PU_amp_oneAmp)+ len(misc_amp_oneAmp))
#        analyzeFitResults([signalBand_twoamp_pulse_dict_list])
                    
def PCAtemplate(list_ofPulse_Dict_Lists, list_ofNoisePSDs = []):
    chis=[]
    peaks=[]
    obs_pulses=[]
    templates=[]
    slow_amps=[]
    fast_amps=[]
    times=[]
    amps=[]
    delays=[]
    errors=[]
    
    obs_fouriers = []
    weighted_fouriers = []
    
    list_ofNoisePSDs_sqrts = []
    for i in list_ofNoisePSDs:
        list_ofNoisePSDs_sqrts.append(np.sqrt(i))
    
    for r in range(len(list_ofPulse_Dict_Lists)):
        run = list_ofPulse_Dict_Lists[r]
        for result in run: #result = a pulse
            chis.append(result.get("chi")    )
            peaks.append(result.get("peak")    )
            obs_pulses.append(result.get("obs pulse") )
            templates.append(result.get("template") )
            slow_amps.append(result.get("slow amplitude") )
            fast_amps.append(result.get("fast amplitude") )
            amps.append(result.get("amplitude") )
            times.append( result.get("time")   )
            delays.append(result.get("delay") )
            errors.append(result.get("amp error") )
            
            obs_fourier = np.fft.fft(result.get("obs pulse"))
#            weighted_fourier = obs_fourier/(list_ofNoisePSDs_sqrts[r])
#            
            obs_fouriers.append(obs_fourier)
#            weighted_fouriers.append(weighted_fourier)
#            
    
    if len(obs_pulses) == 0:
        print "NO PULSES FOUND"
        return
    
    save_principal_0 = []   
    
    print "\n PCA Freq Domain centered\n"
    mean_fourier = np.mean(obs_fouriers, axis=0) 
    centered_fourier_matrix = obs_fouriers - mean_fourier
    (u, s, vh) = np.linalg.svd(centered_fourier_matrix)
    print "hi2"
    print "lengths", len(vh[0]), len(u[:,0])
    print np.shape(u), np.shape(vh)
    one_second_interval = np.linspace(0,1, int(fe))
    print type(u), type(s), type(vh)
    print type(vh[0]), type(vh[0][0])
    
    for i in range(6):
        standard_template = theta_exp(one_second_interval, (0.24, 0.008,1.1))
        standard_template_max = max(standard_template)
        
        vector = np.real(np.fft.ifft(vh[i]))
        zeroed_vector = vector - vector[0]
        if max(np.abs(zeroed_vector)) == max(zeroed_vector):
            scaling = standard_template_max/max(zeroed_vector)
        else:
            scaling = standard_template_max/min(zeroed_vector)
        scaled_vector = scaling*zeroed_vector
#        scaled_vector[int(0.2*fe):] = 0
        try:
            if i==0:
                fit = curve_fit(A_theta_exp, one_second_interval, scaled_vector, (1.0, 0.24, 0.008,1.1))
#                fit = curve_fit(A_theta_exp, one_second_interval, scaled_vector, (1.0, 0.14, 0.008,0.04))
                print fit
                save_principal_0 = scaled_vector
                pl.plot(one_second_interval, scaled_vector, label = 'Principal Component 0', color='b')
                pl.plot(one_second_interval, standard_template, label='Initial Guess')
                pl.plot(one_second_interval, fit[0][0]*theta_exp(one_second_interval, fit[0][1:]) , label='Best Fit', color='r')
                pl.legend()
                pl.show()
                print "RSS error of best fit", sum( (fit[0][0]*theta_exp(one_second_interval, fit[0][1:]) - scaled_vector)**2 )
                print "RSS error of initial guess fit", sum( (standard_template - scaled_vector)**2 )
                print "Best Params", fit[0]
                print "Best rise time: ", fit[0][2], " +- ", np.sqrt(fit[1][2][2])
                print "Best Fall time: ", fit[0][3], " +- ", np.sqrt(fit[1][3][3])
                print "Cov of Rise Time", fit[1][2]
                print "Cov of Fall Time", fit[1][3]
                print "Cov Matrix", fit[1]
                
                fit2 = curve_fit(AB_theta_exp, one_second_interval, scaled_vector, (0.24, 0.008,1.1, 2.5, 0.5, 0.5), 
                bounds=([0.23, 0, 0, 0, 0, 0], [0.25, 0.1 , 2, 3, 1, 1]))
                
                pl.plot(one_second_interval, scaled_vector, label = 'Principal Component 0', color='b')
                pl.plot(one_second_interval, AB_theta_exp(one_second_interval, 0.24, 0.008,1.1, 2.5, 0.5, 0.5), label='Initial Guess 2 Amp')
                pl.plot(one_second_interval, AB_theta_exp(one_second_interval, fit2[0][0], fit2[0][1], fit2[0][2], fit2[0][3], fit2[0][4], fit2[0][5])  , label='Best Fit 2 Amp', color='r')
                pl.legend()
                pl.show()
                print "RSS error of best fit", sum( (AB_theta_exp(one_second_interval, fit2[0][0], fit2[0][1], fit2[0][2], fit2[0][3], fit2[0][4], fit2[0][5]) - scaled_vector)**2 )
                print "Best Params", fit2[0]
                print "Best rise time: ", fit2[0][1], " +- ", np.sqrt(fit2[1][1][1])
                print "Best Fall time 1: ", fit2[0][2], " +- ", np.sqrt(fit2[1][2][2])
                print "Best Fall time 2: ", fit2[0][3], " +- ", np.sqrt(fit2[1][3][3]) 
                print "Cov Matrix", fit2[1]
        except:
            print "Fit failed on principal component 0"
#            
        pl.plot(one_second_interval, scaled_vector, label= "Principal Component " + str(i) , color = 'b')
        pl.plot(one_second_interval, standard_template , label='Template', color = 'r')
        pl.xlabel('Time [s]')
        pl.ylabel('Amplitude [ADU]')
        pl.legend()
        pl.show()
        
#    pl.plot(one_second_interval, np.fft.ifft(vh[0]))
#    pl.show()
        
    lambdas = s**2
    print lambdas[0:20]/(sum(lambdas[0:20]))
#    
#    print "\n PCA Noise Weighted centered\n"
#    mean_fourier = np.mean(weighted_fouriers, axis=0) 
#    centered_fourier_matrix = weighted_fouriers - mean_fourier
#    (u, s, vh) = np.linalg.svd(centered_fourier_matrix)
#    print "hi2"
#    print "lengths", len(vh[0]), len(u[:,0])
#    print np.shape(u), np.shape(vh)
#    one_second_interval = np.linspace(0,1, int(fe))
#    print type(u), type(s), type(vh)
#    print type(vh[0]), type(vh[0][0])
#    
#    save_principal_1 = []
#    for i in range(6):
#        standard_template = theta_exp(one_second_interval, (0.24, 0.008,1.1))
#        standard_template_max = max(standard_template)
#        
#        vector = np.real(np.fft.ifft(vh[i]))
#        zeroed_vector = vector - vector[0]
#        if max(np.abs(zeroed_vector)) == max(zeroed_vector):
#            scaling = standard_template_max/max(zeroed_vector)
#        else:
#            scaling = standard_template_max/min(zeroed_vector)
#        scaled_vector = scaling*zeroed_vector
##        scaled_vector[int(0.2*fe):] = 0
#        try:
#            if i==0:
#                fit = curve_fit(A_theta_exp, one_second_interval, scaled_vector, (1.0, 0.24, 0.008,1.1))
##                fit = curve_fit(A_theta_exp, one_second_interval, scaled_vector, (1.0, 0.14, 0.008,0.04))
#                print fit
#                save_principal_1 = scaled_vector
#                pl.plot(one_second_interval, scaled_vector, label = 'Principal Component 0', color='b')
#                pl.plot(one_second_interval, standard_template, label='Initial Guess')
#                pl.plot(one_second_interval, fit[0][0]*theta_exp(one_second_interval, fit[0][1:]) , label='Best Fit', color='r')
#                pl.legend()
#                pl.show()
#                print "RSS error of best fit", sum( (fit[0][0]*theta_exp(one_second_interval, fit[0][1:]) - scaled_vector)**2 )
#                print "RSS error of initial guess fit", sum( (standard_template - scaled_vector)**2 )
#                print "Best Params", fit[0]
#                print "Best rise time: ", fit[0][2], " +- ", np.sqrt(fit[1][2][2])
#                print "Best Fall time: ", fit[0][3], " +- ", np.sqrt(fit[1][3][3])
#                print "Cov of Rise Time", fit[1][2]
#                print "Cov of Fall Time", fit[1][3]
#                print "Cov Matrix", fit[1]
#                
#                fit2 = curve_fit(AB_theta_exp, one_second_interval, scaled_vector, (0.24, 0.008,1.1, 2.5, 0.5, 0.5), 
#                bounds=([0.23, 0, 0, 0, 0, 0], [0.25, 0.1 , 2, 3, 1, 1]))
#                
#                pl.plot(one_second_interval, scaled_vector, label = 'Principal Component 0', color='b')
#                pl.plot(one_second_interval, AB_theta_exp(one_second_interval, 0.24, 0.008,1.1, 2.5, 0.5, 0.5), label='Initial Guess 2 Amp')
#                pl.plot(one_second_interval, AB_theta_exp(one_second_interval, fit2[0][0], fit2[0][1], fit2[0][2], fit2[0][3], fit2[0][4], fit2[0][5])  , label='Best Fit 2 Amp', color='r')
#                pl.legend()
#                pl.show()
#                print "RSS error of best fit", sum( (AB_theta_exp(one_second_interval, fit2[0][0], fit2[0][1], fit2[0][2], fit2[0][3], fit2[0][4], fit2[0][5]) - scaled_vector)**2 )
#                print "Best Params", fit2[0]
#                print "Best rise time: ", fit2[0][1], " +- ", np.sqrt(fit2[1][1][1])
#                print "Best Fall time 1: ", fit2[0][2], " +- ", np.sqrt(fit2[1][2][2])
#                print "Best Fall time 2: ", fit2[0][3], " +- ", np.sqrt(fit2[1][3][3]) 
#                print "Cov Matrix", fit2[1]
#        except:
#            print "Fit failed on principal component 0"
##            
#        pl.plot(one_second_interval, scaled_vector, label= "Principal Component " + str(i) , color = 'b')
#        pl.plot(one_second_interval, standard_template , label='Template', color = 'r')
#        pl.xlabel('Time [s]')
#        pl.ylabel('Amplitude [ADU]')
#        pl.legend()
#        pl.show()
#    
#    (u, s, vh) = np.linalg.svd(centered_pulse_matrix)
    return save_principal_0
#    print "\n PCA Freq Domain Centered\n"
    

def newFitComparison(list_ofPulse_Dict_Lists, list_ofNoisePSDs , templateFunc):
    chis=[]
    peaks=[]
    obs_pulses=[]
    templates=[]
    slow_amps=[]
    fast_amps=[]
    times=[]
    amps=[]
    delays=[]
    errors=[]
    
    one_sec = np.linspace(0,1,1000)
    
    obs_fouriers = []
    weighted_fouriers = []
    for r in range(len(list_ofPulse_Dict_Lists)):
        run = list_ofPulse_Dict_Lists[r]
        for result in run: #result = a pulse
            chis.append(result.get("chi")    )
            peaks.append(result.get("peak")    )
            obs_pulses.append(result.get("obs pulse") )
            templates.append(result.get("template") )
            slow_amps.append(result.get("slow amplitude") )
            fast_amps.append(result.get("fast amplitude") )
            amps.append(result.get("amplitude") )
            times.append( result.get("time")   )
            delays.append(result.get("delay") )
            errors.append(result.get("amp error") )
            
#            obs_fourier = np.fft.fft(result.get("obs pulse"))
#            weighted_fourier = obs_fourier/np.sqrt(list_ofNoisePSDs[r])
#            
#            obs_fouriers.append(obs_fourier)
#            weighted_fouriers.append(weighted_fourier)
#            
#            for t in range 
#                template_vals ...
#                get_windowed_fft...
#                Coeffs
#                a_hat
#                chi_squared
#                if chi_sq is min:
#                    a_hat0 = ahat
    
    

    


    

#file_name = 'data_run42_dbz1/20180315_14h24'; source = "Fe"; #No source#
file_name = 'data_run42_dbz1/20180313_18h32'; source = "Barium"; #good Ba 14 chunks run
#file_name = 'data_run42_dbz1/20180314_16h12'; source = "Neutron";  #Neutrons 16 chunks
#source = "Neutron" #or "Barium", Neutron, or Fe
chunk_size = 1
chunk_number = 14
if (source == "Barium"): chunk_number = 14 #total number of chunks
elif (source == "Neutron"): chunk_number = 16
    
#folder_header = 'Results/' + file_name + '/' +str(chunk_size) + 'hours_'  #original runs
#folder_header = 'Results/' + 'standardized_starts/'+ file_name + '/' +str(chunk_size) + 'hours_'  #runs w standardized pulse start times
folder_header = 'Results/' + 'test0/'+ file_name + '/' +str(chunk_size) + 'hours_'  #runs w standardized pulse start times
#folder_header = 'Results/' + 'with_spike2/'+ file_name + '/' +str(chunk_size) + 'hours_'  #runs w standardized pulse start times


a1 = 0.02
c1 = 800
a2 = 0.12
c2 = 800
#neutron ss
n0_PU_points_below = (200.,3000., 3000., 2*10**6) 
n0_PU_points_above = (60.,5000., 1000., 3*10**6)
n0_spike_points_below = (60.,5000., 1000., 3*10**6) 
n0_spike_points_above = (10, 1500, 10**3, 1.5*10**7)
#ba ss
ba_spike_points_below = (10.,2000., 100., 1.5*10**5)
ba_PU_points_above = (20.,2000., 200., 2*10**5)
#Fe ss
fe_spike_points_below = (10.,8*10**3., 100., 8*10**5)
fe_PU_points_above = (10.,2000., 100., 1.5*10**5)


list_oneAmpResults = []
list_oneAmpResults_posAmps = []
list_oneAmpResults_signal = []
list_oneAmpResults_spike = []
list_oneAmpResults_PU = []
list_NoisePSDs = []

list_twoAmpResults = []
list_spikeResults = []
list_oneAmpResults_blob = []
list_new_twoAmpResults = []
list_new_oneAmpResults = []


##list_freqs = []
for i in range(chunk_number):
    directory = folder_header + str(i) + '/'
#    directory = 'Results/data_run42_dbz1/test/'
#    with open(directory + 'oneAMP_pulse_processing_results.p', 'rb') as fp1:
#        oneAmp_processed_results = pickle.load(fp1)
    with open(directory + 'full_pulse_processing_results.p', 'rb') as fp1:
        (oneAmp_processed_results, twoAmp_processed_results, spike_processed_results, newT_oneamp, newT_twoamp) = pickle.load(fp1)
#        (oneAmp_processed_results, twoAmp_processed_results, newT_oneamp , spike_processed_results, bq) = pickle.load(fp1)
#        (oneAmp_processed_results, twoAmp_processed_results, spike_processed_results) = pickle.load(fp1)
#        (oneAmp_processed_results, twoAmp_processed_results) = pickle.load(fp1)
    with open(directory + 'NoisePSD.p', 'rb') as fp2:
        (freq, J) = pickle.load(fp2)
     
    #Neutron Cuts    
    #signal band   
    if source=="Neutron":
        oneAmp_posAmps = getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf)  ])
        oneAmp_signal = generalizedCut(oneAmp_posAmps, "amplitude", "chi", cut_parab, (a1,c1) , True)
        #pile-up events?
        oneAmp_PU= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 2000, np.inf)  ])
        oneAmp_PU= generalizedCut(oneAmp_PU, "amplitude", "chi", loglog_line, n0_PU_points_below, veto_greater_than_func= False )
        oneAmp_PU= generalizedCut(oneAmp_PU, "amplitude", "chi", loglog_line, n0_PU_points_above, veto_greater_than_func= True)
        #Spike band events?
    #    print len(oneAmp_processed_results)
        oneAmp_spike= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 2000, np.inf)  ])
        oneAmp_spike= generalizedCut(oneAmp_spike, "amplitude", "chi", loglog_line, n0_spike_points_below, veto_greater_than_func= False )
        oneAmp_spike= generalizedCut(oneAmp_spike, "amplitude", "chi", loglog_line, n0_spike_points_above, veto_greater_than_func= True)
    
    #Ba Cuts    
    #signal band   
    if source=="Barium":
        oneAmp_posAmps = getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf)  ])
        oneAmp_signal = generalizedCut(oneAmp_posAmps, "amplitude", "chi", cut_parab, (a1,c1) , True)
        #pile-up events?
        oneAmp_PU= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 2000, np.inf)  ])
        oneAmp_PU= generalizedCut(oneAmp_PU, "amplitude", "chi", loglog_line, ba_PU_points_above, veto_greater_than_func= True )
        oneAmp_PU = generalizedCut(oneAmp_PU, "amplitude", "chi", cut_parab, (a1,c1) , False)
        #Spike band events?
        oneAmp_spike= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 2000, np.inf)  ])
        oneAmp_spike= generalizedCut(oneAmp_spike, "amplitude", "chi", loglog_line, ba_spike_points_below, veto_greater_than_func= False )
        #blob
        oneAmp_blob = generalizedCut(oneAmp_posAmps, "amplitude", "chi", cut_parab, (a1,c1) , False)
        oneAmp_blob = getPulseSubset(oneAmp_blob, [("amplitude", 7, 10.5), ("chi", 800, 2000)  ])
    if source=="Fe":
        oneAmp_posAmps = getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf)  ])
        oneAmp_signal = generalizedCut(oneAmp_posAmps, "amplitude", "chi", cut_parab, (a2,c2) , True)
        #pile-up events?
        oneAmp_PU= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 2000, np.inf)  ])
        oneAmp_PU= generalizedCut(oneAmp_PU, "amplitude", "chi", loglog_line, fe_PU_points_above, veto_greater_than_func= True )
        oneAmp_PU = generalizedCut(oneAmp_PU, "amplitude", "chi", cut_parab, (a2,c2) , False)
        #Spike band events?
        oneAmp_spike= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 10**4, np.inf)  ])
        oneAmp_spike= generalizedCut(oneAmp_spike, "amplitude", "chi", loglog_line, fe_spike_points_below, veto_greater_than_func= False )


    list_oneAmpResults.append(oneAmp_processed_results)
    list_oneAmpResults_posAmps.append(oneAmp_posAmps)
    list_oneAmpResults_signal.append(oneAmp_signal)
    list_oneAmpResults_PU.append(oneAmp_PU)
    list_oneAmpResults_spike.append(oneAmp_spike)
    if source=="Barium":
        list_oneAmpResults_blob.append(oneAmp_blob)
    list_NoisePSDs.append( J )
    
    list_twoAmpResults.append(twoAmp_processed_results)
    list_spikeResults.append(spike_processed_results)
    
    list_new_oneAmpResults.append(newT_oneamp)
    list_new_twoAmpResults.append(newT_twoamp) 



#Tr_fit = 0.008
#Tf_fit = 1.1
#print "====================================================="
#print "REFIT TEST"
#print "====================================================="
##reFit_optimalFiltering1amp(list_oneAmpResults, list_NoisePSDs )
##assert False

print "====================================================="
print "ALL POSITIVE AMPLITUDE PULSES"
print "====================================================="
analyzeFitResults(list_oneAmpResults_posAmps)
##compareFitResults(list_oneAmpResults, list_twoAmpResults, list_spikeResults)
compareFitResults(list_oneAmpResults, list_new_twoAmpResults, list_spikeResults)
#PCAtemplate(list_oneAmpResults_posAmps)
#spectralClusterPulses(list_oneAmpResults_posAmps)
#print "====================================================="
#print "SIGNAL BAND PULSES ONLY"
#print "====================================================="
#analyzeFitResults(list_oneAmpResults_signal)
#
#
#analyzeFitResults(list_oneAmpResults_signal)
#pc0_signal = PCAtemplate(list_oneAmpResults_signal)
##
##
###print "====================================================="
###print "POTENTIAL PILE UP PULSES ONLY"
###print "====================================================="
###analyzeFitResults(list_oneAmpResults_PU)
###PCAtemplate(list_oneAmpResults_PU)
###t_plot = np.linspace(0,1,int(fe))
###counter = 0
###for r in list_oneAmpResults_PU:
###    for i in r:
###        pl.plot(t_plot, i["obs pulse"], label='Pulse number ' + str(counter) )
###        pl.legend()
###        pl.show()
###        counter+=1
#print "====================================================="
#print "POTENTIAL SPIKE BAND PULSES ONLY"
#print "====================================================="
#analyzeFitResults(list_oneAmpResults_spike)
#pc0_spike = PCAtemplate(list_oneAmpResults_spike)
#with open(source + "_pc0s.p", 'wb') as fp:
#        pickle.dump( (pc0_signal, pc0_spike) , fp )

#print "====================================================="
#print "Blob ONLY"
#print "====================================================="
#analyzeFitResults(list_oneAmpResults_blob)
#PCAtemplate(list_oneAmpResults_blob)
#counter = 0
#for r in list_oneAmpResults_blob:
#    for i in r:
#        pl.plot(t_plot, i["obs pulse"], label='Pulse number ' + str(counter) )
#        pl.legend()
#        pl.show()
#        counter+=1
        
#Compare principal components
#with open('Nu_pc0s.p', 'rb') as fp:
#        (Nu_pc0_signal, Nu_pc0_spike) = pickle.load(fp)
#with open('Ba_pc0s.p', 'rb') as fp:
#        (Ba_pc0_signal, Ba_pc0_spike) = pickle.load(fp)    
#with open('Fe_pc0s.p', 'rb') as fp:
#        (Fe_pc0_signal, Fe_pc0_spike) = pickle.load(fp)   
#one_sec = np.linspace(0, 1, int(fe))
#pl.plot(one_sec, Nu_pc0_signal, label="Neutron Signal PC")
#pl.plot(one_sec, Ba_pc0_signal, label="Barium Signal PC")
#pl.xlabel('Time [s]')
#pl.ylabel('Amplitude')
#pl.legend()
#pl.show()
#pl.plot(one_sec, Nu_pc0_signal, label="Neutron Signal PC")
#pl.plot(one_sec, Ba_pc0_signal, label="Barium Signal PC")
#pl.plot(one_sec, Fe_pc0_signal, label="Fe Signal PC", color='red')
#pl.xlabel('Time [s]')
#pl.ylabel('Amplitude')
#pl.legend()
#pl.show()
#pl.plot(one_sec, Nu_pc0_signal - Ba_pc0_signal, label="Neutron - Barium Signal PC")
#pl.xlabel('Time [s]')
#pl.ylabel('Amplitude')
#pl.legend()
#pl.show()
#
#pl.plot(one_sec, Nu_pc0_signal - Fe_pc0_signal, label="Neutron - Fe Signal PC")
#pl.xlabel('Time [s]')
#pl.ylabel('Amplitude')
#pl.legend()
#pl.show()
#
#pl.plot(one_sec, Fe_pc0_signal - Ba_pc0_signal , label="Fe - Ba Signal PC")
#pl.xlabel('Time [s]')
#pl.ylabel('Amplitude')
#pl.legend()
#pl.show()
#
#fit_test = curve_fit(A_theta_exp, one_sec, Nu_pc0_signal - Ba_pc0_signal , (0.35, 0.24, 0.0007514, 0.0093262  )) 
#pl.plot(one_sec, Nu_pc0_signal - Ba_pc0_signal, label="Neutron - Barium Signal PC")
#pl.plot(one_sec, A_theta_exp(one_sec, 0.35, 0.24, 0.0007514, 0.0093262  ), label = 'Spike Template', color = 'black')
##pl.plot(one_sec, A_theta_exp(one_sec, fit_test[0][0], fit_test[0][1], fit_test[0][2], fit_test[0][3] ), label = 'Best Fit', color= 'red')
#pl.xlabel('Time [s]')
#pl.ylabel('Amplitude')
#pl.legend()
#pl.show()
#print fit_test
#
#pl.plot(one_sec, Nu_pc0_signal/Ba_pc0_signal, label="Neutron/Barium Signal PC")
#pl.xlabel('Time [s]')
#pl.legend()
#pl.show()
#pl.plot(one_sec, Nu_pc0_signal/Ba_pc0_signal, label="Neutron/Barium Signal PC")
#pl.xlabel('Time [s]')
#pl.xlim([0.25,1.0])
#pl.ylim([0,1.5])
#pl.legend()
#pl.show()
#pl.plot(one_sec, Nu_pc0_spike, label="Neutron Spike PC")
#pl.plot(one_sec, Ba_pc0_spike, label="Barium Spike PC")
#pl.xlabel('Time [s]')
#pl.ylabel('Amplitude')
#pl.legend()
#pl.show()
#index_nu = list(Nu_pc0_spike).index(max(Nu_pc0_spike)) 
#index_ba = list(Ba_pc0_spike).index(max(Ba_pc0_spike)) 
#index_fe = list(Fe_pc0_spike).index(max(Fe_pc0_spike)) 
#
#diff = index_ba - index_nu
#modified_nu_spike = diff*[0] + list(Nu_pc0_spike)
#modified_nu_spike = modified_nu_spike[0:int(fe)]
#diff_Fe = index_ba - index_fe
#modified_fe_spike = diff_Fe*[0] + list(Fe_pc0_spike)
#modified_fe_spike = modified_fe_spike[0:int(fe)]
#pl.plot(one_sec, modified_nu_spike, label="Shifted Neutron Spike PC")
#pl.plot(one_sec, Ba_pc0_spike, label="Barium Spike PC")
#pl.plot(one_sec, modified_fe_spike, label="Shifted Fe Spike PC", color='red')
#pl.xlabel('Time [s]')
#pl.ylabel('Amplitude')
#pl.legend()
#pl.show()
#pl.plot(one_sec, Nu_pc0_spike-Ba_pc0_spike, label="Neutron - Barium Spike PC")
#pl.legend()
#pl.show()
#pl.plot(one_sec, np.asarray(modified_nu_spike)-Ba_pc0_spike, label="Shifted Neutron - Barium Spike PC")
#pl.legend()
#pl.show()
#pl.plot(one_sec, np.asarray(modified_nu_spike)-np.asarray(modified_fe_spike), label="Shifted Neutron - Shifted Fe Spike PC")
#pl.legend()
#pl.show()
#pl.plot(one_sec, np.asarray(modified_fe_spike)-Ba_pc0_spike, label="Shifted Fe - Barium Spike PC")
#pl.legend()
#pl.show()



