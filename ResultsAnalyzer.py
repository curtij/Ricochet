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
import pandas as pd
#font =11.5
font =15

#font =22.5
#print pl.rcParams.keys()
pl.rcParams.update({'font.size': font})
#pl.rc('text', usetex=False)
pl.rc('figure' , figsize=(10,5)) #dpi=300
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
l    returns frequencies, FFT
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

def highPassFilter(timeseries, filtersize=2, order = 2, inverse=False):
        '''
        filter size in Hz
        '''
    
        b, a =sgl.butter(order,2*filtersize/fe,btype="highpass")
        if inverse==False:
            z=lfilter(b,a,timeseries)
        
        if inverse==True:
            z=lfilter(a,b,timeseries)
        
        return z
    
def fourier1to2sided(fourieroneside,N):
#    onesidelist=list(fourieroneside)
    
    twosideMemMap = np.memmap('TempFiles/fourier1to2sided.mymemmap', dtype=np.complex128, mode='r+', shape=(N,))
    templen = len(fourieroneside)
    twosideMemMap[0:templen] = fourieroneside[:]
    for k in range(0,int(N-len(fourieroneside))):
        v=int(N-templen)
#        onesidelist.append(np.conjugate(fourieroneside[v]))
        twosideMemMap[templen] = np.conjugate(fourieroneside[v])
        templen += 1
        
#    Fouriertwoside=np.array(onesidelist)
    
#    return Fouriertwoside
    
    return twosideMemMap


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


    amp_hist_2 = pl.hist(amps, range=(0,200), bins = 100, color='black')
    pl.xlabel('Amplitude [A.D.U.]')
    pl.ylabel('Events')
    pl.xlim([0,200])
    pl.show()    
    
    amp_hist_1 = pl.hist(amps, range=(0,200), bins = 200, color='black')
    pl.xlabel('Amplitude [A.D.U.]')
    pl.ylabel('Events')
    pl.xlim([0,200])
    pl.show()   
    
    amp_hist = pl.hist(amps, range=(0,200), bins = 40, color='black')
    pl.xlabel('Amplitude [A.D.U.]')
    pl.ylabel('Events')
    pl.xlim([0,200])
    pl.show()   
    
    amp_hist3 = pl.hist(amps, range=(0,200), bins = 100, color='black', log=True)
    pl.xlabel('Amplitude [A.D.U.]')
    pl.ylabel('Events')
    pl.xlim([0,200])
    pl.show()    
    
    with open(source + "_amp_spectrum.p", 'wb') as fp:
            pickle.dump( (amp_hist_1, amp_hist_2)  , fp )
    
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

def ratio_theta_exp(t, a, t_0, Tr1, Tf1, Tr2, Tf2, ratio):
    return A_theta_exp(t, a, t_0, Tr1, Tf1) +  A_theta_exp(t, a*ratio, t_0, Tr2, Tf2)

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
    '''1 one amp, 2 two amp, 3 spike'''
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

#    pl.scatter(amps, deltaChi, '.', c='black')
#    pl.xlabel('Amplitude')
#    pl.xlim([0,200])
#    pl.ylim([-2000,2000])
#    pl.ylabel('Delta Chi Sq')
#    pl.show()
    
    fig = pl.figure()
    ax = pl.gca()
    ax.plot(amps , deltaChi, '.', c='black', markeredgecolor='none')
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Chi Sq. of Original Fit - Chi Sq. of Fe Fit')
    ax.set_ylim([-2*10**3, 2*10**3])
    ax.set_xlim([-100, 300])
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
#                        signal_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        signal_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
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
#                        spike_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        spike_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
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
#                        PU_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        PU_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        PU_chis_spike.append(pulse_dict_spike["chi"])
                        PU_amp_spike.append(pulse_dict_spike["amplitude"])
                    elif pulse_dict_oneAmp["amplitude"] > 0:
                        #misc
                        misc_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        misc_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        misc_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
#                        misc_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        misc_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        misc_chis_spike.append(pulse_dict_spike["chi"])
                        misc_amp_spike.append(pulse_dict_spike["amplitude"])
                elif source=="Barium":
                    if (pulse_dict_oneAmp["chi"] < cut_parab(pulse_dict_oneAmp["amplitude"], (a1,c1)) 
                        and pulse_dict_oneAmp["amplitude"] > 0):
                    #signal band
                        signal_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        signal_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        signal_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
#                        signal_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        signal_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
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
#                        spike_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        spike_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
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
#                        PU_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        PU_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        PU_chis_spike.append(pulse_dict_spike["chi"])
                        PU_amp_spike.append(pulse_dict_spike["amplitude"])
                    elif pulse_dict_oneAmp["amplitude"] > 0:
                        #misc
                        misc_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        misc_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        misc_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
#                        misc_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        misc_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        misc_chis_spike.append(pulse_dict_spike["chi"])
                        misc_amp_spike.append(pulse_dict_spike["amplitude"])
                elif source=="Fe":
                    if (pulse_dict_oneAmp["chi"] < cut_parab(pulse_dict_oneAmp["amplitude"], (a2,c2)) 
                        and pulse_dict_oneAmp["amplitude"] > 0):
                    #signal band
                        signal_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        signal_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        signal_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
#                        signal_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        signal_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
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
#                        spike_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        spike_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
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
#                        PU_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        PU_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
                        PU_chis_spike.append(pulse_dict_spike["chi"])
                        PU_amp_spike.append(pulse_dict_spike["amplitude"])
                    elif pulse_dict_oneAmp["amplitude"] > 0:
                        #misc
                        misc_chis_oneAmp.append(pulse_dict_oneAmp["chi"])
                        misc_amp_oneAmp.append(pulse_dict_oneAmp["amplitude"])
                        
                        misc_chis_twoAmp.append(pulse_dict_twoAmp["chi"])
#                        misc_slowAmp_twoAmp.append(pulse_dict_twoAmp["slow amplitude"])
#                        misc_fastAmp_twoAmp.append(pulse_dict_twoAmp["fast amplitude"])
                                            
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
        
        fig = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , signal_chis_spike, '.', c='red', markeredgecolor='none')
        ax.plot(spike_amp_oneAmp, spike_chis_spike, '.', c='orange', markeredgecolor='none')
        ax.plot(PU_amp_oneAmp , PU_chis_spike, '.', c='blue', markeredgecolor='none')
        ax.plot(misc_amp_oneAmp , misc_chis_spike, '.', c='black', markeredgecolor='none')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Amplitude')
        ax.set_ylabel('Chi Sq. of Ba/Neutron Fit')
        ax.minorticks_on()
        pl.show()
        
        fig = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , signal_chis_twoAmp, '.', c='red', markeredgecolor='none')
        ax.plot(spike_amp_oneAmp, spike_chis_twoAmp, '.', c='orange', markeredgecolor='none')
        ax.plot(PU_amp_oneAmp , PU_chis_twoAmp, '.', c='blue', markeredgecolor='none')
        ax.plot(misc_amp_oneAmp , misc_chis_twoAmp, '.', c='black', markeredgecolor='none')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Amplitude')
        ax.set_ylabel('Chi Sq. of Fe Fit')
        ax.minorticks_on()
        pl.show()
        
        fig1 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , np.asarray(signal_chis_twoAmp) - np.asarray(signal_chis_oneAmp)  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_amp_oneAmp, np.asarray(spike_chis_twoAmp) - np.asarray(spike_chis_oneAmp), '.', c='orange', markeredgecolor='none')
        ax.plot(PU_amp_oneAmp , np.asarray(PU_chis_twoAmp) - np.asarray(PU_chis_oneAmp), '.', c='blue', markeredgecolor='none')
        ax.plot(misc_amp_oneAmp , np.asarray(misc_chis_twoAmp) - np.asarray(misc_chis_oneAmp), '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
        ax.set_xlabel('Amplitude of Original Fit')
        ax.set_ylabel('Delta Chi Sq. (Fe - One Amp)')
        ax.minorticks_on()
        ax.set_ylim([-500*2, 500*2])
        ax.set_xlim([0, 200])
        pl.show()
        
        fig3 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , np.asarray(signal_chis_twoAmp) - np.asarray(signal_chis_oneAmp)  , '.', c='red', markeredgecolor='none')
        ax.set_xlabel('Amplitude of Original Fit')
        ax.set_ylabel('Delta Chi Sq. (Fe - One Amp)')
        ax.minorticks_on()
        ax.set_ylim([-500*2, 500*2])
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
        ax.set_xlabel('Amplitude of Original Fit')
        ax.set_ylabel('Delta Chi Sq. (Ba/Neutron - One Amp)')
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
        ax.set_xlabel('Amplitude of Original Fit')
        ax.set_ylabel('Delta Chi Sq. (Ba/Neutron - One Amp)')
        ax.minorticks_on()
        ax.set_ylim([-2*10**3, 2*10**3])
        ax.set_xlim([0, 200])
        pl.show()
        
        fig5 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , np.asarray(signal_chis_spike) - np.asarray(signal_chis_twoAmp)  , '.', c='red', markeredgecolor='none')
        ax.plot(spike_amp_oneAmp, np.asarray(spike_chis_spike) - np.asarray(spike_chis_twoAmp), '.', c='orange', markeredgecolor='none')
        ax.plot(PU_amp_oneAmp , np.asarray(PU_chis_spike) - np.asarray(PU_chis_twoAmp), '.', c='blue', markeredgecolor='none')
        ax.plot(misc_amp_oneAmp , np.asarray(misc_chis_spike) - np.asarray(misc_chis_twoAmp), '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
        ax.set_xlabel('Amplitude of Original Fit')
        ax.set_ylabel('Delta Chi Sq. (Ba/Neutron - Fe Fits)')
        ax.minorticks_on()
        ax.set_ylim([-2*10**3, 2*10**3])
        ax.set_xlim([0, 200])
        pl.show()
        
        fig5 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , np.asarray(signal_chis_spike) - np.asarray(signal_chis_twoAmp)  , '.', c='red', markeredgecolor='none')
#        ax.plot(spike_amp_oneAmp, np.asarray(spike_chis_spike) - np.asarray(spike_chis_twoAmp), '.', c='orange', markeredgecolor='none')
#        ax.plot(PU_amp_oneAmp , np.asarray(PU_chis_spike) - np.asarray(PU_chis_twoAmp), '.', c='blue', markeredgecolor='none')
#        ax.plot(misc_amp_oneAmp , np.asarray(misc_chis_spike) - np.asarray(misc_chis_twoAmp), '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
        ax.set_xlabel('Amplitude of Original Fit')
        ax.set_ylabel('Delta Chi Sq. (Ba - Neutron Fits)')
        ax.minorticks_on()
        ax.set_ylim([-0.2*10**3, 0.1*10**3])
        ax.set_xlim([0, 200])
        pl.show()
        
        fig5 = pl.figure()
        ax = pl.gca()
        ax.plot(signal_amp_oneAmp , np.asarray(signal_chis_spike) - np.asarray(signal_chis_twoAmp)  , '.', c='red', markeredgecolor='none')
#        ax.plot(spike_amp_oneAmp, np.asarray(spike_chis_spike) - np.asarray(spike_chis_twoAmp), '.', c='orange', markeredgecolor='none')
#        ax.plot(PU_amp_oneAmp , np.asarray(PU_chis_spike) - np.asarray(PU_chis_twoAmp), '.', c='blue', markeredgecolor='none')
#        ax.plot(misc_amp_oneAmp , np.asarray(misc_chis_spike) - np.asarray(misc_chis_twoAmp), '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
        ax.set_xlabel('Amplitude of Original Fit')
        ax.set_ylabel('Delta Chi Sq. (Ba - Neutron Fits)')
        ax.minorticks_on()
        ax.set_ylim([-2*10**3, 2*10**3])
        ax.set_xlim([0, 200])
        pl.show()
        
        
        delta = np.asarray(signal_chis_spike) - np.asarray(signal_chis_twoAmp)
        with open(source + "_deltachiBa_Fe.p", 'wb') as fp:
            pickle.dump( (signal_amp_oneAmp, delta )  , fp )
#        
        
        
        
#        fig6 = pl.figure()
#        ax = pl.gca()
#        ax.plot(signal_slowAmp_twoAmp , signal_fastAmp_twoAmp  , '.', c='red', markeredgecolor='none')
#        ax.plot(spike_slowAmp_twoAmp, spike_fastAmp_twoAmp, '.', c='orange', markeredgecolor='none')
#        ax.plot(PU_slowAmp_twoAmp ,  PU_fastAmp_twoAmp, '.', c='blue', markeredgecolor='none')
#        ax.plot(misc_slowAmp_twoAmp , misc_fastAmp_twoAmp, '.', c='black', markeredgecolor='none')
##        ax.set_yscale('log')
##        ax.set_xscale('log')
#        ax.set_ylabel('Slow Amplitude in Two Template Fit')
#        ax.set_xlabel('Spike Amplitude in Two Template Fit')
#        ax.minorticks_on()
#        ax.set_ylim([0, 500])
#        ax.set_xlim([0, 200])
#        pl.show()
#        
#        fig6 = pl.figure()
#        ax = pl.gca()
#        ax.plot(signal_slowAmp_twoAmp , signal_fastAmp_twoAmp  , '.', c='red', markeredgecolor='none')
#        ax.plot(spike_slowAmp_twoAmp, spike_fastAmp_twoAmp, '.', c='orange', markeredgecolor='none')
#        ax.plot(PU_slowAmp_twoAmp ,  PU_fastAmp_twoAmp, '.', c='blue', markeredgecolor='none')
#        ax.plot(misc_slowAmp_twoAmp , misc_fastAmp_twoAmp, '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
#        ax.set_ylabel('Slow Amplitude in Two Template Fit')
#        ax.set_xlabel('Spike Amplitude in Two Template Fit')
#        ax.minorticks_on()
##        ax.set_ylim([0, 500])
##        ax.set_xlim([0, 200])
#        pl.show()
#        
#        fig6 = pl.figure()
#        ax = pl.gca()
##        ax.plot(signal_slowAmp_twoAmp , signal_fastAmp_twoAmp  , '.', c='red', markeredgecolor='none')
#        ax.plot(spike_slowAmp_twoAmp, spike_fastAmp_twoAmp, '.', c='orange', markeredgecolor='none')
##        ax.plot(PU_slowAmp_twoAmp ,  PU_fastAmp_twoAmp, '.', c='blue', markeredgecolor='none')
##        ax.plot(misc_slowAmp_twoAmp , misc_fastAmp_twoAmp, '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
#        ax.set_ylabel('Slow Amplitude in Two Template Fit')
#        ax.set_xlabel('Spike Amplitude in Two Template Fit')
#        ax.minorticks_on()
##        ax.set_ylim([0, 500])
##        ax.set_xlim([0, 200])
#        pl.show()
#        
#        
#        fig7 = pl.figure()
#        ax = pl.gca()
#        ax.plot(signal_slowAmp_twoAmp , signal_fastAmp_twoAmp  , '.', c='red', markeredgecolor='none')
#        ax.plot(spike_slowAmp_twoAmp, spike_fastAmp_twoAmp, '.', c='orange', markeredgecolor='none')
#        ax.plot(PU_slowAmp_twoAmp ,  PU_fastAmp_twoAmp, '.', c='blue', markeredgecolor='none')
#        ax.plot(misc_slowAmp_twoAmp , misc_fastAmp_twoAmp, '.', c='black', markeredgecolor='none')
##        ax.set_yscale('log')
##        ax.set_xscale('log')
#        ax.set_ylabel('Slow Amplitude in Two Template Fit')
#        ax.set_xlabel('Spike Amplitude in Two Template Fit')
#        ax.minorticks_on()
#        ax.set_ylim([0, 500])
#        ax.set_xlim([0, 5*10**2])
#        pl.show()
#        
#        signal_slowAmp_twoAmp = np.asarray(signal_slowAmp_twoAmp)
#        signal_fastAmp_twoAmp = np.asarray(signal_fastAmp_twoAmp)
#        spike_slowAmp_twoAmp = np.asarray(spike_slowAmp_twoAmp)
#        spike_fastAmp_twoAmp = np.asarray(spike_fastAmp_twoAmp)
#        PU_slowAmp_twoAmp = np.asarray(PU_slowAmp_twoAmp)
#        PU_fastAmp_twoAmp = np.asarray(PU_fastAmp_twoAmp)
#        misc_slowAmp_twoAmp = np.asarray(misc_slowAmp_twoAmp)
#        misc_fastAmp_twoAmp = np.asarray(misc_fastAmp_twoAmp)
#        
#        fig8 = pl.figure()
#        ax = pl.gca()
#        ax.plot(signal_slowAmp_twoAmp , signal_fastAmp_twoAmp/signal_slowAmp_twoAmp  , '.', c='red', markeredgecolor='none')
#        ax.plot(spike_slowAmp_twoAmp, spike_fastAmp_twoAmp/spike_slowAmp_twoAmp, '.', c='orange', markeredgecolor='none')
#        ax.plot(PU_slowAmp_twoAmp ,  PU_fastAmp_twoAmp/PU_slowAmp_twoAmp, '.', c='blue', markeredgecolor='none')
#        ax.plot(misc_slowAmp_twoAmp , misc_fastAmp_twoAmp/misc_slowAmp_twoAmp, '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
##        ax.set_xscale('log')
#        ax.set_ylabel('Slow Amp/Spike Amp in Two Template Fit')
#        ax.set_xlabel('Spike Amplitude in Two Template Fit')
#        ax.minorticks_on()
#        ax.set_ylim([10**-5, 10**3])
#        ax.set_xlim([0, 5*10**2])
#        pl.show()
#        
##        with open(source + "_ampratio.p", 'wb') as fp:
##            pickle.dump( (signal_slowAmp_twoAmp, (signal_fastAmp_twoAmp/signal_slowAmp_twoAmp ) )  , fp )
#        
#
#        fig9 = pl.figure()
#        ax = pl.gca()
#        ax.plot(signal_slowAmp_twoAmp/signal_fastAmp_twoAmp , signal_fastAmp_twoAmp*1.092 + signal_slowAmp_twoAmp*0.00133  , '.', c='red', markeredgecolor='none')
#        ax.plot(spike_slowAmp_twoAmp/spike_fastAmp_twoAmp, spike_fastAmp_twoAmp*1.092 + spike_slowAmp_twoAmp*0.00133 , '.', c='orange', markeredgecolor='none')
#        ax.plot(PU_slowAmp_twoAmp/PU_fastAmp_twoAmp ,  PU_fastAmp_twoAmp*1.092 + PU_slowAmp_twoAmp*0.00133 , '.', c='blue', markeredgecolor='none')
#        ax.plot(misc_slowAmp_twoAmp/misc_fastAmp_twoAmp , misc_fastAmp_twoAmp*1.092 + misc_slowAmp_twoAmp*0.00133 , '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
#        ax.set_xlabel('Spike Amp/Slow Amp in Two Template Fit')
#        ax.set_ylabel('Integral of Two Template Fit')
#        ax.minorticks_on()
#        ax.set_ylim([10**-1, 10**3])
#        ax.set_xlim([10**-2, 10**4])
#        pl.show()
#        
#        fig89 = pl.figure()
#        ax = pl.gca()
#        ax.plot( signal_fastAmp_twoAmp*1.092 + signal_slowAmp_twoAmp*0.00133, signal_chis_twoAmp  , '.', c='red', markeredgecolor='none')
#        ax.plot( spike_fastAmp_twoAmp*1.092 + spike_slowAmp_twoAmp*0.00133, spike_chis_twoAmp , '.', c='orange', markeredgecolor='none')
#        ax.plot(  PU_fastAmp_twoAmp*1.092 + PU_slowAmp_twoAmp*0.00133 , PU_chis_twoAmp, '.', c='blue', markeredgecolor='none')
#        ax.plot( misc_fastAmp_twoAmp*1.092 + misc_slowAmp_twoAmp*0.00133, misc_chis_twoAmp , '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
#        ax.set_ylabel('Chi Sq. of Two Template Fit')
#        ax.set_xlabel('Integral of Two Template Fit')
#        ax.minorticks_on()
##        ax.set_ylim([10**-1, 10**3])
##        ax.set_xlim([10**-2, 10**4])
#        pl.show()
#        
#        
#        print "now spike study"
#        amp_hist = pl.hist(spike_amp_spike, range=(0,500), bins = 50, color='black')
#        pl.xlabel('Amplitude [A.D.U.]')
#        pl.ylabel('Events')
##        pl.xlim([0,200])
#        pl.show()    
#        
#        amp_hist = pl.hist(spike_amp_spike, range=(0,10**3), bins = 100, color='black')
#        pl.xlabel('Amplitude [A.D.U.]')
#        pl.ylabel('Events')
##        pl.xlim([0,200])
#        pl.show()   
#        
#        amp_hist2 = pl.hist(spike_amp_spike, range=(0,10**3), bins = 50, color='black')
#        pl.xlabel('Amplitude [A.D.U.]')
#        pl.ylabel('Events')
##        pl.xlim([0,200])
#        pl.show()    
#        
#        fig17 = pl.figure()
#        ax = pl.gca()
#        ax.plot(signal_amp_spike , signal_chis_spike, '.', c='red', markeredgecolor='none')
#        ax.plot(spike_amp_spike, spike_chis_spike, '.', c='orange', markeredgecolor='none')
#        ax.plot(PU_amp_spike , PU_chis_spike, '.', c='blue', markeredgecolor='none')
#        ax.plot(misc_amp_spike , misc_chis_spike, '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
#        ax.set_xlabel('Amplitude')
#        ax.set_ylabel('Chi Sq. of Spike Fit')
#        ax.minorticks_on()
#        ax.set_xlim([10**0, 10**4])
#        ax.set_ylim([10**2, 10**7])
#        pl.show()
#        
#        fig18 = pl.figure()
#        ax = pl.gca()
#        ax.plot(signal_amp_spike , signal_chis_spike, '.', c='red', markeredgecolor='none')
#        ax.plot(spike_amp_spike, spike_chis_spike, '.', c='orange', markeredgecolor='none')
#        ax.plot(PU_amp_spike , PU_chis_spike, '.', c='blue', markeredgecolor='none')
#        ax.plot(misc_amp_spike , misc_chis_spike, '.', c='black', markeredgecolor='none')
#        ax.set_yscale('log')
#        ax.set_xscale('log')
#        ax.set_xlabel('Amplitude')
#        ax.set_ylabel('Chi Sq. of Spike Fit')
#        ax.minorticks_on()
##        ax.set_xlim([10**0, 10**4])
##        ax.set_ylim([10**2, 10**7])
#        pl.show()
                
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
            
#            obs_fourier = np.fft.fft(result.get("obs pulse"))
            obs_fourier = np.fft.rfft( highPassFilter(result.get("obs pulse")) )
#            obs_fourier =  getWindowedFFT(result.get("obs pulse"))
            weighted_fourier = obs_fourier/(list_ofNoisePSDs_sqrts[r])
#            
            obs_fouriers.append(obs_fourier)
            weighted_fouriers.append(weighted_fourier)
            
#            
    print "weighted shape", np.shape(weighted_fouriers), np.shape(weighted_fouriers[0])
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
        
        vector = highPassFilter(np.real(np.fft.irfft(vh[i])), inverse=True)
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
                avg_pulse = highPassFilter(np.fft.irfft(mean_fourier), inverse=True)
                avg_pulse = avg_pulse - avg_pulse[0]
                avg_pulse = max(scaled_vector)/max(avg_pulse)*avg_pulse
                pl.plot(one_second_interval, avg_pulse, label='average in freq domain')
#                pl.plot(one_second_interval, AB_theta_exp(one_second_interval, 0.24, 0.008,1.1, 2.5, 0.5, 0.5), label='Initial Guess 2 Amp')
#                pl.plot(one_second_interval, AB_theta_exp(one_second_interval, fit2[0][0], fit2[0][1], fit2[0][2], fit2[0][3], fit2[0][4], fit2[0][5])  , label='Best Fit 2 Amp', color='r')
                pl.legend()
                pl.show()
                
                pl.plot(one_second_interval, scaled_vector - avg_pulse, label = 'Principal Component 0 - avg', color='b')
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
        assert False

#            
        pl.plot(one_second_interval, scaled_vector, label= "Principal Component " + str(i) , color = 'b')
        pl.plot(one_second_interval, standard_template , label='Template', color = 'r')
        pl.xlabel('Time [s]')
        pl.ylabel('Amplitude [ADU]')
        pl.legend()
        pl.show()
        
#    pl.plot(one_second_interval, np.fft.irfft(vh[0]))
#    pl.show()
#    print "==================="
#    COV_SN = np.cov(centered_fourier_matrix, rowvar=False)
#    COV_SN_diagonal = np.diagonal(COV_SN)
#    mean_amp = np.mean(amps)
#    print "mean amplitude", mean_amp
#    std_template_psd = scipy.signal.periodogram( (mean_amp*standard_template), fe)[1]
#    pl.loglog(freq, COV_SN_diagonal, label="diagonal of pulse cov")
#    pl.loglog(freq, std_template_psd, label="Template PSD")
#    pl.loglog(freq, std_template_psd + list_ofNoisePSDs[0], label="Template PSD + noise PSD")
#    pl.loglog(freq, (COV_SN_diagonal*std_template_psd[0]/COV_SN_diagonal[0]) , label="scaled diagonal of pulse cov")
#    pl.legend()
#    pl.ylim([10**-5,10**9])
#    pl.show()
#    pl.loglog(freq, list_ofNoisePSDs[0], label="Noise PSD")
#    pl.legend()
#    pl.ylim([10**-5,10**9])
#    pl.show()
#    
#    pl.loglog(freq, COV_SN_diagonal/list_ofNoisePSDs[0], label="Ratio")
#    pl.legend()
##    pl.ylim([10**-5,10**9])
#    pl.show()
#    
#    pl.loglog(freq, COV_SN_diagonal-list_ofNoisePSDs[0], label="difference")
#    pl.legend()
##    pl.ylim([10**-5,10**9])
#    pl.show()
#    
#    pgrams = []
#    for i in centered_fourier_matrix:
#        psd_temp = (np.abs(i)**2)*2/(fe*(fe+1))
#        pgrams.append(psd_temp)
#    mean_psd_pulse = np.mean(pgrams, axis=0)
#    
#    pl.loglog(freq, mean_psd_pulse, label="Pulse PSD")
#    pl.loglog(freq, std_template_psd + list_ofNoisePSDs[0], label="Template PSD + noise PSD")
#    pl.legend()
#    pl.ylim([10**-5,10**4])
#    pl.show()
#    
#    pl.loglog(freq, mean_psd_pulse -  list_ofNoisePSDs[0], label="Pulse PSD - Noise PSD")
#    pl.loglog(freq, std_template_psd , label="Template PSD")
#    pl.legend()
#    pl.ylim([10**-10,10**4])
#    pl.show()
#    print "==================="
#
#    
#    lambdas = s**2
#    print lambdas[0:20]/(sum(lambdas[0:20]))
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
#        
##        vector = np.real(np.fft.ifft( fourier1to2sided( (list_ofNoisePSDs_sqrts[0]*vh[i]), 1000) ) )
#        vector = np.fft.irfft( (np.mean(list_ofNoisePSDs_sqrts))*vh[i] ) 
#        zeroed_vector = vector - vector[int(0.1*fe)]
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
##    
##    (u, s, vh) = np.linalg.svd(centered_pulse_matrix)
#    pl.plot(one_second_interval, save_principal_0, label='Unweighted PC0')
#    pl.plot(one_second_interval, save_principal_1, label='Weighted PC0')
#    pl.legend()
#    pl.xlabel("Time [s]")
#    pl.ylabel("Amplitude")
#    pl.show()
#    
#    pl.plot(one_second_interval, save_principal_0 - save_principal_1, label='Difference in PC0')
#    pl.legend()
#    pl.xlabel("Time [s]")
#    pl.ylabel("Amplitude")
#    pl.show()
#    lambdas = s**2
#    print lambdas[0:20]/(sum(lambdas[0:20]))
#    return save_principal_0, save_principal_1
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
    
    
    new_amps = []
    new_chis = []
    
    obs_fouriers = []
    weighted_fouriers = []
    for r in range(len(list_ofPulse_Dict_Lists)):
        run = list_ofPulse_Dict_Lists[r]
        J = list_ofNoisePSDs[r]
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
            
            obs_fourier = np.fft.rfft(result.get("obs pulse"))
#            weighted_fourier = obs_fourier/np.sqrt(list_ofNoisePSDs[r])
#            
            obs_fouriers.append(obs_fourier)
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
            
            pulse = result.get("obs pulse")
            original_template = result.get("template")
            original_amp = result.get("amplitude")
            original_chi = result.get("chi")    

            
            one_sec = np.linspace(0,1,1000)
#            temp_hp = (templateFunc)
#            temp_hp = highPassFilter(theta_exp(one_sec, (0.24, 0.008, 1.1)) )
            start = list(original_template).index(filter(lambda x: x>(1-10**-14)*original_template[0], original_template)[0])
            delay = result.get("delay") + 0.25 - 0.1
            delay = (start - 1)/fe
            print "start", start, result.get("delay"), delay
#            templateFunc = theta_exp(one_sec, (delay, 0.008, 1.1))
#            templateFunc = removeDC(templateFunc)
            f, S = getWindowedFFT(highPassFilter(templateFunc))
            f, V = getWindowedFFT(highPassFilter(pulse))
#            f, V = getWindowedFFT((pulse))
#            S = np.fft.rfft(templateFunc)
#            V = np.fft.rfft(pulse)

            
            to_plot = True
            if to_plot:
                pl.loglog(f, np.abs(getWindowedFFT(highPassFilter(original_template))[1])**2, label='og template')
                pl.loglog(f, np.abs(V)**2, label='V')
                pl.loglog(f, J, label='J')
                pl.loglog(f, np.abs(original_amp*S)**2, label='S')
                pl.legend()
                pl.ylim([10**-10, 10**4])
                pl.show()
                
                pl.loglog(f, np.abs(original_amp*S)**2 - np.abs(getWindowedFFT(highPassFilter(original_template))[1])**2, label='S - og temp')
                pl.legend()
                pl.ylim([10**-10, 10**4])
                pl.show()
                
                pl.plot(f, np.angle(S) - np.angle(getWindowedFFT(highPassFilter(original_template))[1]), label = 'phase difference')
                pl.legend()
                pl.show()
                
                pl.plot(one_sec, pulse)
                pl.plot(one_sec, original_template, label='original fit')
                pl.plot(one_sec, original_amp*templateFunc , label='original amp new temp')
    #            pl.plot(np.linspace(0,1,1000), ahat_t*templateFunc , label='new fit')
                pl.legend()
                pl.show()
                
                pl.plot(one_sec, original_template - original_amp*templateFunc - original_template[0] , label='og temp - temp')
                pl.legend()
                pl.show()
                assert False
            
  
            numer = sum(( np.real( np.conjugate(S) * V )/J ))
            denom = sum(np.abs(S)**2/J)
            ahat_t = float(numer)/float(denom)
            
            new_chi = sum((np.abs(V[1:] - ahat_t*S[1:] ))**2 / J[1:])  
            new_amps.append(ahat_t)
            new_chis.append(new_chi)
            

            
            print "amp, chi, new amp, new chi\n",original_amp, ahat_t, original_chi, new_chi
            print sum((np.abs(V[1:] - original_amp*S[1:] ))**2 / J[1:])
            print sum((np.abs(V[1:] - getWindowedFFT(highPassFilter(original_template))[1][1:]   )) **2 / J[1:])
#            assert False
            
            
    fig = pl.figure()
    ax = pl.gca()
    ax.plot(amps , chis, '.', c='black', markeredgecolor='none', label='original fit')
    ax.plot(new_amps , new_chis, '.', c='red', markeredgecolor='none', label= 'new fit')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Chi Sq.')
    ax.minorticks_on()
    pl.legend()
    pl.show()
    
def compare_results_new(Fe_df0, Ba_df0, AmBe_df0):
    Fe_df = Fe_df0[Fe_df0.amplitude > 0]
    Ba_df = Ba_df0[Ba_df0.amplitude > 0]
    AmBe_df = AmBe_df0[AmBe_df0.amplitude > 0]

#    Fe_normal_cut = "chi < @cut_parab(amplitude, (@a2,@c2))"
#    Ba_normal_cut = "chi < @cut_parab(@amplitude, (@a1,@c1))"
#    AmBe_normal_cut = "chi < @cut_parab(@amplitude, (@a1,@c1))"
    
    print list(Fe_df.columns)
#    Fe_normal_events = Fe_df.query(Fe_normal_cut)
#    Ba_normal_events = Ba_df.query(Ba_normal_cut)
#    AmBe_normal_events = AmBe_df.query(AmBe_normal_cut)

    Fe_normal_events = Fe_df[Fe_df.chi < cut_parab(Fe_df.amplitude, (a2,c2))]
    Ba_normal_events = Ba_df[Ba_df.chi < cut_parab(Ba_df.amplitude, (a1,c1))]
    AmBe_normal_events = AmBe_df[AmBe_df.chi < cut_parab(AmBe_df.amplitude, (a1,c1))]
    AmBe_normal_events_upperband = AmBe_normal_events[ (AmBe_normal_events.delta_chi_Ba_AmBe > 0)]
    AmBe_normal_events_lowerband = AmBe_normal_events[ (AmBe_normal_events.delta_chi_Ba_AmBe <= 0)]

    print "upper band events, lower", len(AmBe_normal_events_upperband.index), len(AmBe_normal_events_lowerband.index)
#    Fe_spike_events = Fe_df.loc( not(Fe_df["chi"] < cut_parab(a1,c1)) and  Fe_df["chi"] > 2000 
#                          and Fe_df["chi"] > loglog_line(Fe_df["amplitude"], fe_spike_points_below ) )

    fig5 = pl.figure()
    ax = pl.gca()
    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["delta_chi_Ba_AmBe"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["amplitude"],  Ba_normal_events["delta_chi_Ba_AmBe"] , '.', label='Ba Data', c='red', markeredgecolor='none')
    ax.plot(AmBe_normal_events["amplitude"],  AmBe_normal_events["delta_chi_Ba_AmBe"] , '.', label='AmBe Data', c='black', markeredgecolor='none')
    ax.set_xlabel('Amplitude of Original Fit')
    ax.set_ylabel('Delta Chi Sq. (Ba - AmBe Fits)')
    ax.minorticks_on()
    ax.set_ylim([-0.2*10**3, 0.1*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
#    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["delta_chi_Ba_AmBe"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["amplitude"],  Ba_normal_events["delta_chi_Ba_AmBe"] , '.', label='Ba Data', c='red', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_upperband["amplitude"],  AmBe_normal_events_upperband["delta_chi_Ba_AmBe"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_lowerband["amplitude"],  AmBe_normal_events_lowerband["delta_chi_Ba_AmBe"] , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of Original Fit')
    ax.set_ylabel('Delta Chi Sq. (Ba - AmBe Fits)')
    ax.minorticks_on()
    ax.set_ylim([-0.2*10**3, 0.1*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
#    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["delta_chi_Ba_AmBe"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["Ba_amp"],  Ba_normal_events["delta_chi_Ba_AmBe"] , '.', label='Ba Data', c='red', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_upperband["amplitude"],  AmBe_normal_events_upperband["delta_chi_Ba_AmBe"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_lowerband["amplitude"],  AmBe_normal_events_lowerband["delta_chi_Ba_AmBe"] , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of Ba Fit')
    ax.set_ylabel('Delta Chi Sq. (Ba - AmBe Fits)')
    ax.minorticks_on()
    ax.set_ylim([-0.2*10**3, 0.1*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
    ax.plot(Fe_normal_events["AmBe_amp"],  Fe_normal_events["delta_chi_Ba_AmBe"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["AmBe_amp"],  Ba_normal_events["delta_chi_Ba_AmBe"] , '.', label='Ba Data', c='red', markeredgecolor='none')
    ax.plot(AmBe_normal_events_upperband["AmBe_amp"],  AmBe_normal_events_upperband["delta_chi_Ba_AmBe"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
    ax.plot(AmBe_normal_events_lowerband["AmBe_amp"],  AmBe_normal_events_lowerband["delta_chi_Ba_AmBe"] , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of AmBe Fit')
    ax.set_ylabel('Delta Chi Sq. (Ba - AmBe Fits)')
    ax.minorticks_on()
    ax.set_ylim([-0.2*10**3, 0.1*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["delta_chi_Ba_Fe"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["amplitude"],  Ba_normal_events["delta_chi_Ba_Fe"] , '.', label='Ba Data', c='red', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_upperband["amplitude"],  AmBe_normal_events_upperband["delta_chi_Ba_Fe"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_lowerband["amplitude"],  AmBe_normal_events_lowerband["delta_chi_Ba_Fe"] , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of Original Fit')
    ax.set_ylabel('Delta Chi Sq. (Ba - Fe Fits)')
    ax.minorticks_on()
    ax.set_ylim([-2*10**3, 2*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["delta_chi_Fe_AmBe"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["amplitude"],  Ba_normal_events["delta_chi_Fe_AmBe"] , '.', label='Ba Data', c='red', markeredgecolor='none')
    ax.plot(AmBe_normal_events_upperband["amplitude"],  AmBe_normal_events_upperband["delta_chi_Fe_AmBe"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
    ax.plot(AmBe_normal_events_lowerband["amplitude"],  AmBe_normal_events_lowerband["delta_chi_Fe_AmBe"] , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of Original Fit')
    ax.set_ylabel('Delta Chi Sq. (Fe - AmBe Fits)')
    ax.minorticks_on()
    ax.set_ylim([-2*10**3, 2*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
    ax.plot(Fe_normal_events["delta_chi_Ba_AmBe"],  Fe_normal_events["delta_chi_Fe_AmBe"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["delta_chi_Ba_AmBe"],  Ba_normal_events["delta_chi_Fe_AmBe"] , '.', label='Ba Data', c='red', markeredgecolor='none')
    ax.plot(AmBe_normal_events_upperband["delta_chi_Ba_AmBe"],  AmBe_normal_events_upperband["delta_chi_Fe_AmBe"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
    ax.plot(AmBe_normal_events_lowerband["delta_chi_Ba_AmBe"],  AmBe_normal_events_lowerband["delta_chi_Fe_AmBe"] , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Delta Chi Sq. (Ba - AmBe Fits)')
    ax.set_ylabel('Delta Chi Sq. (Fe - AmBe Fits)')
    ax.minorticks_on()
    ax.set_ylim([-1.5*10**3, 1.5*10**3])
    ax.set_xlim([-0.5*10**3, 0.5*10**3])
    pl.legend()
    pl.show()
    
    
    
    fig5 = pl.figure()
    ax = pl.gca()
    ax.plot(Fe_normal_events["delta_chi_Ba_AmBe"],  Fe_normal_events["delta_chi_Fe_AmBe"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["delta_chi_Ba_AmBe"],  Ba_normal_events["delta_chi_Fe_AmBe"] , '.', label='Ba Data', c='red', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_upperband["delta_chi_Ba_AmBe"],  AmBe_normal_events_upperband["delta_chi_Fe_AmBe"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_lowerband["delta_chi_Ba_AmBe"],  AmBe_normal_events_lowerband["delta_chi_Fe_AmBe"] , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Delta Chi Sq. (Ba - AmBe Fits)')
    ax.set_ylabel('Delta Chi Sq. (Fe - AmBe Fits)')
    ax.minorticks_on()
    ax.set_ylim([-1.5*10**3, 1.5*10**3])
    ax.set_xlim([-0.5*10**3, 0.5*10**3])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
#    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["delta_chi_Fe_AmBe"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
#    ax.plot(Ba_normal_events["amplitude"],  Ba_normal_events["delta_chi_Fe_AmBe"] , '.', label='Ba Data', c='red', markeredgecolor='none')
    ax.loglog(AmBe_normal_events_upperband["amplitude"],  AmBe_normal_events_upperband["chi"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
    ax.loglog(AmBe_normal_events_lowerband["amplitude"],  AmBe_normal_events_lowerband["chi"] , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of Original Fit')
    ax.set_ylabel('Chi Sq. (Original Fit)')
    ax.minorticks_on()
    ax.set_ylim([10**2, 10**7])
    ax.set_xlim([10**0, 3*10**3])
    pl.legend()
    pl.show()
#    num_bins = 100
#    pl.hist(Fe_normal_events["Fe_amp"],color='black', range=(0,200), bins=num_bins, label = 'Fe data')
#    pl.ylabel("Events")
#    pl.xlabel("Amplitude of Fe Fit")
#    pl.xlim([0,200])
#    pl.legend()
#    pl.show()
#    
#    pl.hist(Ba_normal_events["Ba_amp"],color='black', range=(0,200),bins=num_bins, label = 'Ba data')
#    pl.ylabel("Events")
#    pl.xlabel("Amplitude of Ba Fit")
#    pl.xlim([0,200])
#    pl.legend()
#    pl.show()
#    
#    pl.hist(AmBe_normal_events["AmBe_amp"],color='black', range=(0,200),bins=num_bins, label = 'AmBe data')
#    pl.ylabel("Events")
#    pl.xlabel("Amplitude of AmBe Fit")
#    pl.xlim([0,200])
#    pl.legend()
#    pl.show()
#    
#    pl.hist(AmBe_normal_events_upperband["AmBe_amp"],color='black',range=(0,200), bins=num_bins, label = 'AmBe upper band')
#    pl.ylabel("Events")
#    pl.xlabel("Amplitude of AmBe Fit")
#    pl.xlim([0,200])
#    pl.legend()
#    pl.show()
#    
#    pl.hist(AmBe_normal_events_lowerband["AmBe_amp"], color='black',range=(0,200),bins=num_bins, label = 'AmBe lower band')
#    pl.ylabel("Events")
#    pl.xlabel("Amplitude of AmBe Fit")
#    pl.xlim([0,200])
#    pl.legend()
#    pl.show()
#    
#    pl.hist(AmBe_normal_events_lowerband["Ba_amp"],color='black', range=(0,200),bins=num_bins, label = 'AmBe lower band')
#    pl.ylabel("Events")
#    pl.xlabel("Amplitude of Ba Fit")
#    pl.xlim([0,200])
#    pl.legend()
#    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
#    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["Ba_chi"] - Fe_normal_events["chi"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["amplitude"],  Ba_normal_events["Ba_chi"] - Ba_normal_events["chi"], '.', label='Ba Data', c='red', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_upperband["amplitude"], AmBe_normal_events_upperband["Ba_chi"] - AmBe_normal_events_upperband["chi"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_lowerband["amplitude"],  AmBe_normal_events_lowerband["Ba_chi"] - AmBe_normal_events_lowerband["chi"]  , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of Original Fit')
    ax.set_ylabel('Delta Chi Sq. (Ba - Original Fits)')
    ax.minorticks_on()
    ax.set_ylim([-1*10**3, 1*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["Fe_chi"] - Fe_normal_events["chi"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["amplitude"],  Ba_normal_events["Fe_chi"] - Ba_normal_events["chi"], '.', label='Ba Data', c='red', markeredgecolor='none')
    ax.plot(AmBe_normal_events["amplitude"], AmBe_normal_events["Fe_chi"] - AmBe_normal_events["chi"] , '.', label='AmBe Data', c='black', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_upperband["amplitude"], AmBe_normal_events_upperband["Fe_chi"] - AmBe_normal_events_upperband["chi"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_lowerband["amplitude"],  AmBe_normal_events_lowerband["Fe_chi"] - AmBe_normal_events_lowerband["chi"]  , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of Original Fit')
    ax.set_ylabel('Delta Chi Sq. (Fe - Original Fits)')
    ax.minorticks_on()
    ax.set_ylim([-1*10**3, 1*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
#    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["AmBe_chi"] - Fe_normal_events["chi"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
#    ax.plot(Ba_normal_events["amplitude"],  Ba_normal_events["AmBe_chi"] - Ba_normal_events["chi"], '.', label='Ba Data', c='red', markeredgecolor='none')
    ax.plot(AmBe_normal_events_upperband["amplitude"], AmBe_normal_events_upperband["AmBe_chi"] - AmBe_normal_events_upperband["chi"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
    ax.plot(AmBe_normal_events_lowerband["amplitude"],  AmBe_normal_events_lowerband["AmBe_chi"] - AmBe_normal_events_lowerband["chi"]  , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of Original Fit')
    ax.set_ylabel('Delta Chi Sq. (AmBe - Original Fits)')
    ax.minorticks_on()
    ax.set_ylim([-1*10**3, 1*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["AmBe_chi"] - Fe_normal_events["chi"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["amplitude"],  Ba_normal_events["AmBe_chi"] - Ba_normal_events["chi"], '.', label='Ba Data', c='red', markeredgecolor='none')
    ax.plot(AmBe_normal_events["amplitude"], AmBe_normal_events["AmBe_chi"] - AmBe_normal_events["chi"] , '.', label='AmBe Data', c='black', markeredgecolor='none')
#    ax.plot(AmBe_normal_events_lowerband["amplitude"],  AmBe_normal_events_lowerband["AmBe_chi"] - AmBe_normal_events_lowerband["chi"]  , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of Original Fit')
    ax.set_ylabel('Delta Chi Sq. (AmBe - Original Fits)')
    ax.minorticks_on()
    ax.set_ylim([-1*10**3, 1*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    fig5 = pl.figure()
    ax = pl.gca()
    ax.plot(Fe_normal_events["amplitude"],  Fe_normal_events["AmBe_chi"] - Fe_normal_events["chi"] , '.',label='Fe Data', c='blue', markeredgecolor='none')
    ax.plot(Ba_normal_events["amplitude"],  Ba_normal_events["AmBe_chi"] - Ba_normal_events["chi"], '.', label='Ba Data', c='red', markeredgecolor='none')
    ax.plot(AmBe_normal_events_upperband["amplitude"], AmBe_normal_events_upperband["AmBe_chi"] - AmBe_normal_events_upperband["chi"] , '.', label='AmBe Upper Band', c='black', markeredgecolor='none')
    ax.plot(AmBe_normal_events_lowerband["amplitude"],  AmBe_normal_events_lowerband["AmBe_chi"] - AmBe_normal_events_lowerband["chi"]  , '.', label='AmBe Lower Band', c='orange', markeredgecolor='none')  
    ax.set_xlabel('Amplitude of Original Fit')
    ax.set_ylabel('Delta Chi Sq. (AmBe - Original Fits)')
    ax.minorticks_on()
    ax.set_ylim([-1*10**3, 1*10**3])
    ax.set_xlim([0, 200])
    pl.legend()
    pl.show()
    
    time_1sec = np.linspace(0,1,int(fe))
    upper_avg = AmBe_normal_events_upperband["obs pulse"].mean()
    upper_avg = upper_avg - upper_avg[235]
    upper_avg = upper_avg/max(upper_avg)
    lower_avg = AmBe_normal_events_lowerband["obs pulse"].mean()
    lower_avg = lower_avg - lower_avg[235] 
    lower_avg = lower_avg/max(lower_avg) + 0.0025009862143049455
    pl.plot(time_1sec, upper_avg, label='AmBe upper avg pulse')
    pl.plot(time_1sec, lower_avg, label='AmBe lower avg pulse')
    pl.legend()
    pl.show()
    
    pl.semilogy(time_1sec, np.abs(upper_avg), label='AmBe upper avg pulse')
    pl.semilogy(time_1sec, np.abs(lower_avg), label='AmBe lower avg pulse')
    pl.legend()
    pl.show()
    
    pl.semilogy(time_1sec, np.abs(upper_avg)- np.abs(lower_avg), label='AmBe upper avg pulse')
    pl.legend()
    pl.show()
    
    print np.mean( np.abs(upper_avg)[400]- np.abs(lower_avg)[400] )
    
    
#    for i in AmBe_normal_events_upperband["obs pulse"]:
#        pl.plot(np.linspace(0,1,int(fe)), i)
#        pl.show()
#    print "\n PCA upper"
#    PC_upper = PCA_pulselist(AmBe_normal_events_upperband["obs pulse"])
#    print "\n PCA lower"
#    PC_lower = PCA_pulselist(AmBe_normal_events_lowerband["obs pulse"])
#    time_1sec = np.linspace(0,1,int(fe))
#    pl.plot(time_1sec, PC_upper, label='AmBe upper')
#    pl.plot(time_1sec, PC_lower, label='AmBe lower')
#    pl.legend()
#    pl.show()

    
def make_grand_df(folder_header, chunk_number):
    list_of_dfs = []
    for i in range(chunk_number):
        
        directory = folder_header + str(i) + '/'
    #    directory = 'Results/data_run42_dbz1/test/'
    #    with open(directory + 'oneAMP_pulse_processing_results.p', 'rb') as fp1:
    #        oneAmp_processed_results = pickle.load(fp1)
        with open(directory + 'full_pulse_processing_results.p', 'rb') as fp1:
    #        (oneAmp_processed_results, twoAmp_processed_results, spike_processed_results, newT_oneamp, newT_twoamp) = pickle.load(fp1)
            (oneAmp_processed_results, Fe_processed_results, Ba_processed_results, 
             Neutron_processed_results) = pickle.load(fp1)
            
    #        (oneAmp_processed_results, twoAmp_processed_results, newT_oneamp , spike_processed_results, bq) = pickle.load(fp1)
#        (oneAmp_processed_results, twoAmp_processed_results, spike_processed_results) = pickle.load(fp1)
#        (oneAmp_processed_results, twoAmp_processed_results) = pickle.load(fp1)
#    with open(directory + 'NoisePSD.p', 'rb') as fp2:
#        (freq, J) = pickle.load(fp2)
        df_oneamp = pd.DataFrame(oneAmp_processed_results)
        df_Fe_temp = pd.DataFrame(Fe_processed_results)
        df_Ba_temp = pd.DataFrame(Ba_processed_results)
        df_Neutron_temp = pd.DataFrame(Neutron_processed_results)
        
        df_oneamp = df_oneamp.assign(Fe_amp=pd.Series(df_Fe_temp["amplitude"]).values)
        df_oneamp = df_oneamp.assign(Fe_chi=pd.Series(df_Fe_temp["chi"]).values)
        
        df_oneamp = df_oneamp.assign(Ba_amp=pd.Series(df_Ba_temp["amplitude"]).values)
        df_oneamp = df_oneamp.assign(Ba_chi=pd.Series(df_Ba_temp["chi"]).values)
            
        df_oneamp = df_oneamp.assign(AmBe_amp=pd.Series(df_Neutron_temp["amplitude"]).values)
        df_oneamp = df_oneamp.assign(AmBe_chi=pd.Series(df_Neutron_temp["chi"]).values)
        
        df_oneamp = df_oneamp.assign(chunk=pd.Series( [(i+1)]*len(df_oneamp.index) ).values)

        list_of_dfs.append(df_oneamp)
        
    result = pd.concat(list_of_dfs)
    result["delta_chi_Ba_AmBe"] = result.apply (lambda row: row.Ba_chi - row.AmBe_chi, axis=1)
    result["delta_chi_Ba_Fe"] = result.apply (lambda row: row.Ba_chi - row.Fe_chi, axis=1)
    result["delta_chi_Fe_AmBe"] = result.apply (lambda row: row.Fe_chi - row.AmBe_chi, axis=1)

    return result

def PCA_pulselist(list_of_pulses, list_ofNoisePSDs = []):
    obs_fouriers = []
    list_of_pulses = list(list_of_pulses)
    for r in range(len(list_of_pulses)):
        pulse = (list_of_pulses[r])
        pulse = highPassFilter(pulse, inverse=True)
#        obs_fourier = np.fft.rfft( pulse )
        f, obs_fourier =  getWindowedFFT(pulse)

        
        #            
        obs_fouriers.append(obs_fourier)
         
    if len(list_of_pulses) == 0:
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
    
    PCs = []
    for i in range(6):
        standard_template = theta_exp(one_second_interval, (0.24, 0.008,1.1))
        standard_template_max = max(standard_template)
        
        vector = (np.real(np.fft.irfft(vh[i])))
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
                avg_pulse = (np.fft.irfft(mean_fourier))
                avg_pulse = avg_pulse - avg_pulse[0]
                avg_pulse = max(scaled_vector)/max(avg_pulse)*avg_pulse
                pl.plot(one_second_interval, avg_pulse, label='average in freq domain')
#                pl.plot(one_second_interval, AB_theta_exp(one_second_interval, 0.24, 0.008,1.1, 2.5, 0.5, 0.5), label='Initial Guess 2 Amp')
#                pl.plot(one_second_interval, AB_theta_exp(one_second_interval, fit2[0][0], fit2[0][1], fit2[0][2], fit2[0][3], fit2[0][4], fit2[0][5])  , label='Best Fit 2 Amp', color='r')
                pl.legend()
                pl.show()
                
                pl.plot(one_second_interval, scaled_vector - avg_pulse, label = 'Principal Component 0 - avg', color='b')
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
        PCs.append(np.array(scaled_vector))
    
#    pl.plot(one_second_interval, PC[0]) + )
    lambdas = s**2
    print "Percent of variation explained by Principal Components"
    print lambdas[0:20]/(sum(lambdas[0:20]))
    return save_principal_0
    

file_name_Fe = 'data_run42_dbz1/20180315_14h24'; source = "Fe"; #No source#
file_name_Ba = 'data_run42_dbz1/20180313_18h32'; source = "Barium"; #good Ba 14 chunks run
file_name_Nu = 'data_run42_dbz1/20180314_16h12'; source = "Neutron";  #Neutrons 16 chunks
#source = "Neutron" #or "Barium", Neutron, or Fe
chunk_size = 1
chunk_number = 14
if (source == "Barium"): chunk_number = 14 #total number of chunks
elif (source == "Neutron"): chunk_number = 16
#chunk_number = 1

#folder_header = 'Results/' + file_name + '/' +str(chunk_size) + 'hours_'  #original runs
#folder_header = 'Results/' + 'standardized_starts/'+ file_name + '/' +str(chunk_size) + 'hours_'  #runs w standardized pulse start times
#folder_header = 'Results/' + 'test0/'+ file_name + '/' +str(chunk_size) + 'hours_'  #runs w standardized pulse start times
#folder_header = 'Results/' + 'cleancode_test2/'+ file_name + '/' +str(chunk_size) + 'hours_'  #runs w standardized pulse start times

#folder_header = 'Results/' + 'with_spike2/'+ file_name + '/' +str(chunk_size) + 'hours_'  #runs w standardized pulse start times


a1 = 0.02 # for Fe cut parabola
c1 = 800 # 
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

#list_BaNu_oneAmpResults = []
#list_Fe_oneAmpResults = []

list_Ba_oneAmpResults = []
list_Neutron_oneAmpResults = []
list_Fe_oneAmpResults = []





#fe_header = 'Results/' + 'cleancode_test2/'+ file_name_Fe + '/' +str(chunk_size) + 'hours_'  #runs w standardized pulse start times
#Fe_dataframe = make_grand_df(fe_header, 14)        
#Nu_header = 'Results/' + 'cleancode_test2/'+ file_name_Nu + '/' +str(chunk_size) + 'hours_'  #runs w standardized pulse start times
#Neutron_dataframe = make_grand_df(Nu_header, 16)   
#Ba_header = 'Results/' + 'cleancode_test2/'+ file_name_Ba + '/' +str(chunk_size) + 'hours_'  #runs w standardized pulse start times
#Ba_dataframe = make_grand_df(Ba_header, 14)   
#with open("Results/results_dataframes.p", 'wb') as fp:
#        pickle.dump( (Fe_dataframe, Ba_dataframe, Neutron_dataframe) , fp )
with open("Results/results_dataframes.p", 'rb') as fp2:
        (Fe_dataframe, Ba_dataframe, Neutron_dataframe) = pickle.load(fp2)
compare_results_new(Fe_dataframe, Ba_dataframe, Neutron_dataframe)




###list_freqs = []
#for i in range(chunk_number):
#    directory = folder_header + str(i) + '/'
##    directory = 'Results/data_run42_dbz1/test/'
##    with open(directory + 'oneAMP_pulse_processing_results.p', 'rb') as fp1:
##        oneAmp_processed_results = pickle.load(fp1)
#    with open(directory + 'full_pulse_processing_results.p', 'rb') as fp1:
##        (oneAmp_processed_results, twoAmp_processed_results, spike_processed_results, newT_oneamp, newT_twoamp) = pickle.load(fp1)
#        (oneAmp_processed_results, Fe_processed_results, Ba_processed_results, 
#         Neutron_processed_results) = pickle.load(fp1)
#        
##        (oneAmp_processed_results, twoAmp_processed_results, newT_oneamp , spike_processed_results, bq) = pickle.load(fp1)
##        (oneAmp_processed_results, twoAmp_processed_results, spike_processed_results) = pickle.load(fp1)
##        (oneAmp_processed_results, twoAmp_processed_results) = pickle.load(fp1)
#    with open(directory + 'NoisePSD.p', 'rb') as fp2:
#        (freq, J) = pickle.load(fp2)
#     
#        
#    
#    #Neutron Cuts    
#    #signal band   
#    if source=="Neutron":
#        oneAmp_posAmps = getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf)  ])
#        oneAmp_signal = generalizedCut(oneAmp_posAmps, "amplitude", "chi", cut_parab, (a1,c1) , True)
#        #pile-up events?
#        oneAmp_PU= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 2000, np.inf)  ])
#        oneAmp_PU= generalizedCut(oneAmp_PU, "amplitude", "chi", loglog_line, n0_PU_points_below, veto_greater_than_func= False )
#        oneAmp_PU= generalizedCut(oneAmp_PU, "amplitude", "chi", loglog_line, n0_PU_points_above, veto_greater_than_func= True)
#        #Spike band events?
#    #    print len(oneAmp_processed_results)
#        oneAmp_spike= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 2000, np.inf)  ])
#        oneAmp_spike= generalizedCut(oneAmp_spike, "amplitude", "chi", loglog_line, n0_spike_points_below, veto_greater_than_func= False )
#        oneAmp_spike= generalizedCut(oneAmp_spike, "amplitude", "chi", loglog_line, n0_spike_points_above, veto_greater_than_func= True)
#    
#    #Ba Cuts    
#    #signal band   
#    if source=="Barium":
#        oneAmp_posAmps = getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf)  ])
#        oneAmp_signal = generalizedCut(oneAmp_posAmps, "amplitude", "chi", cut_parab, (a1,c1) , True)
#        #pile-up events?
#        oneAmp_PU= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 2000, np.inf)  ])
#        oneAmp_PU= generalizedCut(oneAmp_PU, "amplitude", "chi", loglog_line, ba_PU_points_above, veto_greater_than_func= True )
#        oneAmp_PU = generalizedCut(oneAmp_PU, "amplitude", "chi", cut_parab, (a1,c1) , False)
#        #Spike band events?
#        oneAmp_spike= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 2000, np.inf)  ])
#        oneAmp_spike= generalizedCut(oneAmp_spike, "amplitude", "chi", loglog_line, ba_spike_points_below, veto_greater_than_func= False )
#        #blob
#        oneAmp_blob = generalizedCut(oneAmp_posAmps, "amplitude", "chi", cut_parab, (a1,c1) , False)
#        oneAmp_blob = getPulseSubset(oneAmp_blob, [("amplitude", 7, 10.5), ("chi", 800, 2000)  ])
#    if source=="Fe":
#        oneAmp_posAmps = getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf)  ])
#        oneAmp_signal = generalizedCut(oneAmp_posAmps, "amplitude", "chi", cut_parab, (a2,c2) , True)
#        #pile-up events?
#        oneAmp_PU= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 2000, np.inf)  ])
#        oneAmp_PU= generalizedCut(oneAmp_PU, "amplitude", "chi", loglog_line, fe_PU_points_above, veto_greater_than_func= True )
#        oneAmp_PU = generalizedCut(oneAmp_PU, "amplitude", "chi", cut_parab, (a2,c2) , False)
#        #Spike band events?
#        oneAmp_spike= getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf), ("chi", 10**4, np.inf)  ])
#        oneAmp_spike= generalizedCut(oneAmp_spike, "amplitude", "chi", loglog_line, fe_spike_points_below, veto_greater_than_func= False )
#
#
#    list_oneAmpResults.append(oneAmp_processed_results)
#    list_oneAmpResults_posAmps.append(oneAmp_posAmps)
#    list_oneAmpResults_signal.append(oneAmp_signal)
#    list_oneAmpResults_PU.append(oneAmp_PU)
#    list_oneAmpResults_spike.append(oneAmp_spike)
#    if source=="Barium":
#        list_oneAmpResults_blob.append(oneAmp_blob)
#    list_NoisePSDs.append( J )
#    
##    list_twoAmpResults.append(twoAmp_processed_results)
##    list_spikeResults.append(spike_processed_results)
#    
##    list_BaNu_oneAmpResults.append(BaNu_processed_results)
##    list_Fe_oneAmpResults.append(Fe_processed_results)
#    
#    list_Ba_oneAmpResults.append(Ba_processed_results)
#    list_Neutron_oneAmpResults.append(Neutron_processed_results)
#    list_Fe_oneAmpResults.append(Fe_processed_results)
#
##    list_new_oneAmpResults.append(newT_oneamp)
##    list_new_twoAmpResults.append(newT_twoamp) 
#
##for i in range(len(list_NoisePSDs)):
##    pl.loglog(freq, list_NoisePSDs[i], label="Noise PSD " + str(i))
###pl.legend()
##pl.xlabel("Frequency [Hz]")
##pl.ylabel("Noise PSD")
##pl.ylim([10**-6,10**2])
##pl.xlim([1,500])
##pl.show()    
#



#Tr_fit = 0.008
#Tf_fit = 1.1
#print "====================================================="
#print "REFIT TEST"
#print "====================================================="
##reFit_optimalFiltering1amp(list_oneAmpResults, list_NoisePSDs )
##assert False

#print "====================================================="
#print "ALL POSITIVE AMPLITUDE PULSES"
#print "====================================================="
#analyzeFitResults(list_oneAmpResults_posAmps)
###compareFitResults(list_oneAmpResults, list_twoAmpResults, list_spikeResults)
#compareFitResults(list_oneAmpResults, list_new_twoAmpResults, list_spikeResults)
#PCAtemplate(list_oneAmpResults_posAmps, list_NoisePSDs)
##spectralClusterPulses(list_oneAmpResults_posAmps)
#print "====================================================="
#print "SIGNAL BAND PULSES ONLY"
#print "====================================================="
#analyzeFitResults(list_oneAmpResults)
#analyzeFitResults(list_Fe_oneAmpResults)
#analyzeFitResults(list_BaNu_oneAmpResults)
#compareFitResults(list_oneAmpResults, list_Fe_oneAmpResults, list_Ba_oneAmpResults)
#compareFitResults(list_oneAmpResults, list_Neutron_oneAmpResults, list_Ba_oneAmpResults)
#assert False
#pc0_signal, pc0_signal_weighted = PCAtemplate(list_oneAmpResults_signal, list_NoisePSDs)

#print "====================================================="
#print "PC Fitting "
#print "====================================================="
#with open('Fe_pc0s.p', 'rb') as fp:
#        (Fe_pc0_signal, Fe_pc0_spike) = pickle.load(fp)
#  
#one_sec = np.linspace(0,1,1000)      
#pl.plot(one_sec, Fe_pc0_signal, label='Fe PC0')
#template = ratio_theta_exp
##initial_guess = (1.0, 0.24, 6.777e-3, 1.049, 20e-3, 30e-3, 0.3 )
#initial_guess = [0.60891902, 0.24179725, 0.0043958,  2.0856438,  0.01805898, 0.42211094, 0.71273386]
#fit2 =  curve_fit(template, one_sec, Fe_pc0_signal, initial_guess)
#pl.plot(one_sec, template(one_sec, *initial_guess ), label = 'Initial Guess')
#pl.plot(one_sec, template(one_sec, *fit2[0]), 'r', label = 'Best Fit')
#pl.plot(one_sec, A_theta_exp(one_sec, 1.0,0.24,8e-3, 1.1 ), label = 'OG Temp')
#pl.legend()
#pl.show()
#print fit2[0]
#
#with open('Ba_pc0s.p', 'rb') as fp:
#        (Ba_pc0_signal, Ba_pc0_spike) = pickle.load(fp)
#  
#pl.plot(one_sec, Ba_pc0_signal, label='Ba PC0')
#template = ratio_theta_exp
#initial_guess = fit2[0]
#fit3 =  curve_fit(template, one_sec, Ba_pc0_signal, initial_guess)
#pl.plot(one_sec, template(one_sec, *initial_guess ), label = 'Initial Guess ')
#pl.plot(one_sec, template(one_sec, *fit3[0]), 'r', label = 'Best Fit')
#pl.legend()
#pl.show()
#print fit3[0]
#
#with open('Nu_pc0s.p', 'rb') as fp:
#        (Nu_pc0_signal, Nu_pc0_spike) = pickle.load(fp)
#  
#pl.plot(one_sec, Nu_pc0_signal, label='Neutron PC0')
#template = ratio_theta_exp
#initial_guess = fit2[0]
#fit3 =  curve_fit(template, one_sec, Nu_pc0_signal, initial_guess)
#pl.plot(one_sec, template(one_sec, *initial_guess ), label = 'Initial Guess ')
#pl.plot(one_sec, template(one_sec, *fit3[0]), 'r', label = 'Best Fit')
#pl.legend()
#pl.show()
#print fit3[0]
#
#pl.plot(one_sec, Fe_pc0_signal, label="Fe PC0")
#pl.plot(one_sec, Ba_pc0_signal, label="Ba PC0")
#pl.plot(one_sec, Nu_pc0_signal, label="Nu PC0")
##pl.plot(one_sec, template(one_sec, *fit2[0]), 'r', label = 'Best Fit')
#pl.legend()
#pl.show()
#
#
#pl.semilogy(one_sec, np.abs(Fe_pc0_signal-Ba_pc0_signal), label="Fe PC0-Ba")
#pl.legend()
#pl.show()
#
#pl.semilogy(one_sec, np.abs(Fe_pc0_signal-Nu_pc0_signal), label="Fe PC0-Nu")
#pl.legend()
#pl.show()
#
#pl.semilogy(one_sec, np.abs(Ba_pc0_signal-Nu_pc0_signal), label="Ba PC0-Nu")
#pl.legend()
#pl.show()

#newFitComparison(list_oneAmpResults_signal, list_NoisePSDs, Fe_pc0_signal)
#newFitComparison(list_oneAmpResults_signal, list_NoisePSDs, Fe_pc0_signal)

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

#Compare sources
#with open('Barium_deltachiBa_Fe.p', 'rb') as fp:
#        (sigamps_BaFe, deltaBa_BaFe) = pickle.load(fp)
#with open('Fe_deltachiBa_Fe.p', 'rb') as fp:
#        (sigamps2_BaFe, deltaFe_BaFe) = pickle.load(fp)
#with open('Neutron_deltachiBa_Fe.p', 'rb') as fp:
#        (sigamps3_BaFe, deltaNu_BaFe) = pickle.load(fp)   
#
#with open('Barium_deltachiBa_Nu.p', 'rb') as fp:
#        (sigamps, deltaBa) = pickle.load(fp)
#with open('Fe_deltachiBa_Nu.p', 'rb') as fp:
#        (sigamps2, deltaFe) = pickle.load(fp)
#with open('Neutron_deltachiBa_Nu.p', 'rb') as fp:
#        (sigamps3, deltaNu) = pickle.load(fp)   
##
#fig5 = pl.figure()
#ax = pl.gca()
##ax.plot(sigamps , deltaBa  , '.', c='black', label='Ba Data', markeredgecolor='none')
#ax.plot(sigamps3 , deltaNu  , '.', c='red', label='AmBe Data', markeredgecolor='none')
#ax.plot(sigamps , deltaBa  , '.', c='black', label='Ba Data', markeredgecolor='none')
#ax.plot(sigamps2 , deltaFe  , '.', c='blue', label='Fe Data', markeredgecolor='none')
#ax.set_xlabel('Amplitude of Original Template Fit')
#ax.set_ylabel('Delta Chi Sq. (Ba Template - Neutron Template)')
#ax.minorticks_on()
#ax.set_ylim([-0.1*10**3, 0.1*10**3])
#ax.set_xlim([0, 200])
#pl.legend()
##pl.savefig('deltaChiSq_newTemplates.png', dpi=1200)
#pl.show()
#fig5 = pl.figure()
#ax = pl.gca()
##ax.plot(sigamps , deltaBa  , '.', c='black', label='Ba Data', markeredgecolor='none')
#ax.plot(sigamps3 , deltaNu  , '.', c='red', label='AmBe Data', markeredgecolor='none')
##ax.plot(sigamps , deltaBa  , '.', c='black', label='Ba Data', markeredgecolor='none')
#ax.plot(sigamps2 , deltaFe  , '.', c='blue', label='Fe Data', markeredgecolor='none')
#ax.set_xlabel('Amplitude of Original Template Fit')
#ax.set_ylabel('Delta Chi Sq. (Ba Template - Neutron Template)')
#ax.minorticks_on()
#ax.set_ylim([-0.1*10**3, 0.1*10**3])
#ax.set_xlim([0, 200])
#pl.legend()
##pl.savefig('deltaChiSq_newTemplates.png', dpi=1200)
#pl.show()
#
#fig5 = pl.figure()
#ax = pl.gca()
##ax.plot(sigamps , deltaBa  , '.', c='black', label='Ba Data', markeredgecolor='none')
#ax.plot(sigamps3 , deltaNu  , '.', c='red', label='AmBe Data', markeredgecolor='none')
#ax.plot(sigamps , deltaBa  , '.', c='black', label='Ba Data', markeredgecolor='none')
##ax.plot(sigamps2 , deltaFe  , '.', c='blue', label='Fe Data', markeredgecolor='none')
#ax.set_xlabel('Amplitude of Original Template Fit')
#ax.set_ylabel('Delta Chi Sq. (Ba Template - Neutron Template)')
#ax.minorticks_on()
#ax.set_ylim([-0.1*10**3, 0.1*10**3])
#ax.set_xlim([0, 200])
#pl.legend()
##pl.savefig('deltaChiSq_newTemplates.png', dpi=1200)
#pl.show()

#fig5 = pl.figure()
#ax = pl.gca()
#ax.plot(deltaBa_BaFe , deltaBa  , '.', c='black', label='Ba Data', markeredgecolor='none')
#ax.plot(deltaNu_BaFe , deltaNu  , '.', c='red', label='AmBe Data', markeredgecolor='none')
#ax.plot(deltaFe_BaFe , deltaFe  , '.', c='blue', label='Fe Data', markeredgecolor='none')
#ax.set_xlabel('Delta Chi Sq. (Ba Template - Fe Template)')
#ax.set_ylabel('Delta Chi Sq. (Ba Template - Neutron Template)')
#ax.minorticks_on()
#ax.set_ylim([-250, 200])
#ax.set_xlim([-4e3, 4e3])
#pl.legend()
##pl.savefig('deltaChiSq_newTemplates.png', dpi=1200)
#pl.show()
#
#fig5 = pl.figure()
#ax = pl.gca()
##ax.plot(deltaBa_BaFe , deltaBa  , '.', c='black', label='Ba Data', markeredgecolor='none')
##ax.plot(deltaNu_BaFe , deltaNu  , '.', c='red', label='AmBe Data', markeredgecolor='none')
#ax.plot(deltaFe_BaFe , deltaFe  , '.', c='blue', label='Fe Data', markeredgecolor='none')
#ax.set_xlabel('Delta Chi Sq. (Ba Template - Fe Template)')
#ax.set_ylabel('Delta Chi Sq. (Ba Template - Neutron Template)')
#ax.minorticks_on()
#ax.set_ylim([-250, 200])
#ax.set_xlim([-2e3, 2e3])
#pl.legend()
##pl.savefig('deltaChiSq_newTemplates.png', dpi=1200)
#pl.show()

#Compare sources
#with open('Barium_amp_spectrum.p', 'rb') as fp:
#        (Barium_hist_1, Barium_hist_2) = pickle.load(fp)
#with open('Neutron_amp_spectrum.p', 'rb') as fp:
#        (Neutron_hist_1, Neutron_hist_2) = pickle.load(fp)   
#with open('Fe_amp_spectrum.p', 'rb') as fp:
#        (Fe_hist_1, Fe_hist_2) = pickle.load(fp)
##        
#        
#pl.errorbar(np.linspace(0.5,200,num=100), Neutron_hist_2[0] - (16/14.0)*Fe_hist_2[0], np.sqrt(Neutron_hist_2[0] + ((16/14.0)*Fe_hist_2[0] )), fmt=''  ,  label = 'Neutron - Fe', color = 'red')        
#pl.errorbar(np.linspace(0.5,200,num=100), Neutron_hist_2[0], np.sqrt(Neutron_hist_2[0]), fmt='', label = 'Neutron', color = 'blue')        
#pl.legend()
#pl.show()
#
#pl.errorbar(np.linspace(0.5,200,num=100), Barium_hist_2[0] - Fe_hist_2[0], np.sqrt(Barium_hist_2[0] + Fe_hist_2[0] ) , fmt='', label = 'Barium - Fe', color = 'red')        
#pl.errorbar(np.linspace(0.5,200,num=100), Barium_hist_2[0], np.sqrt(Barium_hist_2[0]), fmt='', label = 'Barium', color = 'blue')        
#pl.legend()
#pl.show()
#
#pl.errorbar(np.linspace(0.5,200,num=200), Neutron_hist_1[0] - (16/14.0)*Fe_hist_1[0], np.sqrt(Neutron_hist_1[0] + ((16/14.0)*Fe_hist_1[0] )), fmt=''  ,  label = 'Neutron - Fe', color = 'red')        
#pl.errorbar(np.linspace(0.5,200,num=200), Neutron_hist_1[0], np.sqrt(Neutron_hist_1[0]), fmt='', label = 'Neutron', color = 'blue')        
#pl.legend()
#pl.show()
#
#pl.errorbar(np.linspace(0.5,200,num=200), Barium_hist_1[0] - Fe_hist_1[0], np.sqrt(Barium_hist_1[0] + Fe_hist_1[0] ) , fmt='', label = 'Barium - Fe', color = 'red')        
#pl.errorbar(np.linspace(0.5,200,num=200), Barium_hist_1[0], np.sqrt(Barium_hist_1[0]), fmt='', label = 'Barium', color = 'blue')        
#pl.legend()
#pl.show()

#pl.plot(np.linspace(0.5,200,num=200), Neutron_hist_1[0], 'v', label = 'Neutron', color = 'red', markeredgecolor='none')        
#pl.plot(np.linspace(0.5,200,num=200), Fe_hist_1[0], '^', label = 'Fe', color = 'black', markeredgecolor='none')        
#pl.plot(np.linspace(0.5,200,num=200), Barium_hist_1[0], 'x', label = 'Barium', color = 'blue')        
#pl.ylabel("Events")
#pl.xlabel("Amplitude [ADU]")
#pl.legend()
#pl.show()
##
#pl.plot(np.linspace(0.5,200,num=100), Neutron_hist_2[0], 'v', label = 'Neutron', color = 'red', markeredgecolor='none')        
#pl.plot(np.linspace(0.5,200,num=100), Fe_hist_2[0], '^', label = 'Fe', color = 'black', markeredgecolor='none')        
#pl.plot(np.linspace(0.5,200,num=100), Barium_hist_2[0], 'x', label = 'Barium', color = 'blue')        
#pl.ylabel("Events")
#pl.xlabel("Amplitude [ADU]")
#pl.legend()
#pl.show()
##
#pl.plot(np.linspace(0.5,200,num=200), Neutron_hist_1[0]/(sum(Neutron_hist_1[0])), 'v', label = 'Neutron', color = 'red', markeredgecolor='none')        
#pl.plot(np.linspace(0.5,200,num=200), Fe_hist_1[0]/(sum(Fe_hist_1[0])), '^', label = 'Fe', color = 'black', markeredgecolor='none')        
#pl.plot(np.linspace(0.5,200,num=200), Barium_hist_1[0]/(sum(Barium_hist_1[0])), 'x', label = 'Barium', color = 'blue')        
#pl.ylabel("Fraction of Total Events")
#pl.xlabel("Amplitude [ADU]")
#pl.legend()
#pl.show()
##
#pl.plot(np.linspace(0.5,200,num=100), Neutron_hist_2[0]/(sum(Neutron_hist_2[0])), 'v', label = 'Neutron', color = 'red', markeredgecolor='none')        
#pl.plot(np.linspace(0.5,200,num=100), Fe_hist_2[0]/(sum(Fe_hist_2[0])), '^', label = 'Fe', color = 'black', markeredgecolor='none')        
#pl.plot(np.linspace(0.5,200,num=100), Barium_hist_2[0]/(sum(Barium_hist_2[0])), 'x', label = 'Barium', color = 'blue')        
#pl.ylabel("Fraction of Total Events")
#pl.xlabel("Amplitude [ADU]")
#pl.legend()
#pl.show()





#        
#with open('Barium_ampratio.p', 'rb') as fp:
#        (Barium_normal_spikeamp, Barium_normal_slowspikeratio) = pickle.load(fp)
#with open('Neutron_ampratio.p', 'rb') as fp:
#        (Neutron_normal_spikeamp, Neutron_normal_slowspikeratio) = pickle.load(fp) 
#with open('Fe_ampratio.p', 'rb') as fp:
#        (Fe_normal_spikeamp, Fe_normal_slowspikeratio) = pickle.load(fp)
#        
#fig8 = pl.figure()
#ax = pl.gca()
#ax.plot(Neutron_normal_spikeamp , Neutron_normal_slowspikeratio  , '.', c='red', markeredgecolor='red', fillstyle='none', label='AmBe')
#ax.plot(Barium_normal_spikeamp, Barium_normal_slowspikeratio, '.', c='black', markeredgecolor='black', fillstyle='none', label = 'Ba')
#ax.plot(Fe_normal_spikeamp ,  Fe_normal_slowspikeratio, '.', c='blue', markeredgecolor='blue', fillstyle='none', label='Fe')
#ax.set_yscale('log')
##   ax.set_xscale('log')
#ax.set_ylabel('Slow Amp/Spike Amp in Two Template Fit')
#ax.set_xlabel('Spike Amplitude in Two Template Fit')
#ax.minorticks_on()
#ax.set_ylim([10**-5, 10**3])
#ax.set_xlim([0, 5*10**2])
#ax.legend()
#pl.savefig('AmplitudeRatio_withFe.png')
#pl.show()
#        
#fig8 = pl.figure(figsize=(20,10))
#pl.rcParams.update({'font.size': 24})
#ax = pl.gca()
#ax.plot(Neutron_normal_spikeamp , Neutron_normal_slowspikeratio  , '.', c='red', markeredgecolor='red', fillstyle='none', label='AmBe')
#ax.plot(Barium_normal_spikeamp, Barium_normal_slowspikeratio, '.', c='black', markeredgecolor='black', fillstyle='none', label = 'Ba')
#ax.plot(Fe_normal_spikeamp ,  Fe_normal_slowspikeratio, '.', c='blue', markeredgecolor='blue', fillstyle='none', label='Fe')
#ax.set_yscale('log')
##   ax.set_xscale('log')
#ax.set_ylabel('Slow Amp/Spike Amp in Two Template Fit')
#ax.set_xlabel('Spike Amplitude in Two Template Fit')
#ax.minorticks_on()
#ax.set_ylim([10**-5, 10**3])
#ax.set_xlim([0, 5*10**2])
#ax.legend()
#pl.savefig('AmplitudeRatio_withFe_big.png')
#pl.show()



#one_sec = np.linspace(0, 1, int(fe))
        
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



