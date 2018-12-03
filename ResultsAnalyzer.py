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
    DoF=len(getWindowedFFT(templates[0])[1])

    amps=np.array(amps)        
    chis_red=np.array(chis)/DoF
        
    total=0
    for i in list_ofPulse_Dict_Lists:
        print "Number of Pulse Fit is ", len(i)
        total += len(i)

    print "Number of Total Pulses", total

#    pl.hist(chis_red)
#    pl.xlabel('Chi Sq (Reduced)')
#    pl.ylabel('Events')
#    pl.show()
    
#    pl.hist(chis, range=(4e2,6e2))
#    pl.xlabel('Chi Sq')
#    pl.ylabel('Events')
#    pl.show()

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


    amp_hist = pl.hist(amps, range=(0,200), bins = 100)
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
#    pl.hist(amps*(1.1 - 0.008), range=(0,200), bins=100)
#    pl.xlabel('Integrated Pulse [A.D.U. * s]')
#    pl.ylabel('Events')
#    pl.xlim([0,200])
#    pl.show()
    
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
    
    pl.scatter(np.abs(amps), chis_red)
    pl.xlabel("Amplitude [A.D.U.]")
    pl.ylabel("Reduced X-square")
    pl.xlim([0, 200])
    pl.ylim([0,5])
    pl.show()
    
    pl.scatter(np.abs(amps), chis)
    pl.xlabel("Amplitude [A.D.U.]")
    pl.ylabel("X-Square")
    pl.xlim([0, 200])
    pl.ylim([0,5*500])
    pl.show()
    
    fig = pl.figure()
    ax = pl.gca()
    ax.plot(amps , chis, '.', c='black', markeredgecolor='none')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Chi Sq.')
    x_plot = np.linspace(0, 10*(max(amps)), 250 )
    ax.plot(x_plot, cut_parab(x_plot, (a1, c1)), color='red')
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
        
        pl.scatter(slow_amps, fast_amps)
        pl.xlabel('Slow Amplitude')
        pl.ylabel('Fast Amplitude')
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
    
    else:
        pl.scatter(amps, chis_red)
        pl.xlabel('Amplitude')
        pl.ylabel('Chi Sq (Reduced)')
        pl.ylim([0.0,4])
        pl.show()
        
        pl.scatter(amps, chis) 
        pl.xlabel('Amplitude')
        pl.ylabel('Chi Sq')
        pl.ylim([-5,200e3])
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

def theta_exp(t, (t_0,Tr,Tf)):
    
    offset=[]
    offset.append(t_0)
    z=len(t)*offset
    z=np.array(z)
    m=expo_temp1(t-z, (Tr,Tf))
        
    return step(t-z)*m
   
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
            
            
            
            
def compareFitResults(list_of_pulse_dict_list1, list_of_pulse_dict_list2 ):
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
    
def PCAtemplate(list_ofPulse_Dict_Lists):
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
            
    pulse_matrix = np.asarray(obs_pulses)
#    mean_pulse = np.mean(obs_pulses, axis=0) 
#    centered_pulse_matrix = pulse_matrix - mean_pulse
    svd_output = np.linalg.svd(pulse_matrix)
    
    print "\n PCA Time Domain Uncentered\n"
    
    (u, s, vh) = svd_output
    print "lengths", len(vh[0]), len(u[:,0])
    print np.shape(u), np.shape(vh)
    one_second_interval = np.linspace(0,1, int(fe))
    for i in range(6):
        standard_template = theta_exp(one_second_interval, (0.24, 0.008,1.1))
        standard_template_max = max(standard_template)
        
        zeroed_vector = vh[i] - vh[i][0]
        zeroed_vector[0:int(0.225*fe)] = 0
        if max(np.abs(zeroed_vector)) == max(zeroed_vector):
            scaling = standard_template_max/max(zeroed_vector)
        else:
            scaling = standard_template_max/min(zeroed_vector)
        scaled_vector = scaling*zeroed_vector
        if i==0:
            fit = curve_fit(A_theta_exp, one_second_interval, scaled_vector, (1.0, 0.24, 0.008,1.1))
            print fit
            pl.plot(one_second_interval, A_theta_exp(one_second_interval, 1.0, 0.24, 0.008,1.1), label='Initial Guess')
            pl.plot(one_second_interval, scaled_vector, label = 'Principal Component 0', color='b')
            pl.plot(one_second_interval, fit[0][0]*theta_exp(one_second_interval, fit[0][1:]) , label='Best Fit', color='r')
            pl.legend()
            pl.show()
            print "Best Params", fit[0]
            print "Best Time Constants", fit[0][2:]
            print "Cov of Rise Time", fit[1][2]
            print "Cov of Fall Time", fit[1][3]
            print "Cov Matrix", fit[1]
            
        pl.plot(one_second_interval, vh[i], label= "Singular Vector " + str(i) , color = 'b')
        pl.plot(one_second_interval, standard_template , label='Template', color = 'r')
        pl.xlabel('Time [s]')
        pl.ylabel('Amplitude [ADU]')
        pl.legend()
        pl.show()
        
    lambdas = s**2
    print lambdas[0:20]
    print lambdas[0:20]/(sum(lambdas[0:20]))
        
#    print "\n PCA Time Domain Centered\n"

    

#    (u, s, vh) = np.linalg.svd(centered_pulse_matrix)
#    
#    print "lengths", len(vh[0]), len(u[:,0])
#    print np.shape(u), np.shape(vh)
#    one_second_interval = np.linspace(0,1, int(fe))
#    for i in range(6):
#        pl.plot(one_second_interval, -10*vh[i], label= "Singular Vector " + str(i) , color = 'b')
#        pl.plot(one_second_interval, theta_exp(one_second_interval, (0.035, 0.008,1.1)), label='Template', color = 'r')
#        pl.xlabel('Time [s]')
#        pl.ylabel('Amplitude [ADU]')
#        pl.legend()
#        pl.show()
#    
#    lambdas = s**2
#    print lambdas[0:20]
#    print lambdas[0:20]/(sum(lambdas[0:20]))
#    pl.plot(one_second_interval, mean_pulse)
#    pl.show()
#    
#    print "\n PCA Freq Domain Uncentered\n"
#    
#    obs_fouriers = []
#    for i in obs_pulses:
#        obs_fouriers.append(np.fft.fft(i))
#    print "hi1"
#    obs_fouriers = np.asarray(obs_fouriers)
#    (u, s, vh) = np.linalg.svd(obs_fouriers)
#    print "hi2"
#    print "lengths", len(vh[0]), len(u[:,0])
#    print np.shape(u), np.shape(vh)
#    one_second_interval = np.linspace(0,1, int(fe))
#    for i in range(6):
#        pl.plot(one_second_interval, np.fft.ifft(vh[i]), label= "Singular Vector " + str(i) , color = 'b')
#        pl.plot(one_second_interval, theta_exp(one_second_interval, (0.035, 0.008,1.1)), label='Template', color = 'r')
#        pl.xlabel('Time [s]')
#        pl.ylabel('Amplitude [ADU]')
#        pl.legend()
#        pl.show()
#        
#    pl.plot(one_second_interval, np.fft.ifft(obs_fouriers[0]))
#    pl.show()
#        
#    lambdas = s**2
#    print lambdas[0:20]
#    print lambdas[0:20]/(sum(lambdas[0:20]))
#    (u, s, vh) = np.linalg.svd(centered_pulse_matrix)
#    
#    print "\n PCA Freq Domain Centered\n"

    

    


    
    

#file_name = 'data_run42_dbz1\\20180315_14h24' #No source#
#file_name = 'data_run42_dbz1\\20180313_18h32' #good Ba 12 chunks
file_name = 'data_run42_dbz1\\20180314_16h12' #Neutrons 15
chunk_size = 1
chunk_number = 1 #total number of chunks
folder_header = 'Results/' + file_name + '/' +str(chunk_size) + 'hours_' 

a1 = 0.02
c1 = 800

list_oneAmpResults = []
list_oneAmpResults_afterCUT = []
list_oneAmpResults_afterCUT2= []

list_NoisePSDs = []
#list_freqs = []


for i in range(chunk_number):
#    directory = folder_header + str(i) + '/'
    directory = 'Results/data_run42_dbz1/test/'
    with open(directory + 'oneAMP_pulse_processing_results.p', 'rb') as fp1:
        oneAmp_processed_results = pickle.load(fp1)
#    with open(directory + 'NoisePSD.p', 'rb') as fp2:
#        (freq, J) = pickle.load(fp2)
        
    oneAmp_afterCUT = getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf)  ])
    oneAmp_afterCUT2 = generalizedCut(oneAmp_afterCUT, "amplitude", "chi", cut_parab, (a1,c1) , True)
#    
    
#    oneAmp_afterCUT = getPulseSubset(oneAmp_processed_results, [("amplitude", 0, np.inf)  ])
#    oneAmp_afterCUT = generalizedCut(oneAmp_afterCUT, "amplitude", "chi", line, (4985/14.0, -477500/7.0) , False)
#    oneAmp_afterCUT2 = generalizedCut(oneAmp_afterCUT, "amplitude", "chi", line, (19940/7.0, -976000) , True)
    
    list_oneAmpResults.append(oneAmp_processed_results)
    list_oneAmpResults_afterCUT.append(oneAmp_afterCUT)
    list_oneAmpResults_afterCUT2.append(oneAmp_afterCUT2)
#    list_NoisePSDs.append( (freq,J) )


#PCAtemplate(list_oneAmpResults)
#analyzeFitResults(list_oneAmpResults)
PCAtemplate(list_oneAmpResults_afterCUT)
print "====================================================="
print "====================================================="
#analyzeFitResults(list_oneAmpResults_afterCUT)
PCAtemplate(list_oneAmpResults_afterCUT2)
#analyzeFitResults(list_oneAmpResults_afterCUT2)
#analyzeFitResults(list_oneAmpResults)
#analyzeFitResults(list_oneAmpResults_afterCUT )
    

#with open('pulses/no_source315/second5h/full_pulse_processing_results_second5h.p', 'rb') as fp2:
#    (oneAmp_processed_results2, twoAmp_processed_results2) = pickle.load(fp2)
##with open('pulses/no_source315/remaining_hours/full_pulse_processing_results_last.p', 'rb') as fp3:
##    (oneAmp_processed_results3, twoAmp_processed_results3) = pickle.load(fp3)
#analyzeFitResults([oneAmp_processed_results1, oneAmp_processed_results2])
#analyzeFitResults([getPulseSubset(oneAmp_processed_results1, [("chi", 300,800), ("amplitude", 0, 1000) ])
#        , getPulseSubset(oneAmp_processed_results2, [("chi", 300,800), ("amplitude", 0, 1000)  ])])           
#assert False
