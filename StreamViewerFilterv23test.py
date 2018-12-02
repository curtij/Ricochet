#!/usr/bin/python

# coding:utf-8



####################################

## StreamSamba_trigger.py

##

## Analys des pulses

####################################



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
#pl.minorticks_on()
#pl.tick_params(axis='both', which='minor', length=2, width=1)
#pl.tick_params(axis='both', which='major', length=5, width=2)    
#


#np.random.choice([])


## Binary file





dx = 1
Gain = 1./4600
voltageRange = 20.
ADCRange = 2**16
Gain = Gain*voltageRange/ADCRange
fe = 1000.
f_Nyq=fe/2.
npoint = fe*dx
dt = int(1*3600) #change time length here #in s
#file_start_time = 2*5*3600+10 #in s
file_start_time = 0


##Fichier binaire



def chisq(V, exp, J):
    return sum((np.abs(V - exp ))**2 / J)  


def expdecayfn(t, (P,r)):
    return P*np.exp(-1*r*t)

def expgrowthfn(t, (P,r)):
    return P*np.exp(1*r*t)

def exprise(t, (f,i,r)):
    return f-(f-i)*np.exp(-1*r*t)

def linfn(t, (m,b)):
    return m*t+b
    
def doubleExp(t, (A, s, w)):
    return A*(np.exp(-t/s) - np.exp(-t/w))



def expo_temp1(t, (Tr,Tf)):
    
    return (-np.exp(-t/Tr)+np.exp(-t/Tf)) 

def theta_exp(t, (t_0,Tr,Tf)):
    
    offset=[]
    offset.append(t_0)
    z=len(t)*offset
    z=np.array(z)
    m=expo_temp1(t-z, (Tr,Tf))
    
    
    
    return step(t-z)*m

def step(x):
    return 1 * (x > 0) 

def A_theta_exp(t, (t_0, Tr, Tf, A, b)):
    shifted_t = np.asarray(t) - t_0    
    m=expo_temp1(shifted_t, (Tr,Tf))
    
    return A*step(shifted_t)*m + b

def AB_theta_exp(t, (t_0, Tr, Tf1, Tf2, A, B, b)):
    shifted_t=np.asarray(t)-t_0
    m=expo_temp1(shifted_t, (Tr, Tf1))
    n=expo_temp1(shifted_t, (Tr, Tf2))
    
    H=step(shifted_t)
    
    return H*(A*m+B*n)+b

def constant(t, (b,a)):
    m=t+b+a-a-t
    return m

def chisqmin(initial_parameter_guess, x_vals, y_vals, errors, fitFunction, meth = None):
    
    def chisqfunc( paramTuple ):
    #Evaluates the chi square value for a given gaussian model with parameters a
        model_values = fitFunction(x_vals, paramTuple)
        chisq = pylab.sum(((y_vals - model_values)/errors)**2)
        return chisq    
    if (meth != None):
        result = opt.minimize(chisqfunc, initial_parameter_guess, method = meth)
    else:
        result = opt.minimize(chisqfunc, initial_parameter_guess)
    print "\n\n"
    print result
#    assert result.success==True
#    print result.x
    #a,b = result.x
    #print "\nBest Fit Coefficient: " + str(a)
    #print "Best Fit Mean: " + str(b)
    
    DoF = len(y_vals) - len(initial_parameter_guess)
    chi = chisqfunc(initial_parameter_guess)
    #print "\nDegrees of Freedom: " + str(DoF)
    #print "Initial Guess's Chi square value: " + str(chi)
    #print "Initial Guess's Reduced Chi Square Value: " + str((chi/float(DoF)))
    #print "p-value " + str(1 - stats.chi2.cdf(chi, DoF))
    
    chisq_ofFit = chisqfunc(result.x)
    #print "\nDegrees of Freedom: " + str(DoF)
    #print "Best Fit's Chi square value: " + str(chisq_ofFit)
    #print "Best Fit's Reduced Chi Square Value: " + str((chisq_ofFit/float(DoF)))
    #print "p-value " + str(1 - stats.chi2.cdf(chisq_ofFit, DoF))
    #print "\n===================================================="
    #(a,b, c) =result.x    
    result_params =[]
    for t in result.x:
        result_params.append(t)
    
    #uncertainty calculations
    uncs = []
    unc_intervals = []
#    for index in range(len(result_params)):
#        copy_params = result_params[:]
#        min_deltaChi1 = 999
#        new_a = -5
#        a = result_params[index]
#        for i in pylab.linspace(a*0.5, a, 10000):
#            copy_params[index] = i
#            if ( abs(chisqfunc( tuple(copy_params) ) - chisq_ofFit - 1) < min_deltaChi1):
#                min_deltaChi1 = chisqfunc( tuple(copy_params) ) - chisq_ofFit - 1
#                new_a = i
#        copy_params[index] = new_a 
##        print "\nNew a such that chi varies by one " + str(new_a)
##        print "Lower Uncertainty on a: " + str(new_a - a)
##        print "Chi Square difference: " + str(chisqfunc( tuple(copy_params) ) - chisq_ofFit)
#        lower_unc = abs(new_a - a)
#        min_deltaChi1 = 999
#        new_a = -5
#        for i in pylab.linspace(a*1.5, a, 10000):
#            copy_params[index] = i
#            if ( abs(chisqfunc( tuple(copy_params) ) - chisq_ofFit - 1) < min_deltaChi1):
#                min_deltaChi1 = chisqfunc( tuple(copy_params) ) - chisq_ofFit - 1
#                new_a = i
#        copy_params[index] = new_a 
##        print "\nNew a such that chi varies by one " + str(new_a)
##        print "Upper Uncertainty on a: " + str(new_a - a)
##        print "Chi Square difference: " + str(chisqfunc( tuple(copy_params) ) - chisq_ofFit)
#        upper_unc = abs(new_a-a)
#        uncs.append((lower_unc+upper_unc)/2.0)
#        unc_intervals.append( (lower_unc, upper_unc) )
#        
    #print "Best Fit Parameters and Best fit uncs"
    #print result_params
    #print uncs
    #print "\n===================================================="
   
    return result.x, uncs, DoF, chisq_ofFit, unc_intervals
    
# Implementation of algorithm from http://stackoverflow.com/a/22640362/6029703
def thresholding_algo(y, lag, threshold, influence):
    print "Thresholding Algo begins at", str(datetime.now())
    #returns signals, avgfilter, stdfilter
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = np.asarray( ([0.0]*len(y)) )
    stdFilter = np.asarray( ([0.0]*len(y)) )
    avgFilter[lag - 1] = np.mean(y[0:lag])
    stdFilter[lag - 1] = np.std(y[0:lag])
    
#    avgFilterv2 = np.asarray( ([0.0]*len(y)) )
#    stdFilterv2 = np.asarray( ([0.0]*len(y)) )
    sum_running = sum(filteredY[0:lag])
    sumSQ_running=sum(filteredY[0:lag]**2)
    length = float(lag)

    new_mean = sum_running/length
    new_std = np.sqrt((sumSQ_running - 2*new_mean*sum_running + length*new_mean**2)/length)        

    for i in range(lag, len(y) - 1):
        if i > lag:
             sum_running = sum_running - filteredY[(i-lag-1)] + filteredY[i-1]
             sumSQ_running= sumSQ_running - filteredY[(i-lag-1)]**2 + filteredY[i-1]**2
        new_mean = sum_running/length
        new_std = np.sqrt((sumSQ_running - 2*new_mean*sum_running + length*new_mean**2)/length)


    
        if abs(y[i] - avgFilter[i-1]) > threshold * stdFilter [i-1]:
            if y[i] > avgFilter[i-1]:
                signals[i] = 1
            else:
                signals[i] = -1

            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
            avgFilter[i] = new_mean
            stdFilter[i] = new_std

        else:
            signals[i] = 0
            filteredY[i] = y[i]
            avgFilter[i] = new_mean
            stdFilter[i] = new_std

    print "Thresholding Algo Complete at", str(datetime.now())
    return dict(signals = signals,
                avgFilter = avgFilter,
                stdFilter = stdFilter)
    
def highPassFilter(timeseries, filtersize=2,  inverse=False):
    
        b, a =sgl.butter(1,2*filtersize/fe,btype="highpass")
        if inverse==False:
            z=lfilter(b,a,timeseries)
        
        if inverse==True:
            z=lfilter(a,b,timeseries)
        
        return z
                

def lowPassFilter(timeseries, filter_size, filter_type="ideal", remove_DC=True): #kills DC too
    
    fourier = np.fft.fft(timeseries)
    print "type of fourier ", type(fourier)
    N=len(fourier)
#    freqs = np.arange(N)*fe/N
    fourier1side=fourier[:int(1+N/2)]
    
    print "type fourier1side ", type(fourier1side)

    print "fe", fe
    if filter_type=="ideal":
        
        if remove_DC:
            fourier1side[0]=0
        k = int(filter_size*N/fe)+1
        fourier1side[k:] = 0
#        for fft_value in fourier1side[k:]:
#            fft_value=0
#            
#            
#        while k < len(fourier1side):
##        for k in range(0, len(fourier1side)):
#            if k*fe/N>filter_size:    
#                fourier1side[k]=0
#            if k%int(1000*fe)==0:
#                print "on second ", k/fe
#            k += 1
            
    output=fourier1to2sided(fourier1side, N)
    
    print "type of output ", type(output)
    output_t=np.real(np.fft.ifft(output))
    
    print "type of output(t) ", type(output_t)
    
    return output_t
    
def removeDC(timeseries):
    fourier = np.fft.fft(timeseries)
    fourier[0]=0
    return np.real(np.fft.ifft(fourier))

def loadVoltageFile(fname, interval_length, start_time=0, TEST=True):
    '''
    start_time is in seconds
    interval length is in seconds, formerly 'dt'
    returns t,V
    '''
    pathtofile = 'DataFiles\\'+fname +'.BIN0'
    here = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(here, pathtofile)
    print filename
    voltages = np.memmap(filename, dtype='int16', mode='r+', offset=0, shape=(int(10*fe*start_time)+int(10*fe*interval_length),1))
    times = np.memmap('TempFiles\\times.mymemmap', dtype='float32', mode='r+', shape=(int(10*fe*interval_length),))
    index=0    
    while index < int(10*fe*interval_length):
        times[index] = index/(10*fe)
        if index%10000000==0:
            print "At time", index/(10*fe)
        index += 1
    times = times + start_time
    voltages=voltages[(int(10*fe*start_time)):]
    (T1, l1) = (times, voltages)
    if TEST:
        print "fe interval_length", fe, interval_length, int(10*fe*interval_length)
        t=T1[:int(10*fe*interval_length)]
        V=l1[:int(10*fe*interval_length),0]*Gain
        print len(t), len(V)
        print "TESTING!!!"
    #    pl.plot(t, (10**6)*V, label="pre downsampled")
    #    pl.legend()
    #    pl.show()
        
        f_cut = fe
        print "filtering before downsample"
        V=lowPassFilter(V, f_cut )
        print "done filtering"
        
        print "fe f_nyq f_cut", fe, f_Nyq, f_cut
        t=t[::10]
        
        V=V[::10]
        print "lengths", len(t), len(V)
        pl.plot(t, V, label='Down Sampled & Low Pass Filtered')
        pl.legend()        
        pl.show()
    
        with open('downsampled_stream_30k_s.p', 'wb') as fp:
            pickle.dump( (t,V) , fp )

    if TEST==False: #then set fe back to 10000 at top 
        t=T1[:int(fe*interval_length)]
        V=l1[:int(fe*interval_length),0]*Gain
        
    return (t, V)

def PulseTrigger(stream, start_timestep = 0, end_timestep = -1, lag = 50000, threshold = 5, influence = 0, writeOutput = None ):
    pulseResults = []
    start_timestep = int(start_timestep)
    end_timestep = int(end_timestep)
    end_time = end_timestep/float(fe)
    if end_timestep >= 0:
        pulseResults = thresholding_algo(stream[start_timestep:end_timestep], lag=lag, threshold=threshold, influence=influence)
    else:
        pulseResults = thresholding_algo(stream[start_timestep:], lag=lag, threshold=threshold, influence=influence)
        end_time = len(stream)/float(fe)
        
        
    fig, ax1 = pl.subplots()
    ax1.plot(np.linspace(start_timestep/fe , end_time, len(stream[ start_timestep:end_timestep ] )), stream[ start_timestep:end_timestep ], 'b-', label = "Low Pass Filtered")
    ax1.set_xlabel('Time [s]')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel("Amplitude from Baseline [uV]")
    
    ax2 = ax1.twinx()
    ax2.step(np.linspace(start_timestep/fe , end_time, len(pulseResults["signals"])), pulseResults["signals"], color="red", lw=2, label = "Trigger")
    ax2.set_ylabel('Trigger Value')
    ax2.set_ylim(-1.5,1.5)

    
    fig.tight_layout()
    
    pl.legend(numpoints = 1, loc = 2, frameon=False)
    pl.savefig('UNrefined_trigger.png')
    pl.show()
    
    pulsePeaks_locations = []
    pulses = []
    pulse_locations=[] #list of lists, each element contains a start and end for pulse
    
    for j in range(1,len((pulseResults["signals"]))-2):
        if (pulseResults["signals"][j] == 1 and pulseResults["signals"][j-1] == 0):
                pulse_start_step = j - 1
                time_start = int(pulse_start_step + start_timestep)
        elif (pulseResults["signals"][j] == 0 and pulseResults["signals"][j-1] == 1):
                pulse_end_step = j
                time_end = int(pulse_end_step + start_timestep)
                pulsePeaks_locations.append( (list(stream[time_start:time_end])).index(max(stream[time_start:time_end])))
                pulse_locations.append([pulse_start_step, pulse_end_step])
                pulses.append(stream[time_start:time_end])
    
    #print "\nNumber of Pulses Found: ", len(pulses)    
    #print "Pulse Locations: ", pulse_locations
#    print "peak locations ",pulsePeaks_locations
    #print "Pulses: ", pulses
    z=pulseResults["signals"]
    signal=z
    return (pulses, pulse_locations, pulsePeaks_locations, signal)


def getTrigger(stream, start_timestep, end_timestep, lag, threshold, influence,names='uVf', read=True):
    triggerDataFile = "Trigger_"+names+"_times" + str(start_timestep/fe) + "-" + str(end_timestep/fe) \
                + "_lag" + str(lag) + "_thresh" + str(threshold)  + "_infl" + str(influence) + ".p"
                
    
    if read==True:
        
        try:
            with open(triggerDataFile, 'rb') as fp:
                triggerOutput = pickle.load(fp)
            
            print "\nRead Trigger Data\n" 
            fig, ax1 = pl.subplots()
            ax1.plot(np.linspace(start_timestep/fe , end_timestep/fe, len(stream[ start_timestep:end_timestep ] )), stream[ start_timestep:end_timestep ], 'b-', label = "Low Pass Filtered")
            ax1.set_xlabel('Time [s]')
            ax1.set_ylabel("Amplitude from Baseline [DU]")
            ax2 = ax1.twinx()
            ax2.step(np.linspace(start_timestep/fe , end_timestep/fe, len(triggerOutput[-1])), triggerOutput[-1], color="red", lw=2, label = "Trigger")
            ax2.set_ylabel('Trigger Value')
            ax2.set_ylim(-1.5,1.5)
            fig.tight_layout()
            pl.legend()
            pl.show()        
        
        except:
            print "\nNo Trigger Data found, running algorithm now\n" 
            triggerOutput = PulseTrigger(stream, start_timestep, end_timestep, lag, threshold, influence)
            with open(triggerDataFile, 'wb') as fp:
                pickle.dump(triggerOutput, fp)
                
    else:
        
        triggerOutput = PulseTrigger(stream, start_timestep, end_timestep, lag, threshold, influence)
    
        with open(triggerDataFile, 'wb') as fp:
                pickle.dump(triggerOutput, fp)
                
    return triggerOutput


def MaxTrigger(signal, stream, pulse_locations):
    #takes a 1/0 signal and stream, refines signal to merge events too close to one another, and then set to spike on a maximum of the stream
    deltat=int(1*fe)
    min_separation=int(0.5*fe)
    
    mod_pulse_locations=[]
    for k in pulse_locations:
        if k[0]+2*deltat<len(stream):
            
            short_stream=stream[k[0]:k[1]]
            short_stream_max=np.max(short_stream)
            minsep_int=stream[k[0]:k[1]+min_separation]
            if short_stream_max==np.max(minsep_int):
                mod_pulse_locations.append((k[0],k[0]+deltat))
            

    ref_signal=np.zeros(len(signal))
    for m in mod_pulse_locations:
        for k in range(m[0],m[1]):
            ref_signal[k]=1
                    
    return (ref_signal, mod_pulse_locations)
    
def PlotPulses(stream, pulses, pulse_locations, start_timestep, ext_window = 3):
    for p in range(len(pulses)):
        window = ext_window
        begin = int(start_timestep + pulse_locations[p][0] - window)
        end = int(start_timestep + pulse_locations[p][1] + window)
        extended_timesteps =   np.asarray(range(begin, end))
        extended_times = extended_timesteps/float(fe)
        pl.plot( extended_times , stream[begin:end] , label = "Observed Pulse", color = "red")
        
        pl.xlabel("Time [s]")
        pl.ylabel("Amplitude from Baseline [DU]")
        pl.legend(numpoints = 1, loc = 2, frameon=False)
        pl.show()   

def PlotFitPulses(stream, pulses, pulse_locations, start_timestep, ext_window , init_params , errors, FitFunc, Method = "Powell"  ):
    for p in range(len(pulses)):
        window = ext_window
        begin = int(start_timestep + pulse_locations[p][0] - window)
        end = int(start_timestep + pulse_locations[p][1] + window)
        extended_timesteps =   np.asarray(range(begin, end))
        extended_times = extended_timesteps/float(fe)
        pl.plot( extended_times , stream[begin:end] , label = "Observed Pulse", color = "red")
        
        pl.xlabel("Time [s]")
        pl.ylabel("Amplitude from Baseline [DU]")
        pl.hold(True)
        
        xguess = (np.linspace(begin, end) + start_timestep)/float(fe)
        errors = (end-begin)*[1]
        results = chisqmin( init_params, np.linspace(0, (end-begin)/float(fe) , len(stream[begin:end]) ) , stream[begin:end], errors, FitFunc, Method)
        
        (best_params, uncs, DoF, finalChi, unc_intervals) = results        
        pl.plot(xguess, doubleExp(xguess - extended_times[0], best_params), label = "Best Fit")
        pl.legend(numpoints = 1, loc = 2, frameon=False)        
        pl.show() 
        print
        print "Degrees of Freedom: ", DoF
        print "Reduced Chi Square of Best Fit: ", (finalChi/float(DoF))
        print "Best Parameters: ", best_params
        print "Uncertainties: ", uncs
        print "Uncertainty Intervals: ", unc_intervals
        print "\n\n"


    
def removePulses(stream, pulse_locations, start_timestep = 0, end_timestep = -1, extra_buffer = 50000):
    pureNoiseRegions = []
    if end_timestep < 0:
        end_timestep = len(stream)
    for i in range(len(pulse_locations)):
        if i==0:
            pureNoiseRegions.append(stream[start_timestep:  pulse_locations[i][0] - extra_buffer  ])
        elif i<= len(pulse_locations) - 1:
            pureNoiseRegions.append(stream[pulse_locations[i-1][0] + extra_buffer:  pulse_locations[i][0] - extra_buffer  ])
        if i==len(pulse_locations):
            pureNoiseRegions.append(stream[ pulse_locations[i][1] + extra_buffer : end_timestep])
    for n in pureNoiseRegions:
        print "Time Duration in seconds of nth Noise region: ", len(n)/float(fe)

Tf_s=1.1   #slow fall
Tf_f=0.5  #fast fall
Tr=0.008 #rise time

def template_fast(t, (t_0,)):
    
    exp=theta_exp(t, (t_0,Tr_2amp_fit,Tf2_2amp_fit))
    
    return exp

def template_slow(t, (t_0,)):

    exp=theta_exp(t, (t_0,Tr_2amp_fit,Tf1_2amp_fit))
    
    return exp    

def template_1amp(t, (t_0,)):
    exp=theta_exp(t, (t_0, Tr_fit, Tf_fit))
    return exp

def zeroNoisy(t, (t_0,)):
    n=len(t)
    return np.random.normal(0, 1e-5,n)
#def Double_theta_exp(t, (t_0, Tr1, Tf1, A1, Tr2, Tf2,  A2, b)):
#    shifted_t = np.asarray(t) - t_0    
#    m1=expo_temp1(shifted_t, (Tr1,Tf2))
#    m2=expo_temp1(shifted_t, (Tr2,Tf2)) 
#    return A1*step(shifted_t)*m1 +  A2*step(shifted_t)*m2 + b

def Double_theta_exp(t, (t_0, A1, A2)):
    b=0
    shifted_t = np.asarray(t) - t_0    
    m1=expo_temp1(shifted_t, (Tr1,Tf1))   #slow
    m2=expo_temp1(shifted_t, (Tr2,Tf2))    #fast
    return A1*step(shifted_t)*m1 +  A2*step(shifted_t)*m2 + b


def optimalFiltering2amp(stream, pulse_locations, J, freq,t0=0.005,temp_1=template_slow,temp_2=template_fast, start_timestep = 0, end_timestep = -1, scan_range = 0.04, bounds = None, window = "boxcar", method = "Powell", plots = True ):
    '''    
    J is the noise PSD
    t0 should be the first parameter in init_guess_params   
    init_guess_params should be a list
    stream should be an unfiltered stream
    init_guess_params[0] should be the delay the start of each pulse trigger you expect the rist to start
    '''
    
    def chisquared_2(V, temp1, temp2, a,b, J):
        
        return sum((np.abs(V[1:] - a*temp1[1:]-b*temp2[1:] ))**2 / J[1:])    
    

 
    min_chis={}
    min_t0s={}
    min_amps={}
    pulses={}
    templates={}
    min_A={}
    pulse_time={}
    chi_zeros={}
    
    
    if end_timestep < 0:
        end_timestep = len(stream)
 

    pulse_no=0
    #results_list = []
    for pulse in pulse_locations:   
        strt=int(pulse[0]-scan_range + start_timestep)
        end=int(pulse[1] + start_timestep)
        
        
        trial_pulse=np.asarray(stream[strt : end])
        b, a =sgl.butter(1,2.0/f_Nyq,btype="highpass")
        pulse_size = len(trial_pulse)
        #rem=np.remainder(pulse_size,2*fe)   #uncomment this to allow variable time lengths for pulses
        rem=pulse_size-1*fe   #makes all pulse sizes =1s
        time=np.arange(strt,end-rem)/float(fe)
        pulse_time[pulse_no]=time
        trial_pulse=np.array(stream[strt:int(end-rem)])
#        pulses[pulse_no]=trial_pulse
#        trial_pulse_hp=lfilter(b,a,trial_pulse) #high pass filtered pulse

        pulse_size=len(trial_pulse)
        
        trial_pulse=removeDC(trial_pulse)
        
        pulses[pulse_no]=trial_pulse

        trial_pulse_hp=highPassFilter(trial_pulse)
        window_vals = sgl.get_window(window, pulse_size)
        w_pulse= window_vals*trial_pulse_hp
#        obs_fourier=np.fft.fft(w_pulse)
##        ENBW=1.5*(fe/pulse_size)     #Hanning Window ENBW
        S1=np.sum(window_vals)
        S2 = np.sum( (window_vals**2)  )
        ENBW = fe*S2/(float(S1**2))   

        f, scaled_obs_FT_oneside=getWindowedFFT(w_pulse)
        
        V=scaled_obs_FT_oneside
        
        
        #init_guess_params[0] = (pulse[0])/float(fe) + t0
        #init_guess = tuple(init_guess_params)
        
        
        t0s= np.arange(t0, t0+scan_range, 0.002)
        for t in t0s:
                
            
            time_delay=t
            shifted_time=np.arange(len(time))/float(fe)-t
            temp_fast=temp_2(shifted_time, (t,))                 #(shifted_time, (Tr2,Tf2,0))
            temp_slow=temp_1(shifted_time, (t,))     

            temp_fast=removeDC(temp_fast)
            temp_slow=removeDC(temp_slow)
                                                    #theta_exp(shifted_time, (Tr1, Tf2,0))
                        
            temp_fast_hp=lfilter(b,a,temp_fast)
            temp_slow_hp=lfilter(b,a,temp_slow)
            
            freqs, S_f=getWindowedFFT(temp_fast_hp)
            freqs, S_s=getWindowedFFT(temp_slow_hp)
            
#            temp_fast_fft=np.fft.fft(temp_fast_hp)
#            #temp_fast_fft=temp_fast[::2]    #downsampled
#            temp_fast_1sideFT=temp_fast_fft[:int(1+len(temp_fast_fft)/2)]
#
#            temp_slow_fft=np.fft.fft(temp_slow_hp)
#            #temp_slow_fft=temp_slow[::2]
#            temp_slow_1sideFT=temp_slow_fft[:int(1+len(temp_slow_fft)/2)]
#            
#            S_f=temp_fast_1sideFT/(S1*np.sqrt(ENBW))
#            S_s=temp_slow_1sideFT/(S1*np.sqrt(ENBW))
            
            S_f[0]=0
            S_s[0]=0
            
            c_fv=np.real(np.sum(np.conjugate(S_f)*V/J))
            c_fs=np.real(np.sum(np.conjugate(S_f)*S_s/J))
            c_ff=np.sum((np.abs(S_f)**2)/J)
            c_sf=np.conjugate(c_fs)
            c_ss=np.sum((np.abs(S_s)**2)/J)
            c_sv=np.real(np.sum(np.conjugate(S_s)*V/J))
            
#            matrix = np.matrix([[c_ff, c_fs],[c_sf,c_ss]])
#            matrix_inv = np.linalg.inv(matrix)
#            obs_vector = [c_fv, c_sv]
#            best_theta = np.dot(matrix_inv,obs_vector) 
#            best_theta=np.asarray(best_theta)
            
#            pl.plot(freq, np.abs(S_f), label="fast")
#            pl.hold(True)
#            pl.plot(freq, np.abs(S_s), label="slow")
#            pl.xlabel("frequency")
#            pl.grid()
#            pl.legend()
#            pl.show()
#            
#    
#            pl.plot(time, temp_fast, label="fast")
#            pl.hold(True)
#            pl.plot(time,temp_slow,label="slow")
#            pl.grid()
#            pl.legend()
#            pl.show()        
#            
#            pl.plot(time, temp_fast, label="fast")
#            pl.hold(True)
#            pl.plot(time,temp_slow,label="slow")
#            pl.grid()
#            pl.legend()
#            pl.show()
            
            
            #print "coefficients are ", c_fv, "and ", c_ss, "and ", c_sf, "and ", c_ff
            
            C_fv=c_fv/c_ff      #normalized matrix values
            C_fs=c_fs/c_ff
            C_sf=c_sf//c_ss
            C_sv=c_sv/c_ss

            #print "coefficients are ", C_fv, "and ", C_sv, "and ", C_sf

                        
            a_f=(C_fv-C_fs*C_sv)/(1-C_fs*C_sf)    #optimal amplitudes
            a_s=C_sv-a_f*C_sf
            
            template=a_s*temp_slow+a_f*temp_fast

            
            chi=chisquared_2(V, S_f, S_s, a_f, a_s, J)
            chi_zero=chisquared_2(V, S_f, S_s, 0, 0,J)
#            print "best theta ", best_theta, best_theta[0]
#            a_f=best_theta[0][0]
#            a_s=best_theta[0][1]
            
            if t==t0s[0]:
                min_chi=chi
                
            if chi<=min_chi:
                min_chi=chi
                min_t=t
                
                min_chis[pulse_no]=min_chi
                min_amps[pulse_no]=(a_s,a_f)
                chi_zeros[pulse_no]=chi_zero
                min_t0s[pulse_no]=t
                templates[pulse_no]=template
                min_A[pulse_no]=np.max(template)-np.min(template)
                
            
            
#            print "fast & slow amplitudes are : ", a_f, " and ", a_s
##            print "ANDD best theta is ", best_theta
#
#            
#            #template=sglTemplateFunc(shifted_time, (0, a_f, a_s))
#            
#            
##            pl.plot(time,trial_pulse, label="Observation" )
##            pl.hold(True)
##            pl.plot(time, template, label="best fit template")
##            pl.legend()
##            pl.show()
##            
##            pl.loglog(freq, np.abs(V)**2, label="Observed FT")
##            pl.hold(True)
##            pl.loglog(freq, np.abs(a_s*S_s+a_f*S_f)**2, label="Best Match")
##            pl.hold(True)
##            pl.loglog(freq, J, label="Error")
##            pl.grid()
##            pl.legend()
##            pl.show()
#            
#            print "this fit had chisquared of ", chi
            
#            best_FT=np.fft.fft(template)
#            #best_FT=best_FT[::2] #downsample
#            best_FT_1side=best_FT[:int(1+len(best_FT)/2)]
#            best_FT_1side=best_FT/(S1*np.sqrt(ENBW))
            
            best_FT=getWindowedFFT(template)
        
        pulse_no+=1
        
    pulsedictlist=[]
    counter = 0
    for j in pulses:
        if counter%30==0:
            pl.plot(pulse_time[j], pulses[j], label="Observed Pulse", color='blue')
            pl.hold(True)   
            pl.xlabel("Time [s]")
            pl.ylabel("Amplitude [DU]")
            pl.plot(pulse_time[j], templates[j],label="Best Fit", color='red')
#            pl.plot(time, sglTemplateFunc(time, results["Best Fit Params"]), label="Best Fitting Template", color='r')
#            pl.hold(True)        
            pl.legend()
            pl.show()
            
            print "chi sq for the fit was ", min_chis[j]
            print "compared to chi zero of ", chi_zeros[j]
            print "with amps of ", min_amps[j]
            print "with optimal time delay of ", min_t0s[j]
        
        pulse_dict={}
        pulse_dict["chi"]=min_chis[j]
        pulse_dict["amplitude"]=min_A[j]
        pulse_dict["obs pulse"]=pulses[j]
        pulse_dict["template"]=templates[j]
        pulse_dict["slow amplitude"]=min_amps[j][0]
        pulse_dict["fast amplitude"]=min_amps[j][1]
        pulse_dict["time"]=pulse_time[j]
        pulse_dict["delay"]=min_t0s[j]
        pulse_dict["window"]=window
    
        pulsedictlist.append(pulse_dict)
        counter+=1

#        print "Bounds: ", str(bounds)
#        if bounds!=None:
##            t_array=np.array(range(int(strt),int(strt+2*scan_range),500)) #timestep units
##            print t_array
##            result_list = []
##            chi_list = []
##            for t in t_array:
##                init_guess_params[0] = t/float(fe)
##                res = chisqmin( tuple(init_guess_params) , time , scaled_obs_FT_oneside, np.sqrt(J), sglTemplateFunc, meth = method , FreqSpace=True, bounds=bounds)        
##                result_list.append(res)
##                chi_list.append(res[ "Chi Sq of Best Fit" ] )
##            
##            min_chi = min(chi_list)
##            index = chi_list.index(min_chi)
##            results = result_list[index]    
#            bounds[0] = tuple(np.asarray(bounds[0]) + init_guess[0])
##            bounds = tuple(bounds)
#            results = chisqmin( init_guess , time , scaled_obs_FT_oneside, np.sqrt(J), sglTemplateFunc, meth = method , FreqSpace=True, bounds=tuple(bounds)  )      
#        else:
#            results = chisqmin( init_guess , time , scaled_obs_FT_oneside, np.sqrt(J), sglTemplateFunc, meth = method , FreqSpace=True, bounds=bounds)        
#        print results
        
        
        
#        if plots:
#            pl.loglog(freq, np.abs(fourierTempForPlotting(time, results["Best Fit Params"], sglTemplateFunc )) ,label="Template Best Fit FT")
#            pl.hold(True)
#            pl.loglog(freq, np.abs(scaled_obs_FT_oneside), label="observed FT") 
#            pl.hold(True)
#            pl.loglog(freq,np.sqrt(J),label="Errors")
#            pl.hold(True)
#            pl.loglog(freq,np.abs(fourierTempForPlotting(time, init_guess, sglTemplateFunc )),label="Initial Guess")
#            pl.xlabel("Frequency [Hz]")
#            pl.ylabel("LPSD [V/sqrt(Hz)]")        
#            pl.legend()
#            pl.show()
#            
#            pl.plot(time, sglTemplateFunc(time, init_guess), label="Initial Guess", color='green')
#            pl.hold(True)   
#            pl.xlabel("Time [s]")
#            pl.ylabel("Amplitude [uV]")
#            pl.plot(time, trial_pulse,label="Observed Pulse", color='blue')
#            pl.plot(time, sglTemplateFunc(time, results["Best Fit Params"]), label="Best Fitting Template", color='r')
#            pl.hold(True)        
#            pl.legend()
#            pl.show()
#        
#        results_list.append(results)
        
    return pulsedictlist


def fourier1to2sided(fourieroneside,N):
#    onesidelist=list(fourieroneside)
    
    twosideMemMap = np.memmap('TempFiles\\fourier1to2sided.mymemmap', dtype=np.complex64, mode='r+', shape=(N,))
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
   

def butterworthLPFilter(stream, wcut):
    fourier=np.fft.fft(stream)
    freq=fe*np.arange(len(stream))/len(stream)
    
    freqoneside=freq[0:int(1+len(stream)/2)]
    fourieroneside=fourier[0:int(1+len(stream)/2)]
        
    butterworth=1/(np.sqrt(1+10*(freqoneside/wcut)**2))
    filteredfourier=butterworth*fourieroneside
    filteredtwoside=fourier1to2sided(filteredfourier,len(stream))
    filteredtwoside[0]=0
    
    filteredstream=np.fft.ifft(filteredtwoside)   
    
    return np.real(filteredstream)
    
    
def OFTransferFn(stream, J,freq , parameters,window="boxcar", sglTemplateFunc=A_theta_exp):
    '''
    stream should be unfiltered
    '''
    
    
    time=np.arange(0,len(stream))/fe
    #two_sec=np.arange(0,2*fe)/fe
    one_sec=np.arange(0,fe)/fe
    
    #assert len(two_sec)%2==0
    
#    freq=np.arange(0,len(J))*fe/(len(J))
    #freqs=np.linspace(0,5000,)
    t_delay=0.5
#    s=sglTemplateFunc(one_sec, (t_delay, 0.005,0.65,1.2,0))
    
    if sglTemplateFunc==AB_theta_exp:
        Tr=parameters[1]
        Tf1=parameters[2]
        Tf2=parameters[3]
        A=0.5
        B=A
        b=0
        s=sglTemplateFunc(one_sec, (t_delay, 0.008,0.5,2.5,0.6,0.4,0))
        tM=2.05*np.pi*t_delay
        

        
    if sglTemplateFunc==A_theta_exp:
        Tr=parameters[1]
        Tf=parameters[2]
        A=1
        b=0
        s=sglTemplateFunc(one_sec, (t_delay, 0.008,1.1,A,0.0))        
        tM=2.05*np.pi*t_delay
    
    window=sgl.get_window(window, len(s))
    
    S1=np.sum(window)
    S2=np.sum(window**2)
#    b, a =sgl.butter(1,0.0004,btype="highpass")
    
    b, a =sgl.butter(1,0.0004,btype="highpass")
    ENBW=fe*S2/(S1**2)
    print "ENBW is in HZ ", ENBW
    
#    b, a =sgl.butter(1,1,btype="highpass")
    
    t_s=np.arange(len(s))/fe

    
##    w, STransfer=sgl.freqz(b,a,worN=(1+len(two_sec)/2))
#    
#    sFT=np.fft.fft(s)
#    sFToneside=sFT[0:int(len(sFT)/2+1)]
#    
#    
#    filteredsFToneside =STransfer*sFToneside
#    filteredsFToneside[0]=0
#    
#    sftwoHztwoside=fourier1to2sided(filteredsFToneside, len(sFT)) #undownsampled, 2Hz filtered, 2 second long interval of template
#
#    
#    filteredsFToneside=filteredsFToneside[0:len(filteredsFToneside):2]
#    sftwoHz=filteredsFToneside
    
#    s_highFiltered_timeDomain = lfilter(b,a,s)*window
    s_highFiltered_timeDomain=highPassFilter(s)
#    stream=lfilter(b,a,stream)
    stream=highPassFilter(stream)    
#    pl.plot(t_s, s, 'r', t_s, s_highFiltered_timeDomain,'b')
#    pl.grid()
#    pl.title("normal and high pass filtered")
#    pl.show()
    
    s_FFT_twoHz_twosided = np.fft.fft(s_highFiltered_timeDomain)
    s_FFT_twoHz_twosided[0]=0
    #s_FFT_twoHz_twosided = s_FFT_twoHz_twosided[0:len(s_FFT_twoHz_twosided):2] #factor of 2 for downsampling
#    fe_ds=0.5*fe
#    ENBW_ds=0.5*ENBW
#    S1_ds=0.5*
    
#    testfreq, s_FFT_twoHz_onesided = getWindowedFFT(s_highFiltered_timeDomain, "density")
    s_FFT_twoHz_onesided = s_FFT_twoHz_twosided[0:1+len(s_FFT_twoHz_twosided)/2]/(S1*np.sqrt(ENBW))  #scaled by the enbw
#    sftwoHztwoside=fourier1to2sided(, len(sFT)) #undownsampled, 2Hz filtered, 2 second long interval of template
#    s_FFT_twoHz_onesided = s_FFT_twoHz_twosided[0:1+len(s_FFT_twoHz_twosided)/2]

#    sftwoHz_t=np.real(np.fft.ifft(sftwoHztwoside))
    
    oneovrh=2*np.sum((np.abs(s_FFT_twoHz_onesided)**2)/J)

    h=1.0/oneovrh
    print "Resolution after OF (h): " , h
    print "Resolution before OF: " , 2*np.sum(J)
    
    
    Nf=(h**2)*(np.abs(s_FFT_twoHz_onesided)**2)/J
    pl.loglog(freq, J, label='Original Noise PSD')
    pl.hold(True)
    pl.loglog(freq, Nf, label='Optimal Filtered Noise PSD')
    pl.legend()
    pl.show()
    
    
    phases = []
    for i in range(len(freq)):
        phases.append(np.complex(0,-freq[i]*tM))
        
    phases=np.array(phases)
#    h=1
    H=h*np.conjugate(s_FFT_twoHz_onesided)*np.exp(phases)/J
    print "the sum of the output signal is ", np.sum(2*np.abs(H*s_FFT_twoHz_onesided))
    H_2sided=fourier1to2sided(H,int(len(s)))
    h_t_ds=np.fft.ifft(H_2sided)
    print "peak to peak gain is ", np.sum(np.abs(h_t_ds))
    pl.plot(np.arange(fe)/len(h_t_ds),h_t_ds,label="Impulse Response Function")
    pl.grid()
    pl.legend()
    pl.show()
    pl.loglog(freq, np.abs(s_FFT_twoHz_onesided)**2, label="Template PSD unscaled" )
    pl.legend()
    pl.show()
    pl.loglog(freq, np.abs(H), label="Transfer Fn Magnitude" )
    pl.legend()
    pl.show()
    Htwoside=fourier1to2sided(H, N=int(fe))
    
    transferred_template_FFT=[]
    

#    for j in range(0,len(Htwoside)):
#        i=2*j
#        transferred_template_FFT.append(Htwoside[j]*s_FFT_twoHz_twosided[i])
#        transferred_template_FFT.append(Htwoside[j]*s_FFT_twoHz_twosided[i+1])
#        
#    transferred_template_FFT=np.array(transferred_template_FFT)
    
    transferred_template_FFT=Htwoside*s_FFT_twoHz_twosided
    
    
    
    trans_s=np.fft.ifft(transferred_template_FFT)
    trans_s=np.real(trans_s)
    
    t_s=np.arange(len(trans_s))/fe
    pl.plot(np.linspace(0,1, int(1*fe)),trans_s, label="transferred template")
    pl.hold(True)
    pl.plot(np.linspace(0,1, int(fe)), highPassFilter(s ), label="high pass filtered template" )
    pl.hold(True)
    pl.plot(np.linspace(0,1, int(fe)), s , label="original template" )
    pl.legend()
    pl.show()
    

    
    h_t=np.fft.ifft(Htwoside)
    
    h_tfirsthalf=h_t[:int(len(h_t)/2)]
    h_tsecondhalf=h_t[int(len(h_t)/2):]
    w_h_tfirsthalf=np.cos(np.linspace(0, np.pi/2.0, len(h_tfirsthalf)))*h_tfirsthalf
    w_h_tsecondhalf=np.cos(np.linspace(3*np.pi/2.0, 2*np.pi, len(h_tsecondhalf) ) )*h_tsecondhalf
    zeros=np.zeros(len(h_t))
    

    
    h_padded=np.concatenate( (w_h_tfirsthalf, zeros, w_h_tsecondhalf),0)

    
    H_padded=np.fft.fft(h_padded)
    
    
    
    H_pad_onesided=H_padded[:int(1+len(H_padded)/2)] #transfer function, length 10k, M=10,000=1s
    
    even_seconds = len(stream)%int(2*fe)==0
    odd_seconds = len(stream[int(1*fe):]) %int(2*fe)==0
    assert even_seconds 
    
    stream_segments1=np.split(stream, int(len(stream)/(2*fe)))
    #stream_segments2=np.split(stream[(1*fe):], int(len(stream)/(2*fe)))
    stream_segments2=np.split(stream[int(1*fe):], np.arange(int(2*fe), len(stream[int(1*fe):]), int(2*fe) ))
    stream_segments2=stream_segments2[:-1]

    
    
    streamcheck1=[]
    
    for i in range(len(stream_segments1)):
        length=len(stream_segments1[i])
        streamcheck1.append(stream_segments1[i][:int(length/2)])
        if i<len(stream_segments2):  
            streamcheck1.append(stream_segments2[i][:int(length/2)])
            
    streamcheck1.append(stream_segments1[-1][length/2:])
    streamcheck1=np.array(streamcheck1)

    np.append
    OF_stream_segments1=[]
    for j in stream_segments1:
        segment_fft=np.fft.fft(j)
        segment_fft_oneside=segment_fft[:len(segment_fft)/2+1]
        OF_fft=H_pad_onesided*segment_fft_oneside
        
        OF_fft_twoside=fourier1to2sided(OF_fft, len(j))
        
        OF_stream1=np.fft.ifft(OF_fft_twoside)
        OF_stream1=np.real(OF_stream1)
        
        OF_stream_segments1.append(OF_stream1)   #list of filtered 2s intervals
        
    
    OF_stream_segments2=[]
    for j in stream_segments2:
        segment_fft=np.fft.fft(j)
        segment_fft_oneside=segment_fft[:len(segment_fft)/2+1]
        OF_fft=H_pad_onesided*segment_fft_oneside
        
        OF_fft_twoside=fourier1to2sided(OF_fft, len(j))
        
        OF_stream2=np.fft.ifft(OF_fft_twoside)
        OF_stream2=np.real(OF_stream2)
        
        OF_stream_segments2.append(OF_stream2)
        
    
    
    OF_full_stream=[]
    segment_length=len(OF_stream_segments1[0])
    zeros=(segment_length/4)*[0]

#    OF_full_stream.append(OF_stream_segments1[0][:segment_length/4])
    OF_full_stream.append(zeros)
    
    assert segment_length%2==0
    count1 = 0
    count2 = 0
    for j in range(len(OF_stream_segments1)):
        
        OF_full_stream.append(OF_stream_segments1[j][segment_length/4:3*segment_length/4])
        if len(OF_stream_segments1[j])!=segment_length:
            print "error!!"
            
        count1+=1
        if j<len(OF_stream_segments2):
            
            OF_full_stream.append(OF_stream_segments2[j][segment_length/4:3*segment_length/4])
            count2+=1
            if len(OF_stream_segments2[j])!=segment_length:
                print "error!!"            
    

    #OF_full_stream.append(OF_stream_segments1[-1][3*segment_length/4:])
    OF_full_stream.append(zeros)
    flat_OF_fullstream =[]
    for sublist in OF_full_stream :
        for item in sublist:
            flat_OF_fullstream.append(item)
            
    flat_OF_fullstream = np.asarray(flat_OF_fullstream)
        
        
    #OF_full_stream=np.array(OF_full_stream)
    
#    OF_full_stream=OF_full_stream.flatten()
#    
#    OF_full_stream=np.array(OF_full_stream)
        
    pl.plot(time, flat_OF_fullstream, label="Optimally Filtered Stream")
    pl.hold(True)
    pl.plot(time, stream, label="Unfiltered Stream")
    pl.legend()
    pl.show()
    pl.plot(time, flat_OF_fullstream, label="Optimally Filtered Stream")
    pl.legend()
    pl.show()
    pl.plot(time, flat_OF_fullstream, label="Optimally Filtered Stream")
    pl.legend()
    pl.xlim([10,14])
    pl.show()
    pl.plot(time, stream, label="Unfiltered Stream")
    pl.legend()
    pl.show()
    
    return flat_OF_fullstream, H_pad_onesided

    
def optimalFiltering1amp(stream, pulse_locations, J, freq,t0=-0.10,temp_1=template_1amp, start_timestep = 0, end_timestep = -1, scan_range = 0.30, bounds = None, window = "boxcar", method = "Powell", plots = True ):
    #t0 was 0.01, scan range 0.35
    
    def chisquared(V, expected, a, J):
        return sum((np.abs(V[1:] - a*expected[1:] ))**2 / J[1:])   
    
    min_chis={}
    min_t0s={}
    min_amps={}
    pulses={}
    templates={}
    min_A={}
    pulse_time={}
    chi_zeros={}
    amp_error={}
    
    
    if end_timestep < 0:
        end_timestep = len(stream)
 

    pulse_no=0
    #results_list = []
    for pulse in pulse_locations:   
        strt=int(pulse[0] + start_timestep)
        end=int(pulse[1] + start_timestep)
        
        trial_pulse=np.asarray(stream[strt : end])
        pulse_size = len(trial_pulse)
        
        #rem=np.remainder(pulse_size,2*fe)   #uncomment this to allow variable time lengths for pulses
        rem=pulse_size-1*fe   #makes all pulse sizes =1s
        time=np.arange(strt,end-rem)/float(fe)
        
        pulse_time[pulse_no]=time
        trial_pulse=np.array(stream[strt:int(end-rem)])
        
#        pulses[pulse_no]=trial_pulse
#        trial_pulse_hp=lfilter(b,a,trial_pulse) #high pass filtered pulse

        pulse_size=len(trial_pulse)
        
#        print "pulse size ", pulse_size
        
        trial_pulse=removeDC(trial_pulse)
        
        pulses[pulse_no]=trial_pulse

#        trial_pulse_hp=highPassFilter(trial_pulse)
        trial_pulse_hp=highPassFilter(trial_pulse)
        window_vals = sgl.get_window(window, pulse_size)
        w_pulse= window_vals*trial_pulse_hp
#        obs_fourier=np.fft.fft(w_pulse)
##        ENBW=1.5*(fe/pulse_size)     #Hanning Window ENBW
        S1=np.sum(window_vals)
        S2 = np.sum( (window_vals**2)  )
        ENBW = fe*S2/(float(S1**2))
#        print "ENBW IN HZ is ", ENBW
#        
#        scaled_obs_FT=obs_fourier/(S1*np.sqrt(ENBW))
#        #scaled_obs_FT=2*scaled_obs_FT[0::int(len(scaled_obs_FT)/fe)]   #downsampling occurs here
#        scaled_obs_FT_oneside=scaled_obs_FT[0:int(1+fe/2)]
#        V=scaled_obs_FT_oneside
        freqs, V=getWindowedFFT(w_pulse)
        V[0]=0
        
        
        #init_guess_params[0] = (pulse[0])/float(fe) + t0
        #init_guess = tuple(init_guess_params)
        
        
        t0s= np.arange(t0, t0+scan_range, 0.001) 
        for t in t0s:
                
            
            time_delay=t
            shifted_time=np.arange(len(time))/float(fe)-t
            temp=temp_1(shifted_time, (t,))                 #(shifted_time, (Tr2,Tf2,0))
                 

            temp=removeDC(temp)
                                                    #theta_exp(shifted_time, (Tr1, Tf2,0))
                        
#            temp_hp=lfilter(b,a,temp)
            temp_hp=highPassFilter(temp)
            
            
#            temp_fft=np.fft.fft(temp_hp)
#            #temp_fast_fft=temp_fast[::2]    #downsampled
#            temp_1sideFT=temp_fft[:int(1+len(temp_fft)/2)]
#
#
#            
#            S=temp_1sideFT/(S1*np.sqrt(ENBW))
            
            freqs, S=getWindowedFFT(temp_hp)
            
            S[0]=0
            
            
            
            numer = 0
            denom = 0
            
            numer = sum(( np.real( np.conjugate(S) * V )/J ))
            denom = sum(np.abs(S)**2/J)
            ahat_t = float(numer)/float(denom)
                        
            
            template=ahat_t*temp

            
            chi=chisquared(V, S, ahat_t, J)
#            print "best theta ", best_theta, best_theta[0]
#            a_f=best_theta[0][0]
#            a_s=best_theta[0][1]
            chi_zero=chisquared(V, S, 0,J)

            
            if t==t0s[0]:
                min_chi=chi
                
            if chi<=min_chi:
                min_chi=chi
                min_t=t
                
                min_chis[pulse_no]=min_chi
                min_amps[pulse_no]=ahat_t
                min_t0s[pulse_no]=t
                templates[pulse_no]=template
                min_A[pulse_no]=np.max(template)-np.min(template)
                chi_zeros[pulse_no]=chi_zero
                amp_error[pulse_no] = denom #algebraic result for error 
                
            
            
#            print "fast & slow amplitudes are : ", a_f, " and ", a_s
##            print "ANDD best theta is ", best_theta
#
#            
#            #template=sglTemplateFunc(shifted_time, (0, a_f, a_s))
#            
#            
##            pl.plot(time,trial_pulse, label="Observation" )
##            pl.hold(True)
##            pl.plot(time, template, label="best fit template")
##            pl.legend()
##            pl.show()
##            
##            pl.loglog(freq, np.abs(V)**2, label="Observed FT")
##            pl.hold(True)
##            pl.loglog(freq, np.abs(a_s*S_s+a_f*S_f)**2, label="Best Match")
##            pl.hold(True)
##            pl.loglog(freq, J, label="Error")
##            pl.grid()
##            pl.legend()
##            pl.show()
#            
#            print "this fit had chisquared of ", chi
            
            best_FT=np.fft.fft(template)
            #best_FT=best_FT[::2] #downsample
            best_FT_1side=best_FT[:int(1+len(best_FT)/2)]
            best_FT_1side=best_FT/(S1*np.sqrt(ENBW))
            
        
        pulse_no+=1
        
    pulsedictlist=[]
    counter = 0
    for j in pulses:
        if counter%30==0:
            pl.plot(pulse_time[j], pulses[j], label="Observed Pulse", color='blue')
            pl.hold(True)   
            pl.xlabel("Time [s]")
            pl.ylabel("Amplitude [DU]")
            pl.plot(pulse_time[j], templates[j],label="Best Fit", color='red')
#            pl.plot(time, sglTemplateFunc(time, results["Best Fit Params"]), label="Best Fitting Template", color='r')
#            pl.hold(True)        
            pl.legend()
            pl.show()
        
            print "chi sq for the fit was ", min_chis[j]
            print "compared to chi_zero of ", chi_zeros[j]
            print "with amps of ", min_amps[j]
            print "with optimal time delay of ", min_t0s[j]
        
        pulse_dict={}
        pulse_dict["chi"]=min_chis[j]
        pulse_dict["peak"]=min_A[j]
        pulse_dict["obs pulse"]=pulses[j]
        pulse_dict["template"]=templates[j]
        pulse_dict["amplitude"]=min_amps[j]
        pulse_dict["time"]=pulse_time[j]
        pulse_dict["delay"]=min_t0s[j]
        pulse_dict["window"]=window
        pulse_dict["amp error"]=amp_error[j] 
    
        pulsedictlist.append(pulse_dict)
        counter+=1
            
    return pulsedictlist

    
def getBufferedSignal(signal, buff = 7, plots = True ):
    '''
    signals = Trigger output and buffer_size in seconds 
    returns extended trigger values
    '''
    print "length in buffered sig",  len(signal) 
    buff=4  
    buffsignal=[]
    j=0        
    while j<(int(len(signal)/(fe))):
        z=signal[int(j*fe):int((j+1)*fe)]
        maximum=max(z)
        minimum=min(z)
        if maximum > 0.25 or minimum < -0.25:
            for l in range(0, buff + 1):
                buffsignal.append(1)
                j+=1
        else:
            buffsignal.append(0)
            j+=1
            
    if plots:
        pl.title("Buffered Signal Indicating Pulses")
        pl.plot(np.arange(len(signal))/fe, signal)
        pl.xlabel("Seconds")
        pl.ylabel("Trigger with Buffer [On/Off]")
        pl.ylim([-1.5,1.5])
        pl.show()
    return buffsignal    
    
def getPureNoiseIntervals(stream, buffered_signal, plots = False):
    '''
    N is length of stream
    Finds pure noise intervals that do not contain any pulse triggers
    returns list of lists, each element is the stream values for a noise interval
    '''
    N=len(stream)
    print len(buffered_signal)
    print "N", N, fe
    purenoiseints=[]
    k=0
    while k < int(N/(fe)):
        intrvl=stream[int(k*fe):int((k+1)*fe)]
        if buffered_signal[k]==0:
            purenoiseints.append(intrvl)
            k+=1
        else:
            k+=1
    if plots and len(purenoiseints)>0:
        #plot pure noise interval
        x = np.random.choice(range(len(purenoiseints)))
        pl.plot(np.arange(0,1,1/fe),purenoiseints[x])
        pl.title("Pure Noise Interval " + str(x))
        pl.xlabel("Seconds")
        
        pl.show()
    
    return purenoiseints

def timeUncertainty(stream, buffered_signal):
    
    mean_sigma=0
    
    purenoiseints=getPureNoiseIntervals(stream, buffered_signal)
    
    for j in purenoiseints:
        
        int_lp=j
#        int_lp=(butterworthLPFilter(j,100))
        sigma=np.std(int_lp)
        mean_sigma+=sigma
    
    mean_sigma=mean_sigma/(len(purenoiseints))

    return mean_sigma    
       
def noisePSDrandomWelch( purenoiseints , window = "boxcar", samplesize = -1, sigmaCut = 5, cutOrder=  1,plots = True):
    '''
    purenoiseints - takes a set of pure noise streams
    Computes the noise psd averaging a random subset of the noise intervals
    returns the freqs and the resulting averaged noise PSD
    '''
    if samplesize < 0:
        samplesize = int( len(purenoiseints ) /2 )
    if samplesize == "all":
        samplesize = len(purenoiseints)
                    
    sampleindices = np.random.choice( range(len(purenoiseints) ), samplesize , replace=False)
    samplePSDs=[]
    count = 0
    PSDdict = {}
    for j in sampleindices:
        #compute psds for random subset
        interval = purenoiseints[j]
        interval=removeDC(interval)
#        f, PSD = sgl.periodogram(interval,fs=fe, window = window)
#        print "PSD IS ", PSD
        f, int_FT=getWindowedFFT(interval, scaling="density")
        PSD=np.abs(int_FT)**2
        
#        f=np.arange(0,len(PSD))*fe/(len(PSD))        
        

        
        
        PSDdict[j] = PSD
        samplePSDs.append(PSD)
        if count==0:
            freqarray=f
        count+=1
        
    avgPSD=[]
    numfreqs=len(freqarray)
    
    for m in range(numfreqs): #compute avg
        k=0
        for j in range(samplesize):
            k+=samplePSDs[j][m]/float(samplesize)
        avgPSD.append(k)
        
    avgPSD=np.array(avgPSD)
    
        
    DoF = len(avgPSD)
    
    
    ChiSq_sigma = np.sqrt(2*DoF)
    thresholdValue = ChiSq_sigma*sigmaCut + DoF
    newindices=sampleindices[:]           
    newindices=list(newindices)
    PSD_orders_list = [avgPSD]
    list_indicesLists = [newindices]
    size_list = [len(sampleindices)]
    if sigmaCut > 0: #applies chi sq cut and recomputes the avg with outliers removed
        for order in range(cutOrder):     
            print "Doing Chi Sq "+str(order+1) +" order cut now"
            starting_length = len(newindices)
            values_to_remove = []
            for i in range(len(newindices)):
                noiseint = removeDC(purenoiseints[newindices[i]])
                frequencies, FT = getWindowedFFT(noiseint, scaling="density", window=window, oneSided=True)         
                chi=chisq(FT[1:], 0, PSD_orders_list[order][1:])

                if chi > thresholdValue:
                    values_to_remove.append(newindices[i])
#                    newindices.remove(newindices[i])
#                    print "Removing interval: " + str(newindices[i])
#                    print "Which had Chi Square value: " + str(chi) 
#                    print "Threshold value :" + str(thresholdValue )
                    
            for val in values_to_remove:
                newindices.remove(val)
#            print "Starting Length: " + str(starting_length)
#            print "Removed: " + str(starting_length - len(newindices))
            
            NewSamplePSDs=[]
            for k in newindices:
                NewSamplePSDs.append(PSDdict[k])
            avgPSDv2=[]
            numfreqs=len(freqarray)
            for m in range(numfreqs):
                k=0
                for j in range(len(newindices)):
                    k+=NewSamplePSDs[j][m]/float(len(newindices))
                avgPSDv2.append(k)
            avgPSDv2=np.array(avgPSDv2)
            PSD_orders_list.append(avgPSDv2)
            list_indicesLists.append(newindices)
            size_list.append(len(newindices))
                
#    LPSD=np.sqrt(avgPSD)
    if plots and len(purenoiseints) > 0 and samplesize > 0:
        pl.title('PSD averaged by Welch with sample size  ' + str(samplesize))
        pl.xlabel("Frequency[Hz]")
        pl.ylabel("PSD DU^2/Hz]")
        pl.loglog(freqarray,avgPSD, label="Pulses removed")
#        pl.hold(True)
#        f, Pwelch_spec = sgl.welch(uV[start_timestep:end_timestep], fe, nperseg=fe)
#        pl.loglog(f, Pwelch_spec,label="with Pulses")
#        pl.ylim([10**-6,10**2])            
#        pl.legend(frameon = False)
        pl.show()
        
        if sigmaCut > 0:
            pl.xlabel("Frequency[Hz]")
            pl.ylabel("PSD [V^2/Hz]")
            for k in range(len(PSD_orders_list)):
                pl.loglog(freqarray, PSD_orders_list[k], label="Chi Sq Cut Order " + str(k) + "\nSample Size " + str(size_list[k]))
                pl.hold(True)
            pl.ylim([10**-6,10**2])            
            pl.legend()
            pl.savefig('NoisePSD.png')
            pl.show()
            
    
    J=PSD_orders_list[-1]
#        print "Sample Indices" + str(sampleindices)
#        print "New sample indices" + str(newindices)
    return freqarray, J, list_indicesLists[-1]
    
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



def NoiseChiSqHist(purenoiseints, noisePSD, indices_used, window = "boxcar", plots = True, bins=10):
    '''
    Computes and plots the chi square distribution of pure noise intervals
    Returns list of Chi Sq values for each interval used
    '''
    #CALCULATE S2
    PSDchisqvals=[]
    print "Indices used - NoiseChiSqHist " + str(indices_used)
    for m in indices_used:
        z = removeDC(purenoiseints[m])
        FTfreq, FT = getWindowedFFT(z, scaling="density", window=window, oneSided=True)
        chi=chisq(FT,0,noisePSD)
        PSDchisqvals.append(chi)

    if plots:
        pl.hist(PSDchisqvals, bins=bins)
        pl.xlabel("Chi Squared")
        pl.ylabel("No. of Noise Intervals")
        pl.show()

    print "mean chi is ", np.mean(PSDchisqvals)
    print "max chi is ", np.max(PSDchisqvals)
    print "sample size is ", len(indices_used)
    return PSDchisqvals

def getNoisePSD(stream, triggerSignal, win = "boxcar", buffer_length = 4, sample_size = "all", chiSq_sigmaCut = 5, cut_order = 1, plots = True, allowRead = True ):
#    assert len(triggerSignal)==len(stream)
    print "Length of triggerSignal and length of stream ", len(stream), len(triggerSignal)
    noiseFile =  "NoisePSDdata_"  + "buffer" + str(buffer_length) + "_samplesize" + str(sample_size) + "_sigma" + str(chiSq_sigmaCut) +  ".p"
    if allowRead:
        try:
            with open(noiseFile, 'rb') as fp:
                (freqarray, J) = pickle.load(fp)
            
            print "\nRead Noise PSD Data\n" 
            pl.loglog(freqarray,J, label="Noise PSD")
            pl.xlabel("Frequency[Hz]")
            pl.ylabel("PSD [V^2/Hz]")
            pl.legend()
            pl.show()
                 
        except:
            print "\nNo Noise PSD Data found, running algorithm now\n" 
            sig = getBufferedSignal(triggerSignal, buffer_length, plots=plots)   
            purenoiseints =  getPureNoiseIntervals(stream, sig, plots=plots)
            freqarray, J, indicesused = noisePSDrandomWelch(purenoiseints, window=win, samplesize=sample_size, sigmaCut = chiSq_sigmaCut, cutOrder=cut_order)
            NoiseChiSqHist(purenoiseints, J, indicesused, window=win,bins=15)
            with open(noiseFile, 'wb') as fp:
                pickle.dump( (freqarray, J), fp)    
    else:
        print "\nAllowing reading of saved file set to False, running Noise PSD algorithm now\n" 
        sig = getBufferedSignal(triggerSignal, buffer_length, plots=plots)   
        purenoiseints =  getPureNoiseIntervals(stream, sig,  plots=plots)
        freqarray, J, indicesused = noisePSDrandomWelch(purenoiseints, window=win, samplesize=sample_size, sigmaCut = chiSq_sigmaCut, cutOrder=cut_order)
        NoiseChiSqHist(purenoiseints, J, indicesused, window=win, bins=10)
        with open(noiseFile, 'wb') as fp:
            pickle.dump( (freqarray, J), fp) 
        
    return  freqarray, J

def time_constant_calculator(pulse_locations,stream,template=AB_theta_exp):
    
    
    
    if template==expo_temp1:
        
        mean_params=np.array([0,0])
        
        
    if template==theta_exp:
        
        mean_params=np.array([0,0,0])
    
    if template==A_theta_exp:
        
        mean_params=np.array([0,0,0,0,0])
        
    if template==AB_theta_exp:
        
        mean_params=np.array([0,0,0,0,0,0,0])
        
    k=0
    param_list=[]
    for pulse in pulse_locations:
        trial_pulse=stream[pulse[0]:pulse[1]]
        start_val=np.mean(trial_pulse[:int(0.03*fe)])
        trial_pulse=trial_pulse-start_val
        amp=np.max(trial_pulse)
        max_t=np.argmax(trial_pulse)
        max_time=max_t/fe
        ds_factor=int(len(trial_pulse)/fe)
        
        print "ds factor: ", ds_factor
        max_t_ds=int(max_t/ds_factor)
        
        trial_pulse_ds=trial_pulse[::ds_factor]
       
        m=0
        for j in range(len(trial_pulse_ds[max_t_ds:])):
            m+=1
            if trial_pulse_ds[j]<0.8*amp:
                print "breaking the loop :) " 
                break
        
        print " j is ", j
        print "m is ", m
        half_t=ds_factor*m
        halving_time=half_t/(fe)
        Tf_guess=halving_time/(np.log(1/0.8))
        
        print "my guess for the fall time for this pulse is ", Tf_guess
        
        time=np.arange(0,pulse[1]-pulse[0])/fe
        
        if template==expo_temp1:
            
            z=chisqmin((0.04,1),time,trial_pulse,np.array((pulse[1]-pulse[0])*[3]),expo_temp1)
            
        if template==theta_exp:

            z=chisqmin((0.04,1,0),time,trial_pulse,np.array((pulse[1]-pulse[0])*[3]),theta_exp,meth="Powell")
          
        errors=np.array((pulse[1]-pulse[0])*[dV_t])
        if template==A_theta_exp:

            z=chisqmin((0.65*max_time,0.005,Tf_guess,amp,0),time,trial_pulse,np.array((pulse[1]-pulse[0])*[dV_t]),A_theta_exp)
            chi_exp=z[3]
            zero=constant(time, (np.mean(trial_pulse),0))
            chi_zero=chisq(trial_pulse,zero,errors)
            
        if template==AB_theta_exp:
            
            z=chisqmin((0.1,0.005,1.1,1.2,5,40,0), time, trial_pulse,np.array((pulse[1]-pulse[0])*[dV_t]),AB_theta_exp )
            chi_exp=z[3]
            zero=constant(time,(np.mean(trial_pulse),0))
            chi_zero=chisq(trial_pulse,zero,errors)
            
        params=[]
        
        for j in range(len(z[0])):
            params.append(z[0][j])

        pl.plot(time,trial_pulse,label="trial pulse")
        pl.hold(True)
        pl.plot(time, template(time,params),label="best fit time constant result")
        pl.grid()
        pl.legend()
        pl.show()
        
        
        print "with chi squared of ", chi_exp
        print "compared to chi_zero of ", chi_zero
        
#        print "chi of exp fit and constant fit are:   ", chi_exp, " and ", chi_const
            
        params=np.array(params)
        print "these are the params ", params
        
        if template==A_theta_exp:
            
            param_tuple=(params[1],params[2])
#        if template==expo_temp1:
#            
#        
#            pl.plot(time,expo_temp1(time,params),label="best fit")
#            pl.hold(True)
#            pl.plot(time, trial_pulse, label="signal")
#            
#            pl.title("best fit time constant value plots")
#            pl.show()
#            
#            print "respective best rise and fall time ", params[0], " and ", params[1]
            
        
        if template==A_theta_exp and chi_exp<3*len(time) and chi_zero>2*len(time) and params[2]>0:
            mean_params=mean_params+params

            param_list.append(param_tuple)
            k+=1
        if template==AB_theta_exp and chi_exp<3*len(time) and chi_zero>2*len(time) and params[2]>0:
                mean_params=mean_params+params
                k+=1    
                param_list.append(params)
        
    if k!=0:
        
        mean_params=mean_params/k
        
    elif k==0 and template==A_theta_exp:
        
        mean_params=(0,0.008,1.1,40,0)
        
    elif k==0 and template==AB_theta_exp:
        
        mean_params=(0,0.008,2.4,0.5,40,20,0)
            
    return mean_params,params

def gainCalculator(OF_stream, stream, mod_pulse_locations):
    'returns gain signal to plot, and a gain list for each pulse found'
    assert len(OF_stream)==len(stream)
    g=np.zeros(len(stream))
    gains_list=[]
    for k in mod_pulse_locations:
        short_stream_OF=OF_stream[k[0]:k[1]]
        short_stream=stream[k[0]:k[1]]
        OF_max=np.max(short_stream_OF)
        stream_max=np.max(short_stream)
        gain=OF_max/stream_max
        gains_list.append(gain)
        for l in range(k[0],k[1]):
            g[l]=gain
        
        
    return g, gains_list

def noiseGenerator(J, window="boxcar"):
    'creates a randomized noise interval using the Noise PSD (no DC, high-pass filtered) '
    phases=np.zeros(len(J))
    for k in range(len(phases)):
        phase=np.random.uniform(low=0.0, high=2*np.pi)
        phases[k]=np.exp(np.complex(0,phase))
        
    J_phased=J*phases
    J_full_phased=fourier1to2sided(J_phased)
    
    S1=np.sum(sgl.get_windows(window, len(J_full_phased)))
    S2
    


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

    pl.hist(chis_red)
    pl.xlabel('Chi Sq (Reduced)')
    pl.ylabel('Events')
    
    pl.show()
    
    pl.hist(chis, range=(4e3,6e3))
    pl.xlabel('Chi Sq')
    pl.ylabel('Events')
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


    amp_hist = pl.hist(amps, range=(0,200), bins = 100)
    pl.xlabel('Amplitude [A.D.U.]')
    pl.ylabel('Events')
    pl.xlim([0,200])
    pl.show()    
    
    bin_values_list = list(amp_hist[0])
    peak_1_location = bin_values_list.index(max(bin_values_list))
    print "1st peak location", amp_hist[1][peak_1_location]
    
    pl.hist(amps, range=(0,200), bins = 100)
    pl.xlabel('Amplitude [A.D.U.]')
    pl.ylabel('Events')
    pl.yscale('log', nonposy='clip')
    pl.xlim([0,200])
    pl.show()  
    
    pl.hist(amps, range=(0,200), bins = 100)
    pl.xlabel('Amplitude [A.D.U.]')
    pl.ylabel('Events')
    pl.yscale('log', nonposy='clip')
    pl.xscale('log', nonposx='clip')
    pl.xlim([0,200])
    pl.show()   
    
    pl.hist(amps*(1.1 - 0.008), range=(0,200), bins=100)
    pl.xlabel('Integrated Pulse [A.D.U. * s]')
    pl.ylabel('Events')
    pl.xlim([0,200])
    pl.show()
    
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
                break #skip event is any of cuts not met
        if not veto:
            subset.append(result)
    return subset
            
def compareFitResults(pulse_dict_list1, pulse_dict_list2 ):
    chis1=[]
    chis2=[]
    slow_amps=[]
    fast_amps=[]
    amps=[]
    assert len(pulse_dict_list1)==len(pulse_dict_list2)
    for k in range(len(pulse_dict_list1)):
        chis1.append(pulse_dict_list1[k].get("chi")    )
        chis2.append(pulse_dict_list2[k].get("chi")    )
        slow_amps.append(pulse_dict_list2[k].get("slow amplitude") )
        fast_amps.append(pulse_dict_list2[k].get("fast amplitude") )
        amps.append(pulse_dict_list1[k].get("amplitude") )
        
    print "Number of Pulse Fit is ", len(pulse_dict_list1)
    chis1 = np.asarray(chis1)
    chis2 = np.asarray(chis2)
    deltaChi = chis2 - chis1
    pl.hist(deltaChi)
    pl.xlabel('Delta Chi Sq')
    pl.ylabel('Events')
    pl.show()

    pl.scatter(amps, deltaChi)
    pl.xlabel('Amplitude')
    pl.ylabel('Delta Chi Sq')
    pl.show()
    
def processStream(V, t, directory):
    uV=10**6*V  
    Vf = lowPassFilter(V, 500,filter_type="ideal")
    uVf=10**6*Vf
    #Vf=butterworthLPFilter(V, 100)
    VfB=butterworthLPFilter(V, 100)
    uVf=10**6*Vf
    uVfB=10**6*VfB
    uV_noDC = (removeDC(uV))
    Vs=removeDC(V)
    uVs=10**6*Vs
    
    # Settings: lag = 5, threshold = 6, influence = 0
    lag = int(5*fe)
    threshold = 6
    influence = 0
    start_timestep = int(0*fe)
    end_timestep = len(V)  
    time=np.arange(start_timestep,end_timestep)/fe
    pulsePeaks_locations = []
    pulses = []
    pulse_locations=[] #list of lists, each element contains a start and end for pulse
    window = int(0*fe) #if you want to extend pulse window
    
    stream_active=uVs[start_timestep:end_timestep]
    
    
    #stream_active=uVs[start_timestep:end_timestep]
    
    cutOfffreq = 2 #Hz
    cutOfffrac = cutOfffreq*2/fe
    b, a =sgl.butter(1,cutOfffrac,btype="highpass")
    stream_active_hp = highPassFilter(stream_active)
    
    
    
    stream_lpf=butterworthLPFilter(uVs[start_timestep:end_timestep], 100)   #butterworth filtered stream
    
    
    #f, testFT = getWindowedFFT(uVs[start_timestep:end_timestep], None)
    #f,testFT2 = getWindowedFFT(uVs[start_timestep:end_timestep], None, oneSided=False)
    #twosidedTest = fourier1to2sided(testFT, len(uVs[start_timestep:end_timestep]) )
    #diff1 = np.fft.ifft(twosidedTest) - np.asarray(uVs[start_timestep:end_timestep])
    #diff2 = np.fft.ifft(testFT2) - np.asarray(uVs[start_timestep:end_timestep])
    #pl.plot(np.linspace(start_timestep/fe, end_timestep/fe, len(uVs[start_timestep:end_timestep])   ) , diff1, label = "diff1")
    #pl.plot(np.linspace(start_timestep/fe, end_timestep/fe, len(uVs[start_timestep:end_timestep])   ) , diff2, label = "diff2")
    #pl.legend()
    #pl.show()
    
    
    (pulses, pulse_locations, pulsePeaks_locations, signal) = getTrigger(uVfB, start_timestep, end_timestep, lag, threshold, influence, read=False)
    with open(directory + 'triggerData.p', 'wb') as fp:
                pickle.dump((pulses, pulse_locations, pulsePeaks_locations, signal), fp)
    
    
    #time_k=time_constant_calculator(pulse_locations,stream_lpf,template=AB_theta_exp)
    
    #params=(0.01,0.95,0)
    #
    ##print "best params are ", time_k
    #
    #
    #
    #Tr=params[0]   #fall and rise times
    #Tf=params[1]
    #avgdelay=params[2]
    ##amp=params[3]
    #
    #print "Best rise time ", Tr
    #print "Best fall time ", Tf
          
    #sig = getBufferedSignal(signal, 4, plots = True)   
    ##z=PulseDiscriminator(uVs)
    #purenoiseints =  getPureNoiseIntervals(sig, len(signal), plots = True)
    #freqarray, J = noisePSDrandomWelch(purenoiseints,window="hann")
    #NoiseChiSqHist(purenoiseints)
    
    freqarray, J = getNoisePSD(stream_active_hp, signal, win="boxcar", buffer_length=4, sample_size="all", chiSq_sigmaCut = 6, cut_order=3, plots = True, allowRead=False)
    freq=freqarray
    with open(directory + 'NoisePSD.p', 'wb') as fp:
        pickle.dump( (freqarray, J) , fp )   
        
    
    #Compare with
    Fig2=pl.figure('Fig2')
    pl.title('Noise PSD Comparison  ')
    pl.xlabel("Frequency[Hz]")
    pl.ylabel("PSD [V^2/Hz]")
    pl.loglog(freqarray, J, label="Pulses Removed")
    pl.hold(True)
    f, Pwelch_spec = sgl.welch(uV[start_timestep:end_timestep], fe, nperseg=fe)
    pl.loglog(f, Pwelch_spec,label="With Pulses")
    pl.legend()
    pl.show()

    OF_stream, H_trans=OFTransferFn(stream_active_hp, J, freq, (0.01, 0.008,2.8,0.8,0.6,0.4,0))
    
    OF_thresh=10
    (pulses, pulse_locations, pulsePeaks_locations, signal) = getTrigger(OF_stream, start_timestep, end_timestep, lag, OF_thresh, influence,names="OF_stream", read=False)
        
    MT=MaxTrigger(signal, OF_stream, pulse_locations)
    ref_sig=MT[0]
    mod_pulse_locations=MT[1] #change this l8r to screen out NTD events with MC simulations maybe???
    
    
    buffered_ref_sig=getBufferedSignal(ref_sig)
    dV_t=timeUncertainty(stream_lpf,buffered_ref_sig)   #time_uncertainty_ estimate for time chi squared calculation
    
    gain_signal, gains_list =gainCalculator(OF_stream, stream_active, mod_pulse_locations)
    #pl.plot(time, gain_signal, label="Gain Signal")
    #pl.grid()
    #pl.legend()
    #pl.show()
    
    current_template=A_theta_exp
    useCalc=False
    
    if useCalc==True:
        time_k=time_constant_calculator(mod_pulse_locations,stream_lpf,template=current_template)
       
    
    if current_template==AB_theta_exp and useCalc==True:
        
        Tr_2amp_fit=time_k[0][1]
        Tf1_2amp_fit=time_k[0][2]
        Tf2_2amp_fit=time_k[0][3]
        
    #Tf1_fit=1.1
    #Tf2_fit=2.5
    
    if current_template==A_theta_exp and useCalc==True:
        
        Tr_fit=time_k[0][1]
        Tf_fit=time_k[0][2]    #Athetaexp, 1 amp
        
    
    #print "LIST PARAMS: ", time_k[1]
    
    
    #pl.plot(time, 100*ref_sig, label="Refined Signal w/ Max Trigger")
    pl.plot(time, stream_lpf, label="Low Pass Filt Stream")
    pl.legend()
    pl.savefig(directory + 'LPF_stream.png')
    pl.show()
    
    
    if current_template==A_theta_exp:
        
        oneAmp_processed_results=optimalFiltering1amp(stream_active,mod_pulse_locations,J,freq)
        analyzeFitResults([oneAmp_processed_results])
    with open(directory + 'oneAMP_pulse_processing_results.p', 'wb') as fp:
        pickle.dump( (oneAmp_processed_results) , fp )    
    
    #analyzeFitResults(getPulseSubset(oneAmp_processed_results, [("chi", 300,700)]))
    twoAmp_processed_results = []
    current_template=AB_theta_exp
    if current_template==AB_theta_exp:
        twoAmp_processed_results=optimalFiltering2amp(stream_active,mod_pulse_locations,J,freq)
        analyzeFitResults([twoAmp_processed_results])
        
    with open(directory + 'full_pulse_processing_results.p', 'wb') as fp:
        pickle.dump( (oneAmp_processed_results, twoAmp_processed_results) , fp )
    
    #twoAmp_processed_results=optimalFiltering2amp(stream_active,pulse_locations,J,freq)
    #print "PARAMS: ", time_k
    #with open('pulse_processing_results.p', 'wb') as fp:
    #    pickle.dump( (oneAmp_processed_results, time_k) , fp )
        
    #pl.plot(t_test, stream, label="Before Test")
    ##pl.hold(True)
    #pl.plot(t_test, stream_post_test, label="After Test")
    #pl.grid()
    #pl.legend()
    #pl.show()
    #
    #
    #
    #
    #optimal=optimalFiltering1amp(stream_active ,pulse_locations,J,f,window="boxcar",parameters=params,sglTemplateFunc=theta_exp,scan_range=int(1500))
    #
    #
    #
    #chi=optimal[0]
    #t_0=optimal[1]
    #a=optimal[2]
    #DOF=optimal[3]
    #
    #
    #
    #print "chi is ", chi
    #print "t_0 is ", t_0
    #print "mean chi is ",np.mean(chi)
    #print "mean a is ", np.mean(a)
    ##
    ##
    ##print "done filtering, reduced chi sq is ", chi/DOF
    ##print DOF
    ##print "optimal t_0 is ", t_0
    ##print "amplitude a is ", a
    #
    ##pl.xlabel('Frequency [Hz]')
    ##pl.ylabel('PSD')
    ##pl.title("Welch on whole time-series")
    ##
    ##pl.show()
    #
    #
    ###sgl.welch()
    ##f, Pwelch_spec = sgl.welch(uV, fe)
    ##pl.semilogy(f, Pwelch_spec)
    ##pl.xlabel('Frequency [Hz]')
    ##pl.ylabel('PSD')
    ##
    ##pl.show()
    #
    ##(windowed_stream, weights, N, S1, S2) = HanningWindow(uV_noDC, 0, -1)
    ##print N
    ##print S1
    ##print S2
    ##(fourier, PSD, ENBW) = PowerSpectra(windowed_stream, 0, -1, weights, N, S1, S2)
    ##print ENBW
    #
    #

    
    pl.plot(t[start_timestep:end_timestep],uVfB[start_timestep:end_timestep], label = "Low Pass Filtered")
    pl.ylabel("Amplitude from Baseline (DU)")
    pl.xlabel("Time (s)")
    pl.show()
    #
    #pl.legend( frameon=False )
    #pl.show()
    #
    ##pl.plot(t,uVf)
    ##pl.ylabel("Amplitude from Baseline (uV)")
    ##pl.xlabel("Time (s)")
    ##
    ##pl.legend("Filtered")
    ##pl.show()
    #
    #
    #
    ##pl.plot(t[10**6:2*10**6],uV[10**6:2*10**6])
    ##pl.ylabel("Amplitude (uV)")
    ##pl.xlabel("Time (s)")
    ##
    ##pl.legend("Unfiltered")
    ##pl.show()
    #
    return

def processFile(file_stream, chunk_size, start_chunk=0):
    '''
    str file_stream: is the file name containing the data
    chunk_size: number of hours to process at a time (i.e. if chunk_size is 1, break up stream into 1 hour chunks)
    Completely analyzes an entire stream of data
    start_chunk: optional int, allows you to start processing from a later chunk of the stream
    '''
    total_length_str = file_stream[-5:]
    print file_stream
    print "hi", total_length_str, total_length_str[0:2], total_length_str[3:]
    total_length_in_hours = int(total_length_str[0:2]) + (int(total_length_str[3:]))/(60.0)
    chunks = float(total_length_in_hours)/chunk_size
#    full_chunks = int(chunks)
    chunks_needed = int(math.ceil(chunks))
    
    print "Processing File", file_stream
    for i in range(start_chunk, chunks_needed):
        current_directory = 'Results/' + file_stream + '/' +str(chunk_size) + 'hours_' + str(i) + '/'
        print "Processing chunk", i , "out of", chunks_needed
        print "Current directory is", current_directory
        chunk_start_time = i*chunk_size*3600 + i*10 #add ten seconds 
        chunk_end_time = chunk_start_time + chunk_size*3600
        if i == (chunks_needed - 1): #if last chunk
            chunk_end_time = int(total_length_in_hours*3600 - 1*60) #end 1 minute before stream end
        chunk_interval_length = chunk_end_time - chunk_start_time
        
        print "Chunk time values", chunk_start_time, chunk_end_time, chunk_interval_length
#        load this section of stream
        (t, V) = loadVoltageFile(file_stream, interval_length = chunk_interval_length, start_time=int(chunk_start_time), TEST = True)
        processStream(V,t, current_directory)
    
#Begin Main Code!

#file1 = 'data_run42_dbz1\\20180314_16h12' #original neutron
#file1 = 'data_run42_dbz1\\20180314_10h11' #Ba, fails after a couple hours
#file1 = 'data_run42_dbz1\\20180313_18h32' #good Ba
#file1 = 'data_run42_dbz1\\20180315_14h24' #No source original#
#file1 = 'data_run42_dbz1\\20180315_09h55' #SANS source post calibatration neutron# fails
#print "\n\nFile: ", file1
print ""
voie = int(0)  #0-NTD1 ; 1-NTD2

Tr_fit=0.008
Tf_fit=1.1

Tr_2amp_fit=0.008
Tf1_2amp_fit=1.1
Tf2_2amp_fit=2.5

processFile(file1, chunk_size=1)
#file2 = 'data_run42_dbz1\\20180314_16h12' #original neutron
#processFile(file2, chunk_size=1)


##Decoupage du stream

#with open('pulses/no_source315/first5h/full_pulse_processing_results_first5h.p', 'rb') as fp1:
#    (oneAmp_processed_results1, twoAmp_processed_results1) = pickle.load(fp1)
#with open('pulses/no_source315/second5h/full_pulse_processing_results_second5h.p', 'rb') as fp2:
#    (oneAmp_processed_results2, twoAmp_processed_results2) = pickle.load(fp2)
##with open('pulses/no_source315/remaining_hours/full_pulse_processing_results_last.p', 'rb') as fp3:
##    (oneAmp_processed_results3, twoAmp_processed_results3) = pickle.load(fp3)
#analyzeFitResults([oneAmp_processed_results1, oneAmp_processed_results2])
#analyzeFitResults([getPulseSubset(oneAmp_processed_results1, [("chi", 300,800), ("amplitude", 0, 1000) ])
#        , getPulseSubset(oneAmp_processed_results2, [("chi", 300,800), ("amplitude", 0, 1000)  ])])           
#assert False



#load=False
#TEST=True #Preprocess downsampling ON if TRUE
#if load==False:
#    (t, V) = loadVoltageFile(file1, interval_length = dt, start_time=int(file_start_time), TEST = TEST)
#if load==True:
#    with open('DownsampledProcessedResults\\downsampled_stream_30k_s.p', 'rb') as fp:
#        (t, V) = pickle.load(fp)


#processStream(V, t, "test/")


        
    
    
        
