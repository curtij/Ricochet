#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 22:05:18 2019

@author: josephcurti
"""

import numpy as np
import os
import math
import matplotlib.pyplot as pl
import scipy.signal as sgl
import pandas as pd
import pickle
from scipy.optimize import curve_fit


Gain = 1./4600
voltageRange = 20.
ADCRange = 2**16
Gain = Gain*voltageRange/ADCRange

font =11.5
#print pl.rcParams.keys()
pl.rcParams.update({'font.size': font})
#pl.rc('text', usetex=False)
pl.rc('figure' , figsize=(10,5))
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


def process_file(file_stream, chunk_size, output_folder, start_chunk=0, 
                 end_chunk=-1, offset = 0, downsampled = True, 
                 total_duration=-1, pre_hp_filter = True):
    '''
    processes a data file by breaking it up into smaller chunks
    
    str file_stream: is the file name containing the raw data
    chunk_size: number of hours to process at a time (i.e. if chunk_size is 1, break up stream into 1 hour chunks)
    Completely analyzes an entire stream of data
    start_chunk: optional int, allows you to start processing from a later chunk of the stream
    offset: integer offset for the beginning of the data within the file
    downsample: if True downsamples the data by a factor of 10 (useful for Zinc data sampled at 10 kHz)
    total_duration: specify total length if not available from last characters of file name
    '''
    print "Processing File", file_stream
    total_length_in_hours = total_duration
    if total_duration==-1: #if not given, gets length from filename
        total_length_str = file_stream[-5:] #last characters of filename are the duration
        print file_stream
        print total_length_str, total_length_str[0:2], total_length_str[3:]
        total_length_in_hours = int(total_length_str[0:2]) + (int(total_length_str[3:]))/(60.0)
    chunks = float(total_length_in_hours)/chunk_size
    chunks_needed = int(math.ceil(chunks))
    if end_chunk > 0:
        chunks_needed = end_chunk
    
    for i in range(start_chunk, chunks_needed):
        current_directory = output_folder + file_stream + '/' +str(chunk_size) + 'hours_' + str(i) + '/'
        print "Processing chunk", (i + 1) , "out of", chunks_needed
        print "Current directory is", current_directory
        chunk_start_time = i*chunk_size*3600 + i*10 #add ten seconds 
        chunk_end_time = chunk_start_time + chunk_size*3600
        if i == (chunks_needed - 1) and end_chunk ==  -1: #if last chunk
            chunk_end_time = int(total_length_in_hours*3600 - 1*60) #end 1 minute before stream end
        chunk_interval_length = chunk_end_time - chunk_start_time
        
        print "Chunk time values", chunk_start_time, chunk_end_time, chunk_interval_length
        
        (V_raw, t, V_hpf) = load_voltage_file(file_stream, interval_length = chunk_interval_length,
                                    file_offset = offset ,start_time=int(chunk_start_time)
                                    , downsample = downsampled, preHPFilter = pre_hp_filter)
        
        process_stream(V_hpf, t, current_directory, V_raw)
        
def load_voltage_file(fname, interval_length, file_offset = 0, start_time=0, downsample=True, preHPFilter = True):
    '''
    start_time is in seconds
    interval length is in seconds, formerly 'dt'
    returns raw stream, times, hpf stream
    
    if taking too long, don't build t
    '''
    pathtofile = 'DataFiles/'+fname +'.BIN0'
    if fname=='Germanium/se25g000_000':
        pathtofile = 'DataFiles/' + fname
    here = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(here, pathtofile)
    print filename
       
    if downsample: #downsamples by factor of 10
        voltages = np.memmap(filename, dtype='int16', mode='r+', offset=file_offset, shape=(int(10*fe*start_time)+int(10*fe*interval_length),1))
        voltages=voltages[(int(10*fe*start_time)):]
        
        print "fe interval_length", fe, interval_length, int(10*fe*interval_length)
        V=voltages[:int(10*fe*interval_length),0]*Gain
#        print len(t), len(V)
#        t = np.arange(len(V)) #in timestep units
#        pl.plot(t/(10*fe), (10**6)*V, label="pre downsampled")
#        pl.xlabel("Timesteps")
#        pl.legend()
#        pl.show()      
#        freqs, Vpsd = signal.periodogram(10**6*V, 10*fe)
#        pl.loglog(freqs, Vpsd, label='pre downsampled PSD') 
#        pl.legend()        
#        pl.show()
        
#        t = t[::10]
        V = sgl.decimate(V, 10) #downsamples by 10x
        V = 10**6*V 

        t = np.arange(len(V)) #in timestep units
#        pl.plot(t/fe, V, label='scipy downsampled ') 
        pl.plot(t/fe, V, label='Raw Data') 
        pl.xlabel('Time [s]')
        pl.ylabel('Amplitude [ADU]')
        pl.legend()        
        pl.show()   
           
#        freqs2, Vpsd2 = signal.periodogram(10**6*V, fe)        
#        pl.loglog(freqs, Vpsd, label='pre downsampled PSD') 
#        pl.loglog(freqs2, Vpsd2, label='scipy downsampled PSD') 
#        pl.legend()        
#        pl.show()  
          
    else: #then set fe back to 10000 at top for Zinc
        voltages = np.memmap(filename, dtype='int16', mode='r+', offset=file_offset, shape=(int(fe*start_time)+int(fe*interval_length),1))
        voltages = voltages[(int(fe*start_time)):]
        t = np.arange(len(voltages))
        V = voltages[:int(fe*interval_length),0]*(Gain*10**4)
        pl.plot(t/fe, V, label='raw data') 
        pl.legend()        
        pl.show() 
        
    V_raw = V
    V_hpf = []
    if preHPFilter:
        V_hpf = high_pass_filter(V)
    
    return (V_raw, (t/fe), V_hpf)

def process_stream(stream_active_hp, t, directory, raw_stream):
    '''
    Processes a data stream:
    1. Calculates the Noise PSD
    2. Optimally filters the data stream
    3. Runs a trigger to find pulses in the optimally filtered data
    4. Fits the found pulses
    5. Returns the found pulses and their fit information
    '''
    
    pl.plot(t, stream_active_hp, label='pre high pass filtered data')
    pl.legend()
    pl.show()
    
    V_lpf = low_pass_filter(raw_stream, 100, 1)
    pl.plot(t, V_lpf, label='Low Pass Filter applied to raw stream')
    pl.legend()
    pl.show()
    
    freq_array, J, noise_events_df = get_noise_PSD_new(stream_active_hp)
    with open(directory + 'NoisePSD_and_DF.p', 'w+') as fp:
        pickle.dump( (freq_array, J, noise_events_df) , fp )
        
    template_func = theta_exp_pulse
    
    rise = 8e-3
    fall = 1.1
    standard_params = (PULSE_START, rise, fall)
    if target=="Germanium":
        rise = 6e-3
        fall = 16.3e-3
        standard_params = (PULSE_START, rise, fall)

    
    OF_stream, H_trans, h = OFTransferFn(stream_active_hp, J, freq_array, standard_params,  template_func)
    trigger_df = absolute_trigger(stream_active_hp, 50*np.sqrt(h))
    raw_pulses=[]
    hpf_pulses=[]
    for index, row in trigger_df.iterrows():
        timestep_of_peak = row['timestep_of_peak']
        start_timestep = int(timestep_of_peak) - int(PULSE_START*fe)
        end_timestep = start_timestep + int(PULSE_WINDOW*fe)
        raw_pulse = raw_stream[start_timestep:end_timestep]
        hpf_pulse = stream_active_hp[start_timestep:end_timestep]
        raw_pulses.append(raw_pulse)
        hpf_pulses.append(hpf_pulse)
        
#        if index%10==0:
#            pl.plot(np.linspace(start_timestep/fe,end_timestep/fe,int(PULSE_WINDOW*fe)), raw_pulse, label='pulse'+str(index) )
#    #        pl.plot(np.linspace(start_timestep/fe,end_timestep/fe,int(PULSE_WINDOW*fe)), hpf_pulse, label='hpf pulse'+str(index) )
#            pl.legend()
#            pl.show()
#        
    pulse_df = trigger_df
    pulse_df['raw_pulse'] = raw_pulses
    pulse_df['hpf_pulse'] = hpf_pulses

    print len(pulse_df.index), "Pulses found"
#    for index, pulse in pulse_df.iterrows():
#        timestep_of_peak = row['timestep_of_peak']
#        start_timestep = int(timestep_of_peak) - int(PULSE_START*fe)
#        end_timestep = start_timestep + int(PULSE_WINDOW*fe)
#        time_pulse = np.linspace(start_timestep/fe, end_timestep/fe, int(PULSE_WINDOW*fe))
#        pl.plot(time_pulse, pulse['raw_pulse'], label = 'raw pulse')
#        pl.plot(time_pulse, pulse['hpf_pulse'], label = 'hpf pulse')
#        pl.legend()
#        pl.show()
    
    result_og = fit_1amp(pulse_df, J, freq_array, standard_params)
    if target=="Zinc":
                    #t_0, Tr1, Tf1, Tr2, Tf2, ratio
        Fe_params = (PULSE_START, 0.00336984, 1.57837824, 0.02332444, 0.35925711, 0.41811186 )
        Ba_params = (PULSE_START, 0.0043958, 2.08564392, 0.01805898, 0.42211096, 0.71273393 )
        AmBe_params = (PULSE_START, 0.00456286, 2.08542989, 0.01883717, 0.43836557, 0.70508009)
        
        result_Fe_Temp = fit_1amp(pulse_df, J, freq_array, Fe_params, temp_1 = two_theta_exp_oneamp)
        result_Ba_Temp = fit_1amp(pulse_df, J, freq_array, Ba_params, temp_1 = two_theta_exp_oneamp)
        result_AmBe_Temp = fit_1amp(pulse_df, J, freq_array, AmBe_params, temp_1 = two_theta_exp_oneamp)
        
        t_plot = np.linspace(0, PULSE_WINDOW, int(PULSE_WINDOW*fe))
        pl.plot(t_plot, two_theta_exp_oneamp(t_plot, *Fe_params), label = 'Fe based temp')
        pl.plot(t_plot, two_theta_exp_oneamp(t_plot, *Ba_params), label = 'Ba based temp')
        pl.plot(t_plot, two_theta_exp_oneamp(t_plot, *AmBe_params), label = 'AmBe based temp')
        pl.legend()
        pl.show()
        
        t_plot = np.linspace(0, 2*PULSE_WINDOW, int(2*PULSE_WINDOW*fe))
        pl.plot(t_plot, two_theta_exp_oneamp(t_plot, *Fe_params), label = 'Fe based temp')
        pl.plot(t_plot, two_theta_exp_oneamp(t_plot, *Ba_params), label = 'Ba based temp')
        pl.plot(t_plot, two_theta_exp_oneamp(t_plot, *AmBe_params), label = 'AmBe based temp')
        pl.legend()
        pl.show()

        with open(directory + 'oneAMP_' + str(int(10*PULSE_WINDOW)) + 'ds_window_results_df.p', 'w+') as fp:
            pickle.dump( (result_og, result_Fe_Temp, result_Ba_Temp, result_AmBe_Temp) , fp )  
       
    if target=="Germanium":
        test_pulse = result_og["raw_pulse"].iloc[-1]
        t_plot = np.linspace(0, PULSE_WINDOW, int(PULSE_WINDOW*fe))
        initial_guess = (PULSE_START, 17, 6e-3, 16e-3, -0.1)
        def func(t, t0, a, R, F, b):
            return a*theta_exp_pulse(t, t0, R, F) + b
#        bound = ()
        fit = curve_fit(func, t_plot, test_pulse, initial_guess )
        pl.plot(t_plot, test_pulse, 'b', label = 'pulse')
        pl.plot(t_plot, func(t_plot, *fit[0]), 'r', label = 'Fit')
        pl.legend()
        pl.show()
        print fit
        print fit[0]
        
        bound = ((PULSE_START-50e-3, PULSE_START+50e-3), (3,100), (0, 1), (0,1), (-0.05, 0.05))
        fit = curve_fit(func, t_plot, test_pulse, initial_guess, bounds=bound )
        pl.plot(t_plot, test_pulse, 'b', label = 'pulse')
        pl.plot(t_plot, func(t_plot, *fit[0]), 'r', label = 'bound Fit')
        pl.legend()
        pl.show()
        print fit
        print fit[0]
        with open(directory + 'oneAMP_' + str(int(10*PULSE_WINDOW)) + 'ds_window_results_df.p', 'w+') as fp:
            pickle.dump( (result_og) , fp )  
       

        
    return




def high_pass_filter(timeseries, filtersize=-1, order = 2, inverse=False):
        '''
        Butterworth High Pass Filter
        filter size in Hz
        '''
        if filtersize == -1:
            filtersize = 2.0/PULSE_WINDOW
    
        b, a =sgl.butter(order,2*filtersize/fe,btype="highpass")
        if inverse==False:
            z=sgl.lfilter(b,a,timeseries)
        else:            
            z=sgl.lfilter(a,b,timeseries)
        return z
    
def low_pass_filter(timeseries, filtersize, order, inverse=False):
        '''
        Butterworth Low Pass Filter
        filter size in Hz
        '''
        b, a =sgl.butter(order,2*filtersize/fe,btype="lowpass")
        if inverse==False:
            z=sgl.lfilter(b,a,timeseries)
        else:
            z=sgl.lfilter(a,b,timeseries) 
        return z
    
def getWindowedFFT(x, scaling="density", window = "boxcar", removeDC = True, oneSided = True):
    '''
    set scaling to "density" for comparison with LPSD
    returns frequencies, FFT
    '''
    win_vals = sgl.get_window(window, len(x))
    windowed_x = win_vals*x
    
    FTfreq=np.arange(0,len(x))*fe/len(x)
    FTfreq=FTfreq[0:int(1+len(FTfreq)/2)]
    FT=np.fft.rfft(windowed_x)
    
    if scaling == "density":
        S1 = np.sum(win_vals)
        S2 = np.sum(win_vals**2)
        ENBW = fe*S2/(float(S1**2))
        FT = FT/(S1*np.sqrt(ENBW))
#    if oneSided:
#        FTfreq=FTfreq[0:int(1+len(FTfreq)/2)]
#        FT=FT[0:int(1+len(FT)/2)]
#        for j in range(1, 1+len(FT)/2):
#            FT[j]=2*FT[j]
    if removeDC:
        FT[0]=0
            
    return FTfreq, FT

def getInverseFFT(x, scaling="density", window = "boxcar", oneSided = True):
    '''
    unscales: scaling to "density" for comparison with LPSD
    returns time domain
    '''
    win_vals = sgl.get_window(window, 2*len(x))
#    windowed_x = win_vals*x
    inverseFT=np.fft.irfft(x)
    if scaling == "density":
        S1 = np.sum(win_vals)
        S2 = np.sum(win_vals**2)
        ENBW = fe*S2/(float(S1**2))
        inverseFT = inverseFT*(S1*np.sqrt(ENBW))

    return inverseFT

def step(x):
    return 1 * (x > 0) 

def expo_temp1(t, (Tr,Tf)):
    return (-np.exp(-t/Tr) + np.exp(-t/Tf)) 

def theta_exp_pulse(t, t_0, Tr, Tf):
    shifted_t = np.asarray(t) - t_0    
    m=expo_temp1(shifted_t, (Tr,Tf))
    
    return step(shifted_t)*m 

def two_theta_exp_oneamp(t, t_0, Tr1, Tf1, Tr2, Tf2, ratio):
    return theta_exp_pulse(t, t_0, Tr1, Tf1) + ratio*theta_exp_pulse(t, t_0, Tr2, Tf2)
    
def get_noise_PSD_new(stream_hpf, method = "uniform_sample", rms_threshold = -1):
    '''
    Calculates the noise PSD 
    
    f: array of sample frequencies
    J: noise psd
    noise_dataframe is a dataframe that has noise interval start time, end time,
                    raw interval, pre-hp filtered interval, fit with good template,
                    fit with no template
    '''
    if method == "uniform_sample":
        total_timesteps = len(stream_hpf)
        stream_end_time = total_timesteps/fe
         
        sample_fraction = 0.5 #0.15
        sample_size = int(sample_fraction*stream_end_time/PULSE_WINDOW)
        print sample_size
        
        np.random.seed(1996)
        sample_start_times = (stream_end_time - 4.0*PULSE_WINDOW)*np.random.sample(sample_size)\
                                + 2.0*PULSE_WINDOW
        sorted_sample_start_times = np.sort(sample_start_times)
        print sorted_sample_start_times
        
        sample_start_timesteps = fe*sorted_sample_start_times
        sample_start_timesteps = (sample_start_timesteps.astype(int))
        sample_end_timesteps = sample_start_timesteps + int(fe*PULSE_WINDOW)
                
        periodograms = []
        hpf_intervals = []
        rms_amplitudes = []
        for i in range(len(sample_start_timesteps)):
            hpf_interval = stream_hpf[sample_start_timesteps[i]:sample_end_timesteps[i]]
            f, periodogram = sgl.periodogram(hpf_interval, fe)
            
#            onesec = np.linspace(0,1,int(PULSE_WINDOW*fe))
#            pl.plot(onesec, hpf_interval, label = 'noise interval '+str(i))
#            pl.legend()
#            pl.show()
            
            rms_amp = np.sqrt(np.mean((hpf_interval**2)))
            rms_amplitudes.append(rms_amp)
            hpf_intervals.append(hpf_interval)
            periodograms.append(periodogram)
            
        initial_noise_dataframe = pd.DataFrame(data={'start timestep': sample_start_timesteps, 
                                                'end timestep': sample_end_timesteps,
                                                'periodogram': periodograms,
                                                'HPF Interval': hpf_intervals,
                                                'RMS_amplitude': rms_amplitudes})
        cut_order = 10
        sigma_cut = 5
        chi_sq_mean = fe*PULSE_WINDOW/2.0
        chi_sq_variance = 2*chi_sq_mean
        chi_sq_std =np.sqrt(chi_sq_variance)
        chi_threshold = chi_sq_mean + chi_sq_std*sigma_cut
        
        print "chi mean, std, threshold", chi_sq_mean, chi_sq_std, chi_threshold
        noise_PSD_list = []
        
        #Initial rms amplitude cut
        pl.hist(initial_noise_dataframe['RMS_amplitude'], bins = 30)
        pl.ylabel('Noise Events')
        pl.xlabel('RMS Amplitude')
        pl.show()
        print initial_noise_dataframe['RMS_amplitude'].describe()
        if rms_threshold==-1:
            #rms_threshold = 2*10**0
            rms_threshold = 3
            if target=="Zinc":
                rms_threshold = 30
        initial_noise_dataframe = initial_noise_dataframe[initial_noise_dataframe.RMS_amplitude < rms_threshold]
        
#        noise_events = list(initial_noise_dataframe['HPF Interval'])
#        t_interval = np.linspace(0,PULSE_WINDOW, int(fe*PULSE_WINDOW))
#        for k in noise_events:
#            pl.plot(t_interval, k)
#            pl.xlabel('Time [s]')
#            pl.show()
        
        for j in range(cut_order):
            noise_PSD_temp = initial_noise_dataframe['periodogram'].mean()
#            pl.loglog(f, noise_PSD_temp, label = "noise psd " +  str(j))
#            pl.xlabel('Frequency [Hz]')
#            pl.legend()
#            pl.show()
#            avg_noise_interval = initial_noise_dataframe['HPF Interval'].mean()
#            pl.plot(np.arange(len(avg_noise_interval))/fe, avg_noise_interval)
#            pl.show()
            noise_PSD_list.append(noise_PSD_temp)
            pre_cut_length = len(initial_noise_dataframe.index) 
            if j>0  : #drop intervals with too large chi-squared
                initial_noise_dataframe = initial_noise_dataframe[initial_noise_dataframe.chi < chi_threshold]
                if pre_cut_length == len(initial_noise_dataframe.index):
#                    print "removing max"
                    initial_noise_dataframe = initial_noise_dataframe[initial_noise_dataframe.chi < max(initial_noise_dataframe["chi"])]

            initial_noise_dataframe["chi"] = initial_noise_dataframe.apply(
                    lambda row: sum( row.periodogram[1:]/noise_PSD_temp[1:] ), axis=1)
            print "\n", initial_noise_dataframe['chi'].describe()
            pl.loglog(f, noise_PSD_list[j], label = "noise psd " +  str(j) + 
                      " sample size " + str(len(initial_noise_dataframe.index)))
            
            if (initial_noise_dataframe["chi"].mean() == chi_sq_mean and 
                initial_noise_dataframe["chi"].std()== chi_sq_std):
                print "Stopping recursive cut since mean chi sq is good"
                break
        pl.legend()
        if target=="Germanium":
            pl.ylim([10**-10,10**-2]) #Ge
        if target=="Zinc":
            pl.ylim([10**-6,10**2]) #Zn
        pl.show()
        
        pl.hist(initial_noise_dataframe["chi"], bins = 30)
        pl.xlabel('Chi Sq')
        pl.ylabel('Noise Events')
        pl.show()
        pl.hist(initial_noise_dataframe["chi"], bins = 30, range=(0,fe))
        pl.xlabel('Chi Sq')
        pl.ylabel('Noise Events')
        pl.show()
        pl.hist(initial_noise_dataframe["chi"], bins = 50, range=(chi_sq_mean-5*chi_sq_std,chi_sq_mean+5*chi_sq_std))
        pl.xlabel('Chi Sq')
        pl.ylabel('Noise Events')
        pl.show()
        print "\n", initial_noise_dataframe['chi'].describe()

#        for j in range(len(noise_PSD_list)):
#            pl.loglog(f, noise_PSD_list[j], label = "noise psd " +  str(j))
#            pl.legend()
#        pl.ylim([10**-1,10**5])
#        pl.show()
            
        
        J = noise_PSD_list[-1]
        noise_dataframe = initial_noise_dataframe
        print noise_dataframe.describe()
            
    
#    sig = getBufferedSignal(triggerSignal, buffer_length, plots=plots)   
#    purenoiseints =  getPureNoiseIntervals(stream, sig,  plots=plots)
#    freqarray, J, indicesused = noisePSDrandomWelch(purenoiseints, window=win, samplesize=sample_size, sigmaCut = chiSq_sigmaCut, cutOrder=cut_order)
#    NoiseChiSqHist(purenoiseints, J, indicesused, window=win, bins=10)
#    with open(noiseFile, 'wb') as fp:
#        pickle.dump( (freqarray, J), fp) 
        
    return  f, J, noise_dataframe

def OFTransferFn(stream, J, freq , parameters, template_function):
    '''
    stream should be unfiltered
    Computes the Optimal Filter transfer function based on the supplied template
    Returns the optimally filtered stream, the OF transfer function, and the resolution
    '''
    
    time=np.arange(0,len(stream))/fe
#    b, a =sgl.butter(1,0.0004,btype="highpass")

    time_window = np.linspace(0, PULSE_WINDOW, int(fe*PULSE_WINDOW))
    s = template_function(time_window, *parameters)
    s_highpassfiltered_timedomain = high_pass_filter(s)
 
    #old way
    win_vals = sgl.get_window("boxcar", len(s))
    S1 = np.sum(win_vals)
    S2 = np.sum(win_vals**2)
    ENBW = fe*S2/(float(S1**2))
        
    s_FFT_twoHz_twosided = np.fft.fft(s_highpassfiltered_timedomain)
    s_FFT_twoHz_twosided[0]=0
    s_FFT_twoHz_onesided = s_FFT_twoHz_twosided[0:1+len(s_FFT_twoHz_twosided)/2]/(S1*np.sqrt(ENBW))  

    oneovrh=2*np.sum((np.abs(s_FFT_twoHz_onesided)**2)/J)
    h=1.0/oneovrh
    print "h old method", h
    #end old way
    
    
    f, s_hpf_freqdomain = getWindowedFFT(s_highpassfiltered_timedomain, scaling="density")
    s_hpf_freqdomain[0]=0
    one_over_h =  np.sum( (np.abs(s_hpf_freqdomain)**2)/J ) 
    h = (1.0/6.49)*1.0/one_over_h
    
    
    print "\nResolution after OF (sqrt(h)): " , np.sqrt(h)
    print "Resolution before OF: " , np.sqrt(2*np.sum(J))
    
    
    Nf=(h**2)*(np.abs(s_hpf_freqdomain)**2)/J
    pl.loglog(freq, J, label='Original Noise PSD')
    pl.loglog(freq, Nf, label='Optimal Filtered Noise PSD')
    pl.legend()
    pl.show()
    
#    t_max=2.05*np.pi*t_delay #time position of template max
    t_max=1.025*PULSE_START #desired time position of template max within pulse window
    phases = []
    for i in range(len(freq)):
        phases.append(np.complex(0,-freq[i]*2*np.pi*t_max))  
    phases=np.array(phases)
    print "type phases", type(phases), type(phases[0])
    
    H=h*np.conjugate(s_hpf_freqdomain)*np.exp(phases)/J #one sided transfer function

    pl.loglog(freq, np.abs(s_hpf_freqdomain), label="Template LPSD after hpf" )
    pl.loglog(freq, np.abs(s_hpf_freqdomain*H), label="OF filtered template LPSD" )
    pl.legend()
    pl.show()
            
    pl.loglog(freq, np.abs(H), label="Transfer Fn Magnitude" )
    pl.legend()
    pl.show()
    
    pl.loglog(freq, np.abs(s_hpf_freqdomain), label='abs s')
    pl.loglog(freq, np.sqrt(J), label='LPSD')
    pl.loglog(freq, np.abs(H), label="Transfer Fn Magnitude" )
    pl.legend()
    pl.show()
    
#    Htwoside=fourier1to2sided(H, N=int(fe))
#    
#    transferred_template_FFT=[]
#    
#    
#    transferred_template_FFT=Htwoside*s_FFT_twoHz_twosided
    
    
    OF_filtered_template_freqdomain = H*s_hpf_freqdomain

    OF_filtered_template_timedomain = getInverseFFT(OF_filtered_template_freqdomain, scaling="density")

    
    pl.plot(time_window, OF_filtered_template_timedomain, label="OF template")
    pl.plot(time_window, s_highpassfiltered_timedomain, label="high pass filtered template" )
    pl.plot(time_window, s , label="original template" )
    pl.legend()
    pl.show()


    even_seconds = len(stream)%int(2*fe)==0
    assert even_seconds 

    
#    h_t=np.fft.ifft(Htwoside) #impulse response function
    
    h_t = np.fft.irfft(H) #impulse response function of length M
    
    h_tfirsthalf=h_t[:int(len(h_t)/2)] 
    h_tsecondhalf=h_t[int(len(h_t)/2):] 
    w_h_tfirsthalf=np.cos(np.linspace(0, np.pi/2.0, len(h_tfirsthalf)))*h_tfirsthalf
    w_h_tsecondhalf=np.cos(np.linspace(3*np.pi/2.0, 2*np.pi, len(h_tsecondhalf) ) )*h_tsecondhalf
    zeros=np.zeros(len(h_t))
    
    h_padded=np.concatenate( (w_h_tfirsthalf, zeros, w_h_tsecondhalf),0)

#    H_padded=np.fft.fft(h_padded)
#    H_pad_onesided=H_padded[:int(1+len(H_padded)/2)] #transfer function, length 10k, M=10,000=1s
    H_pad_onesided = np.fft.rfft(h_padded)
    
    stream_segments1=np.split(stream, int(len(stream)/(2*PULSE_WINDOW*fe)))
    stream_segments2=np.split(stream[int(1*PULSE_WINDOW*fe):], np.arange(int(2*PULSE_WINDOW*fe), len(stream[int(1*PULSE_WINDOW*fe):]), int(2*PULSE_WINDOW*fe) ))
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
#        segment_fft=np.fft.fft(j)
#        segment_fft_oneside=segment_fft[:len(segment_fft)/2+1]
        
        segment_fft_oneside = np.fft.rfft(j)
        OF_fft = H_pad_onesided*segment_fft_oneside
        
#        OF_fft_twoside=fourier1to2sided(OF_fft, len(j))
        
#        OF_stream1=np.fft.ifft(OF_fft_twoside)
#        OF_stream1=np.real(OF_stream1)
        OF_stream1 = np.fft.irfft(OF_fft)
        OF_stream_segments1.append(OF_stream1)   #list of filtered 2s intervals
        
    
    OF_stream_segments2=[]
    for j in stream_segments2:
#        segment_fft=np.fft.fft(j)
#        segment_fft_oneside=segment_fft[:len(segment_fft)/2+1]
        
        segment_fft_oneside = np.fft.rfft(j)
        OF_fft = H_pad_onesided*segment_fft_oneside
        OF_stream2 = np.fft.irfft(OF_fft)
        OF_stream_segments2.append(OF_stream2)
        
#        OF_fft=H_pad_onesided*segment_fft_oneside
        
#        OF_fft_twoside=fourier1to2sided(OF_fft, len(j))
#        
#        OF_stream2=np.fft.ifft(OF_fft_twoside)
#        OF_stream2=np.real(OF_stream2)
#        
#        OF_stream_segments2.append(OF_stream2)

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
    
    pl.plot(np.arange(36*fe)/fe, stream[:int(36*fe)], label="pre hpf Stream")  
    pl.legend()
    pl.show()
    
    pl.plot(np.arange(36*fe)/fe, flat_OF_fullstream[:int(36*fe)], label="Optimally Filtered Stream")  
    pl.legend()
    pl.show()
        
    pl.plot(time, flat_OF_fullstream, label="Optimally Filtered Stream")  
    pl.legend()
    pl.show()

    pl.plot(time, stream, label="pre hpf Stream")
    pl.legend()
    pl.show()
    
    return flat_OF_fullstream, H_pad_onesided, h


def absolute_trigger(stream, threshold):
    '''
    Trigger g
    '''
    dead_time = 0.5*PULSE_WINDOW
#    dead_time = 1.0*PULSE_WINDOW
    dead_timesteps = int(dead_time*fe)
    print "dead time", dead_time, dead_timesteps
    print "threshold", threshold
    
    total_time = len(stream)/fe
#    stream = stream[:int(64*fe)] #remove this later
    
    time = np.arange(len(stream))/fe
    mask = np.diff(1 * (stream > threshold) != 0)
    pl.plot(time, stream)
    pl.plot(time[:-1][mask], stream[:-1][mask], 'go', markersize=5)
    pl.show()

    try:
        pulse_intervals = np.split(time[:-1][mask], len(time[:-1][mask])/2 )
    except:
        print "unequal division, removing first element"
        pulse_intervals = np.split( time[:-1][mask][1:], len(time[:-1][mask][1:])/2 )

    pulse_peaks = []
    pulse_peak_timesteps = []
    for i in pulse_intervals:
        start_timestep = int(i[0]*fe)
        end_timestep = int(i[1]*fe)+1
        if ((start_timestep/fe < 2*PULSE_WINDOW) or (
                (end_timestep/fe > total_time - 2*PULSE_WINDOW))):
            #skip events within 2 pulse windows of start or end of stream
            continue
        
        interval = list(stream[start_timestep:end_timestep])
        peak = max(interval)
        peak_timestep = int((interval.index(peak)) + start_timestep)
        pulse_peaks.append(peak)
        pulse_peak_timesteps.append(peak_timestep)


    pulse_peak_times =  np.array(pulse_peak_timesteps)/fe #in seconds
    
    pulse_df = pd.DataFrame(data={'timestep_of_peak': pulse_peak_timesteps, 
                                  'peak_value': pulse_peaks,
                                  'time_of_peak': pulse_peak_times
                                  })
    print pulse_df.describe()
    #now to impose dead time
    sorted_by_peaks = pulse_df.sort_values(by='peak_value', ascending=False)
    indices_to_remove = []
    
    inital_length = len(sorted_by_peaks.index)
    for index, row in sorted_by_peaks.iterrows():
#        print index
#        print row
        
        if index in indices_to_remove:
            #if pulse already selected for removal: skip
#            print "skipping index", index, "because already selected for removal"
            continue
        
        inside_dead_window = True
        n = 1
        while inside_dead_window and index+n < inital_length:
            if ((pulse_df.loc[index + n]['timestep_of_peak'] - row['timestep_of_peak'])
                    <= dead_timesteps):
#                print pulse_df.loc[index + n]['timestep_of_peak'], row['timestep_of_peak']
#                print "Removing index", index+n, pulse_df.loc[index + n]['time_of_peak']
#                print "because it is within dead time ahead of index", index,  row['time_of_peak']
#                print 
                indices_to_remove.append(index+n)
            else:
                inside_dead_window = False
            n += 1
            
        inside_dead_window = True
        n = 1
        while inside_dead_window and index - n >= 0:
            #Check pulses behind to see if they are within dead time of larger pulse
            if ( (row['timestep_of_peak'] - pulse_df.loc[index - n]['timestep_of_peak'] )
                    <= dead_timesteps):
                indices_to_remove.append(index-n)
#                print row['timestep_of_peak'] - pulse_df.loc[index - n]['timestep_of_peak']
#                print "Removing index", index-n, pulse_df.loc[index - n]['time_of_peak']
#                print "because it is within dead time behind of index", index, row['time_of_peak']
#                print
            else:
                inside_dead_window = False
            n += 1

#    print "indices to remove", indices_to_remove
    processed_pulse_df = pulse_df.drop(indices_to_remove)
    
    pl.plot(time, stream)
    pl.plot(time[:-1][mask], stream[:-1][mask], 'go', markersize=5)
    pl.plot(processed_pulse_df['time_of_peak'], processed_pulse_df['peak_value'], 'ro', markersize=5)
    pl.show()

    return processed_pulse_df 

def fit_1amp(pulse_df, J, freq, parameters, temp_1=theta_exp_pulse, t0=0.0, scan_range = 0.04, window = "boxcar"):
    #t0 was 0.01, scan range 0.35
    
    def chisquared(V, expected, a, J):
        return sum((np.abs(V[1:] - a*expected[1:] ))**2 / J[1:])   
    
    min_chis=[]
    min_t0s=[]
    min_amps=[]
#    pulses={}
    templates=[]
#    min_A={}
#    pulse_time={}
    chi_zeros=[]
    amp_errors=[]
    
    #results_list = []
#    time_window = np.linspace(0,PULSE_WINDOW, int(PULSE_WINDOW*fe))
#    template = temp_1(time_window, *parameters)
#    freqs, S = getWindowedFFT(template)
#    S[0] = 0
#            
    time_window = np.linspace(0,PULSE_WINDOW, int(PULSE_WINDOW*fe))

    t0s= np.concatenate( [np.arange(t0, t0+scan_range, 0.001), np.arange(t0 - scan_range, t0, 0.001)] )
   
    S_t0s = []
    for i in range(len(t0s)):
        time_window = np.linspace(0,PULSE_WINDOW, int(PULSE_WINDOW*fe))
        template = temp_1(time_window, parameters[0] + t0s[i] , *(parameters[1:]) )
        template = high_pass_filter(template)
        freqs, S = getWindowedFFT(template)
        S[0] = 0
        S_t0s.append(S)

    for index, pulse in pulse_df.iterrows():   
        hpf_pulse = pulse['hpf_pulse']
        freqs, V = getWindowedFFT(hpf_pulse)
        V[0]=0
        
        #made the t0s like this because t0 ~= 0 should be close to minimum
        for j in range(len(t0s)):
            t0 = t0s[j]
#            phases = (-2*np.pi*freq*t0)*np.complex256(1j)
#            S_t0 = S*np.exp(phases)
          
            S_t0 = S_t0s[j]
       
            numer = 0
            denom = 0            
            numer = sum(( np.real( np.conjugate(S_t0) * V )/J ))
            denom = sum(np.abs(S_t0)**2/J)
            ahat_t = float(numer)/float(denom)
                        
            
#            template_t0= ahat_t*temp

            
            chi=chisquared(V, S_t0, ahat_t, J)
#            print "best theta ", best_theta, best_theta[0]
#            a_f=best_theta[0][0]
#            a_s=best_theta[0][1]
            chi_zero=chisquared(V, S_t0, 0, J)

            if t0==t0s[0]:
                min_chi=chi
                best_amp = ahat_t
                best_t0 = t0
                best_template = ahat_t*temp_1(time_window, PULSE_START+best_t0, *parameters[1:])
                best_chizero = chi_zero
                amp_error = denom
                
            if chi<=min_chi:
                min_chi=chi
                best_amp = ahat_t
                best_t0 = t0
                best_template = ahat_t*temp_1(time_window, PULSE_START+best_t0, *parameters[1:])
                best_chizero = chi_zero
                amp_error = denom

                
        min_chis.append(min_chi)
        min_amps.append(best_amp)
        min_t0s.append(best_t0)
        templates.append( best_template  )
        chi_zeros.append(best_chizero)
        amp_errors.append(amp_error) #algebraic result for error 

#    fit_results = {'chi': min_chis, 'amplitude': min_amps, 'delay': min_t0s, 
#                   'template': templates, 'chi_zero': chi_zeros, 'amp_error': amp_errors}
#    fit_results_df = pd.DataFrame(fit_results)
    
#    result_df = pd.concat([pulse_df, fit_results_df], axis=1)
#    result_df2 = pulse_df.assign(chi = min_chis)
    result_df = pulse_df.assign(chi = min_chis, amplitude = min_amps, delay = min_t0s, 
                   template = templates, chi_zero = chi_zeros, amp_error = amp_errors)


    counter = 0
        
    print result_df['timestep_of_peak'].describe()
    for j, row in result_df.iterrows():
#        if counter==0:
#            fit = curve_fit(A_theta_exp2, np.arange(fe)/fe, pulses[0][0:int(fe)], (0.24, 0.008,0.04, 350, 0))
#            print fit
#            pl.plot(np.arange(fe)/fe, pulses[0][0:int(fe)])
#            pl.plot(np.arange(fe)/fe, A_theta_exp2(np.arange(fe)/fe, *fit[0]))
#            pl.show()
            
        if counter%10==0:
            timestep_of_peak = row['timestep_of_peak']
#        try:
            start_timestep = int(timestep_of_peak) - int(PULSE_START*fe)
            end_timestep = start_timestep + int(PULSE_WINDOW*fe)
            times = np.linspace(start_timestep/fe, end_timestep/fe, int(PULSE_WINDOW*fe))
            
            pl.plot(times, row['raw_pulse'], label="Observed Pulse", color='blue')
            pl.xlabel("Time [s]")
            pl.ylabel("Amplitude [ADU]")
            pl.plot(times, row['template'] + np.mean( row['raw_pulse'][0:int((PULSE_START-0.05)*fe)] )
                    ,label="Best Fit", color='red')
#            pl.plot(time, sglTemplateFunc(time, results["Best Fit Params"]), label="Best Fitting Template", color='r')
            pl.legend()
            pl.show()
            
#            pl.plot(times, row['hpf_pulse'], label="HPF Pulse", color='blue')
#            pl.xlabel("Time [s]")
#            pl.ylabel("Amplitude [ADU]")
#            pl.plot(times, high_pass_filter(row['template']) + np.mean( row['hpf_pulse'][0:int((PULSE_START-0.01)*fe)] )
#                    ,label="Best Fit", color='red')
##            pl.plot(time, sglTemplateFunc(time, results["Best Fit Params"]), label="Best Fitting Template", color='r')
#            pl.legend()
#            pl.show()
        
            print "chi sq for the fit was ", row['chi']
            print "compared to chi_zero of ", row['chi_zero']
            print "with amp of ", row['amplitude']
            print "with optimal time delay of ", row['delay']
            print "peak time", timestep_of_peak/fe

#        except:
#            print "\n\nI think timestep_of_peak is nan\n\n"
#            print timestep_of_peak
        
        counter+=1
        
    pl.hist(result_df[result_df.amplitude > 0]['amplitude'], bins=100)
    pl.xlabel('Amplitude')
    pl.show()
    
    pl.hist(result_df['delay'], bins=100)
    pl.xlabel('delay')
    pl.show()
          
    fig = pl.figure()
    ax = pl.gca()
    ax.plot(result_df['amplitude'] , result_df['chi'], '.', c='black', markeredgecolor='none')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Chi Sq.')
    ax.minorticks_on()
    pl.show()
    print result_df['chi'].describe()

    return result_df

'''  
===================================================  
Begin Main Code!
===================================================
'''
#file1 = 'data_run42_dbz1/20180314_16h12'; target = "Zinc"; #original neutron
#file1 = 'data_run42_dbz1/20180314_10h11'; target = "Zinc"; #Ba, fails after a couple hours
#file1 = 'data_run42_dbz1/20180313_18h32'; target = "Zinc"; #good Ba
file1 = 'data_run42_dbz1/20180315_14h24'; target = "Zinc"; #No source original#
#Ge_file = 'Germanium/se25g000_000'; Ge_offset = 7077; target = "Germanium";
#file1 = 'data_run42_dbz1/20180315_09h55' #SANS source post calibatration neutron# fails
#print "\n\nFile: ", file1, "\n"


fe = 1000. #use for Zinc data
#fe = 400. #Germanium data and turn off downsampling
if target=="Germanium": #checks for germanium
    assert fe==400
if target=="Zinc": #checks for zinc
    assert fe==1000
f_Nyq=fe/2.
PULSE_WINDOW = 2.0 #length in seconds of pulse windows
PULSE_START = 0.25 #time in seconds within pulse window that each pulse starts at
if target=="Germanium":
    assert PULSE_WINDOW == 1.0
    assert PULSE_START == 0.5
PRE_HP_FILTER = True
WINDOW = "boxcar"
#templates = ( (temp1, name1), ...  )

process_file(file1, output_folder = 'Results/' + 'cleancode_test3/', chunk_size=1, start_chunk=0, end_chunk=-1)
#process_file(Ge_file, output_folder = 'Results/' + 'Ge_Test/', chunk_size=1, start_chunk=0, end_chunk=1, offset=Ge_offset, downsampled=False, total_duration=24)
#process_file(file1, output_folder = 'Results/' + 'cleancode_test3/', chunk_size=0.01, start_chunk=0, end_chunk=1)
