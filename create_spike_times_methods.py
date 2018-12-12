#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 15:19:33 2018

@author: Max
"""
from pathlib import Path
import numpy as np
import nest

class SimulationParameters:
    def __init__(self, P = 1, resolution = 1.0, interval = 1.0, length_ts = 2000, data_folder = "./data/"):
        self.P = P                                          # number of time series
        self.resolution = resolution                        # delta_t to solve the ODE's
        self.interval = interval                            # sampling rate if the multi meter in ms
        self.length_ts = length_ts                          # number of samples of each time series
        self.nr_samples =self.length_ts * self.P            # number of samples that need to be generated
        self.sim_length = self.nr_samples * self.interval   # simulation length in ms
        self.data_folder = Path(data_folder)                # Defime the folder where the data is stored
        self.overwrite_files = "no"                         # "no" or "yes" Default is no.
        
    def show_parameters(self):
        print( '----------------------Simulation parameters----------------------' )
        print( 'Number of time series                           : ' + str(self.P) )
        print( 'Length of each time series                      : ' + str(self.length_ts) )
        print( 'Resolution of the solver                        : ' + str(self.resolution))
        print( 'Sampling interval of the multimeter data        : ' + str(self.interval) )
        print( 'Overwrite the files                             : ' + self.overwrite_files )
        print( 'Save the data in the following folder           : ' + str(self.data_folder) )
        print( '-----------------------------------------------------------------' )
        
    
class NeuralModelParameters:
    def __init__(self, nr_neurons = 16, tau_m = 20.0 , tau_s = 2.0, alpha_int = 6.0, 
                 alpha_ext = 6.0, U = 0.3, t_delay = 0.5, tau_rec = 500.0, 
                 poisson_spike_rate = 1.6, V_thres = 20.0):  
        self.nr_neurons = nr_neurons                            # 100 neurons or 16 neurons.
        self.tau_m      = tau_m                                 # membrane time constant 
        self.tau_s      = tau_s                                 # synaptic time constant
        self.alpha_int  = alpha_int                             # Constant, which regulates the connection strength between the neurons
        self.alpha_ext  = alpha_ext                             # Constant, which regulates the influence of the Poisson spikes
        self.U          = U                                     # fraction of the neurons in recovering state enter the active state
        self.t_delay    = t_delay                               # time delay in ms
        self.tau_rec    = tau_rec                               # Recovering time scale in ms
        self.poisson_spike_rate = poisson_spike_rate            # rate at which Poisson spikes are emitted
        self.V_thres    = V_thres                               # value at which the membrane potential is thresholded
        
        
    def show_parameters(self):
        print( '---------------------Neural network parameters---------------------' )
        print( 'Number of neurons                   : ' + str(self.nr_neurons) )
        print( 'tau_m                               : ' + str(self.tau_m) )
        print( 'tau_s                               : ' + str(self.tau_s) ) 
        print( 'alpha_int                           : ' + str(self.alpha_int) )
        print( 'alpha_ext                           : ' + str(self.alpha_ext) )
        print( 'U                                   : ' + str(self.U) )
        print( 't_delay                             : ' + str(self.t_delay) )
        print( 'tau_rec                             : ' + str(self.tau_rec) )
        print( 'poisson spike rate                  : ' + str(self.poisson_spike_rate) )
        print( 'V_thres                             : ' +str(self.V_thres) )
        print( '----------------------------------------------------------------' )
        
        
    
def get_data(multimeter,spike_detector,nr_neurons,nr_samples):
        potential =     np.zeros((nr_samples, nr_neurons))
        currents_ex =   np.zeros((nr_samples, nr_neurons))
        potential =     np.zeros((nr_samples,nr_neurons))
        for i in range(0,nr_neurons):
            potential[...,i] = nest.GetStatus(multimeter, "events")[i]["V_m"]
            currents_ex[...,i]=nest.GetStatus(multimeter, "events")[i]["I_syn_ex"]
        times = nest.GetStatus(multimeter, "events")[0]["times"]
        spike_event = nest.GetStatus(spike_detector,"events")
        spike_time = spike_event[0]['times']
        spike_sender = spike_event[0]['senders']-1
        spike_time = spike_time[...,np.newaxis]
        spikes = spike_sender
        spikes = spikes[...,np.newaxis]
        spikes = np.append(spikes,spike_time,axis=1)
        return [potential, currents_ex, spikes, times]
    
    
    
def time_series_split(data,N):
    
    if data.shape[0]%N == 0:
       data_split = np.dstack(np.split(data,N))
       return [data_split]
    else:
        print("time series can not be divided in N equisized chuncks.")
        
        

def spike_times_split(spikes, sim_length, P):
    
    categorie = np.linspace(0, sim_length, num = P+1)
    digits = np.digitize(spikes[...,1], bins = categorie)
    
    for i in range(1,P+1):
        if i==1:
            string="N%o"%(i)
            spikesdict={}
            spikesdict={string : spikes[digits == i,...]}    
        else:
            string="N%d"%(i)
            spikesdict[string] = spikes[digits == i,...]
            spikesdict[string][...,1] = spikesdict[string][...,1] - categorie[i-1] 
            
        
            
    return spikesdict
