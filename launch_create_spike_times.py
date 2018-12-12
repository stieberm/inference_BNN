#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 11:27:02 2018
Define all the hyperparameters and launch the create_spike_times script.
Description:	This python script allows to create the spike times, the membrane potential and the current.

Input :	    None

Output:    	spikes.mat :                matlab file containing a struct of the 
                                        spiking neurons and times
                                        
            spikes.split.mat :	        matlab file contains a struct of the spiking 
                                        neurons and times splitted into P different chunks.
                                        
			hyperparameters.mat: 	    matlab file contains struct holding the 
                                        hyperparameters used in the simulation
                                        
			position.mat: 			   	matlab struct containing the position 
                                        of each neuron.
                                        
			multimeter_data.mat: 	 	matlab struct containing the potential 
                                        and the current of each neuron
                                        
			connection.txt:			   	text file contains the connections of 
                                        the neural network in the following order:
        									source || target || weight || delay || position

@author: Max
"""
from pathlib import Path
from create_spike_times import create_spike_times
from create_spike_times_methods import SimulationParameters
from create_spike_times_methods import NeuralModelParameters


# Define hyperparameters for simulation
Sim = SimulationParameters()                         #Initialize the class containing the simulation relevant parameters
Sim.P                   = 1                                         # number of time series
Sim.resolution          = 0.01                                      # delta_t to solve the ODE's
Sim.interval            = 0.5                                       # sampling rate if the multi meter in ms
Sim.length_ts           = 2000                                      # number of samples of each time series
Sim.nr_samples          = Sim.length_ts * Sim.P                     # number of samples that need to be generated
Sim.sim_length          = Sim.nr_samples * Sim.interval             # simulation length in ms
Sim.data_folder         = Path("./data/100x1_network/")
Sim.overwrite_files     = "yes"                                      # "no" or "yes"

# Define hyperparameters for biological neural network.
NModel = NeuralModelParameters()
NModel.nr_neurons       = 100                           # 100 neurons or 16 neurons.
NModel.tau_m            = 20.0                          # membrane time constant 
NModel.tau_s            = 2.0                           # synaptic time constant
NModel.alpha_int        = 5.74                          # Constant, which regulates the connection strength between the neurons
NModel.alpha_ext        = 5.0                           # Constant, which regulates the influence of the Poisson spikes
NModel.U                = 0.3                           # fraction of the neurons in recovering state enter the active state
NModel.t_delay          = 0.5                           # time delay in ms
NModel.tau_rec          = 5000.0                         # Recovering time scale in ms
NModel.poisson_spike_rate = 1.6                         # rate at which Poisson spikes are emitted
NModel.V_thres          = 20.0                          # potential threshold

Sim.show_parameters()
NModel.show_parameters()

create_spike_times(Sim, NModel)
