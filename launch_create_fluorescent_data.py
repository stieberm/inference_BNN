#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 13:53:57 2018
Description: Launch the function create fluorescent data from spike times.
From the spike times data we can create the fluorescent time series.

Input:		   spikes.mat:		        matlab struct containing the spiking neurons
                                        and times
                                        
	           hyperparameters.mat	    matlab struct containing the hyperparameters
                                        used in the simulation
                                        
		       position.mat: 			matlab struct containing the position 
                                        of each neuron.
Output:        fluorescent_data.mat     matlab struct containing the fluorescent
                                        time series and the calcium concentrations
                                        

@author: Max Stieber
ETH Zürich
stieberm@ethz.ch
"""
from pathlib import Path
from create_fluorescent_data_methods import FluorescentParameters
from create_fluorescent_data_from_spike_times import create_fluorescent_data


# Define the hyperparameters
fluoro = FluorescentParameters()
fluoro.A_Ca         = 50                       # step amount for each action potential
fluoro.tau_Ca       = 0.01                     # decay constant
fluoro.noise        = 0.7                      # noise level for the network
fluoro.k_d          = 300                      # saturation concentration
fluoro.A_sc         = 0.15                     # constant defining the amount of influence of the neighboring neurons from fluorescent light
fluoro.overwrite_files  = 'yes' #'yes' or 'no'
fluoro.input_folder     = Path("./data/100x1_network/")
fluoro.output_folder    = Path("./data/100x1_network/")


fluoro.show_parameters

create_fluorescent_data(fluoro)