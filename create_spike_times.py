#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 17:49:10 2018
This script allows to create a neural network, save the spike times and the 
data recorded by a "digital" multimeter. 
@author: Max
"""
#into console: %reset -f
print(__doc__)
from pathlib import Path
import sys
import nest
import numpy as np
import nest.topology as topp
from create_spike_times_methods import get_data
from create_spike_times_methods import time_series_split
from create_spike_times_methods import spike_times_split
#from create_fluorescent_data_functions import SimulationParameters
#from create_fluorescent_data_functions import NeuralModelParameters
import scipy.io as sio

def create_spike_times(Sim, NModel):
    
    if Sim.overwrite_files == "no":
        if Sim.data_folder.is_dir()==True:
            print("You are not allowed to overwrite the files.")
            sys.exit("Error message")
        else:
            Sim.data_folder.mkdir(parents=True)
    elif Sim.overwrite_files == "yes":
        if Sim.data_folder.is_dir()==True:
            print("You have overwritten old files.")
        else:
            Sim.data_folder.mkdir(parents=True)
    else:
        print("You have two options for overwrite files")
        
    connections_path        = Sim.data_folder / "connections.txt"
    spikes_path             = Sim.data_folder / "spikes.mat"
    spikes_split_path       = Sim.data_folder / "spikes_split.mat"
    multimeter_data_path    = Sim.data_folder / "multimeter_data.mat"
    hyperparameters_path    = Sim.data_folder / "hyperparameters.mat"
    position_path           = Sim.data_folder / "position.mat"
    
    # Define important simulation parameters
    nest.ResetKernel()
    seed=1008.0
    nest.SetKernelStatus({"resolution": Sim.resolution,
                          "print_time": True,
                          "overwrite_files":True,
                          "grng_seed": int(seed),
                          "rng_seeds": [int(seed)]
                          })
    
    # Construct the position grid of the neural network(NN)
    jit = 0.03
    if NModel.nr_neurons==100:
        xs = np.arange(-0.45,.451,0.1) # defines the size of the network
    elif NModel.nr_neurons==16:
        xs = np.arange(-0.15,.151,0.1)
    else:
        print("Current network can constitute of 16 or 100 neurons.")
    np.random.seed(int(seed))
    pos = [[x,y] for y in xs for x in xs]
    pos = [[p[0]+np.random.uniform(-jit,jit),p[1]+np.random.uniform(-jit,jit)] for p in pos]
    
    # Construct the neurons on the grid and establish connections between them.
    # The probabilty of connection varies with the distance between the neurons
    # Define synapse connections
    nest.SetDefaults("tsodyks_synapse",{"delay": NModel.t_delay, #1.5 in Stetter's code
                                        "tau_rec": NModel.tau_rec,
                                        "tau_fac":0.0,
                                        "U": NModel.U
                                        })
    conn1 = {  "connection_type":"divergent",
                "mask": {"circular":{"radius":0.75}},
                "kernel": {"gaussian":{"p_center":1.,"sigma":0.15}}, #0.15 for 100 neurons
                "allow_autapses":False,
                "synapse_model":"tsodyks_synapse",
                "weights": NModel.alpha_int
                }
    
    # specify the neural model
    neuron_param=  {
                    #"I_e"       : 0.0,
                    "C_m"       : 1.0,
                    "tau_m"     : NModel.tau_m,
                    "t_ref"     : NModel.tau_s, #refactory periods in ms 2.0 is default
                    "E_L"       : 0.0,
                    "V_th"      : NModel.V_thres,
                    "V_m"       : 0.0,
                    "V_reset"   : 0.0
                    }
    nest.SetDefaults("iaf_psc_alpha", neuron_param)
    
    layer_dict_ex = {"positions": pos,
            "extent" : [1.1,1.1],
            "elements" : "iaf_psc_alpha"}
    layer = topp.CreateLayer(layer_dict_ex)
    
    topp.ConnectLayers(layer,layer,conn1)
    
    # Plot layer
    topp.PlotLayer(layer)
    
    # change the seed for different Poisson spike trains
    nest.SetKernelStatus({
            'grng_seed': int(seed),
            'rng_seeds': [int(seed)]
            })
    
    # Creation of a poisson generator
    nest.CopyModel('poisson_generator', 'PG',
                   params={'rate': NModel.poisson_spike_rate}) #1.6 in the paper, I don't know why they changed it in the programm
    pg = topp.CreateLayer({ 'rows' : 1,
                           'columns' : 1,
                          'elements' : 'PG'})
    cdict_stim = {'connection_type' : 'divergent',
                  'weights': NModel.alpha_ext}
    topp.ConnectLayers(pg,layer,cdict_stim)
    
    # create multimeter
    nrns=nest.GetLeaves(layer,local_only=True)[0]
    multimeter = nest.Create("multimeter", NModel.nr_neurons)
    nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m","I_syn_ex"],"interval":Sim.interval}) #, "input_currents_ex","input_currents_in"
    nest.Connect(multimeter,nrns,"one_to_one")
    
    # Create spike detector
    sd1 = nest.Create('spike_detector')
    nest.SetStatus(sd1,{'precise_times':True})
    nest.Connect(nrns,sd1)
    
    # Simulate
    nest.Simulate(Sim.sim_length + Sim.interval)
    
    # Retrieve the generated data
    [potential,currents_ex,spikes,time]=get_data(multimeter,sd1,NModel.nr_neurons,Sim.nr_samples)
    
    # Save the spikes to a file
    spikesdict={'N1':spikes}
    sio.savemat(str(spikes_path),spikesdict)
    
    if Sim.P>1:
        potential_fin = time_series_split(potential,Sim.P)[0]
        currents_ex_fin = time_series_split(currents_ex,Sim.P)[0]
        spikesdict_split = spike_times_split(spikes, Sim.sim_length, Sim.P)
    else:
        potential_fin=potential
        currents_ex_fin=currents_ex
        spikesdict_split = {'N1':spikes}
        
    # Pass important hyperparameters to create_fluorescent_data_from_spike_times.py
    hyperdict = {
            'P' : Sim.P,
            'length_ts' : Sim.length_ts,
            'nr_samples' : Sim.nr_samples,
            'nr_neurons' : NModel.nr_neurons,
            'interval' : Sim.interval
            
            }
    sio.savemat(str(hyperparameters_path),hyperdict)
    # save the position
    sio.savemat(str(position_path), {'position' : pos})

    # save the results
    sio.savemat(str(multimeter_data_path),mdict={ 'potential' : potential_fin,
                               'input_currents_ex': currents_ex_fin
                               })
    sio.savemat(str(spikes_split_path),spikesdict_split)
    
    topp.DumpLayerConnections(layer,'tsodyks_synapse',str(connections_path))