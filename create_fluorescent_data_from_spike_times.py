#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 09:16:18 2018
This file takes the spike times and some parameters and outputs fluroescent times series
@author: Max
"""

#into console: %reset -f
print(__doc__)
#from sklearn.svm import LinearSVC
#from scipy.special import erf
#import pylab
from pathlib import Path
import scipy.io as sio
import time
import sys
from create_fluorescent_data_methods import get_distance
from create_fluorescent_data_methods import make_n
from create_fluorescent_data_methods import Ca_model
from create_fluorescent_data_methods import F_model
#from create_fluorescent_data_methods import F_sc_model_opti
from create_fluorescent_data_methods import time_series_split
from create_fluorescent_data_methods import F_sc_model_blas
import warnings
warnings.filterwarnings('ignore')



def create_fluorescent_data(fluoro_params):
    tic_total = time.clock()
    
    # Check overwrite file
    if fluoro_params.overwrite_files == "no":
        if fluoro_params.output_folder.is_dir() == True:
            if fluoro_params.input_folder != fluoro_params.output_folder:
                print("You are not allowed to overwrite the files.")
                sys.exit("Error message")
        else:
            fluoro_params.output_folder.mkdir(parents=True)
    elif fluoro_params.overwrite_files == "yes":
        if fluoro_params.output_folder.is_dir()==True:
            print("You have overwritten old files.")
        else:
            fluoro_params.output_folder.mkdir(parents=True)
    else:
        print("You have two options for overwrite files")
        
    
    # load data paths
    spikes_path             = fluoro_params.input_folder / "spikes.mat"
    hyperparameters_path    = fluoro_params.input_folder / "hyperparameters.mat"
    position_path           = fluoro_params.input_folder / "position.mat"
    
    # save data paths
    fluorescent_data = fluoro_params.output_folder / "fluorescent_data.mat"
    
    # load data
    spikes          = sio.loadmat(str(spikes_path))['N1']
    hyperparameters = sio.loadmat(str(hyperparameters_path))
    position        = sio.loadmat(str(position_path))
    P               = hyperparameters['P'][0][0]
    length_ts       = hyperparameters['length_ts'][0][0]
    nr_neurons      = hyperparameters['nr_neurons'][0][0]
    nr_samples      = hyperparameters['nr_samples'][0][0]
    delta_t         = hyperparameters['interval'][0][0]
    pos             = position['position']
    
    # Built a network with the nest simulator
    n   = make_n( spikes,delta_t,nr_samples,nr_neurons )
    
    Ca  = Ca_model( n,fluoro_params.tau_Ca, fluoro_params.A_Ca )
    
    F   = F_model( Ca, fluoro_params.k_d, fluoro_params.noise )
    
    d   = get_distance(pos,nr_neurons)
    
#    tic_function_F_sc_opti = time.clock()
#    F_sc = F_sc_model_opti(F, d, fluoro_params.A_sc)
#    toc_function_F_sc_opti = time.clock()
    
    tic_function_F_sc_blas = time.clock()
    F_sc = F_sc_model_blas(F, d, fluoro_params.A_sc)
    toc_function_F_sc_blas = time.clock()
    
#    timer_function_F_sc_opti = toc_function_F_sc_opti - tic_function_F_sc_opti
    timer_function_F_sc_blas = toc_function_F_sc_blas - tic_function_F_sc_blas
    
    # split the data
    if P>1:
        Ca_fin = time_series_split(Ca,P)[0]
        F_sc_fin = time_series_split(F_sc,P)[0]
    else:
        Ca_fin = Ca
        F_sc_fin = F_sc
    
    # save the results  
    sio.savemat(str(fluorescent_data),mdict = {
                               'Ca' : Ca_fin,
                               'F_sc' : F_sc_fin
                               })
        
    toc_total = time.clock()
    timer_total = toc_total - tic_total
    print( 'Total simulation time                           = ' + str(timer_total) + '  seconds')
#    print( 'Numpy matrix multiplication simulation time     = ' + str(timer_function_F_sc_opti) + '  seconds')
    print( 'BLAS matrix multiplication simulation time      = ' + str(timer_function_F_sc_blas) + '  seconds')
    