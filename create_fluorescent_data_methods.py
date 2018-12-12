#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 09:18:27 2018
The functions we use to create synthetic fluorescent data are written in this file
@author: Max
"""
import numpy as np
import math
from pathlib import Path
from scipy.linalg.blas import dsymm

class FluorescentParameters:
    def __init__(self, A_Ca = 50, tau_Ca = 0.01, noise = 0.4, k_d = 300, A_sc = 0.7,
                 overwrite_files = 'no', input_folder = "./data/", output_folder = "./data/"):  
        # Define the hyperparameters
        self.A_Ca       = A_Ca                       # step amount for each action potential
        self.tau_Ca     = tau_Ca                   # decay constant
        self.noise      = noise                    # noise level for the network
        self.k_d        = k_d                      # saturation concentration
        self.A_sc       = A_sc                    # constant defining the amount of influence of the neighboring neurons from fluorescent light
        self.overwrite_files    = 'no' #'yes' or 'no'
        self.input_folder       = Path(input_folder)
        self.output_folder      = Path(output_folder)
        
    def show_parameters(self):
        print( '-------------------Fluorescent model parameters-------------------' )
        print( 'A_Ca                    : ' + str(self.A_Ca) )
        print( 'tau_Ca                  : ' + str(self.tau_Ca) )
        print( 'noise                   : ' + str(self.noise) ) 
        print( 'k_d                     : ' + str(self.k_d) )
        print( 'A_sc                    : ' + str(self.A_sc) )
        print( 'overwrite_files         : ' + self.overwrite_files )
        print( 'input folder            : ' + str(self.input_folder) )
        print( 'output folder           : ' +str(self.output_folder) )
        print( '-----------------------------------------------------------------' )

    

# Calculate the distance between the neurons and write them in a matrix
def get_distance(pos,size):
    
    d = np.zeros((size,size))
    
    for i in range(0, size):
        for j in range(0,size):
            d[i,j] = math.sqrt( (pos[i][0] - pos[j][0])**2 + ( pos[i][1] - pos[j][1] )**2 )
            
    return d


def make_n(spikes, interval, nr_samples, nr_neurons):
    
    spike_time_round = np.round(spikes[...,1]/interval)
    spike_time_round = spike_time_round.astype(int)
    n = np.zeros( (nr_samples,nr_neurons) )
    
    for i in range(0,len(spike_time_round)):
        n[spike_time_round[i] - 1,int(spikes[i,0])-1] = n[spike_time_round[i] - 1,int(spikes[i,0])-1] + 1.0
        
    return n


def Ca_model(n,tau_Ca,A_Ca):
    
    Ca = np.zeros(n.shape)
    
    for i in range(0,n.shape[1]):
        for j in range(0,n.shape[0]-1):
            Ca[j+1,i] = -Ca[j,i]*tau_Ca + A_Ca*n[j+1,i] + Ca[j,i]
            
    return Ca


def F_model(Ca,k_d,noise):
    
    F=np.zeros(Ca.shape)
    
    for i in range(0,Ca.shape[1]):
        for j in range(0,Ca.shape[0]):
            F[j,i] = Ca[j,i]/(Ca[j,i] + k_d) + noise*np.random.normal(0,0.03)
            
    return F
        
        
def F_sc_model(F,d,A_sc):
    
    F_sc = np.zeros(F.shape)
    sum_fluoro = np.zeros(F.shape)
    
    for i in range(0,F.shape[1]):
        for k in range(0,F.shape[0]):
            for j in range(0,F.shape[1]):
                if i!=j:
                    sum_fluoro[k,i] += F[k,j]*np.exp(-(d[i,j]/0.15)**2)
            F_sc[k,i] = F[k,i] + A_sc*sum_fluoro[k,i]
            
    return F_sc
    

def F_sc_model_opti(F, d, A_sc):
    
    F_sc = np.zeros(F.shape)
    C = np.zeros( (F.shape[1], F.shape[1]) )
    
    for i in range(0, F.shape[1]):
        for j in range(0, i):
            C[i,j] = np.exp(-(d[i,j]/0.15)**2)
            C[j,i] = C[i,j]
            if C[i,j] < math.pow(10,-5):
                C[i,j] = 0

    S = np.dot(C, F.T)
    
    F_sc = F + A_sc * S.T
    
    return F_sc


def F_sc_model_blas(F, d, A_sc):
    
    F_sc = np.zeros(F.shape)
    C = np.zeros( (F.shape[1], F.shape[1]) )
    
    for i in range(0, F.shape[1]):
        for j in range(0, i):
            C[i,j] = np.exp(-(d[i,j]/0.15)**2)
#            if C[i,j] < math.pow(10,-5):
#                C[i,j] = 0
                
    S = dsymm(1.0, C, F.T,lower=True)
    F_sc = F + A_sc * S.T
    
    return F_sc
                

def time_series_split(data,N):
    
    if data.shape[0]%N == 0:
       data_split = np.dstack(np.split(data,N))
       return [data_split]
    else:
        print("time series can not be divided in N equisized chuncks.")


    