==README==
This folder contains the necessary files to
  1. build a biological neural network,
  2. simulate the spike time and
  3. create synthetic fluorescent data.
The folder also includes a Matlab file, which allows to explore the resulting time series.


==DATA==
fluorescence_iNet1_Size100_CC01inh.txt
This file contains fluorescent data from the Connectomics challenge (http://connectomics.chalearn.org). We want our data to be
similar to the time series that have been created for the challenge.

==Instructions==
1) Run script ‘launch_create_spike_times.py’ while specifying the parameters and the path called data_folder. (works only with the nest simulator)
2) Run script ‘launch_create_fluorescent_data_from_spike_times.py’  while making sure that the path called input_folder is equal to the data_folder specified in 1).
3) Run script ‘Compare_data.m’ to visualise the datasets you have created in the previous steps.

==Remark==
We have included some spike time files with their corresponding hyperparameters in case you want to rush to step 2) and avoid an installation
of the nest simulator.
