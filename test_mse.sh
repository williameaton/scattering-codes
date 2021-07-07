#!/bin/bash


# USER DEFINE VARIABLES: 

# First define the variables for downloading the mars data: 
data_fname='example_mars_data'
stime='2019-05-23T02:22:00.000'
etime='2019-05-23T02:35:00.000'

# Now define the MSE stuff: 
r=0.1 
tau=30
scale=10
input_path='./example_mars_data.ascii'
output_folder='./example/'
output_file='example_MSE_output.txt'



# Download data and process using python function: 
python3 mars_download.py $data_fname $stime $etime


# Calculate MSE for trace: 
matlab -nodisplay -r "mars_MSE($r, $tau, $scale, '$data_fname', '$output_folder', '$output_file');exit"
