#!/bin/bash


# USER DEFINE VARIABLES: 

# First define the variables for slicing the mars data: 
raw_data_file='example_mars_data.MSEED'
slice_len=20000
overlap=50
output_name='mmse_ex'


# Now define the MSE stuff: 
r=0.1 
tau=30
scale=20

output_folder='./example/'
output_file='example_MMSE_output.txt'



# Download data and process using python function: 
python3 gen_MSE_slices.py $raw_data_file $slice_len $overlap $output_name




# Calculate Moving MSE for trace: 
matlab -nodisplay -r "moving_mse($r, $tau, $scale, '$output_name', '$output_folder', '$output_file');exit"
