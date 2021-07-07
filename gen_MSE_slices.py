import numpy as np
from obspy import read
import sys

def gen_mse_slices(input_fname, slice_length, overlap, out_name=''):
    # Import trace/data:
    data = read(input_fname)
    print('Data loaded...')


    # Loop for each channel
    for ch in range(len(data)):

        # Calculating number of slices
        trace = data[ch].data
        data_len = len(trace)
        no_slices = int(np.floor((data_len-overlap)/(slice_length - overlap)))

        ## Initialise slices array:
        slices_array = np.zeros((no_slices, slice_length))

        # Normalise the trace:
        trace_norm = trace/np.amax(np.abs(trace))

        # Slice up the array:
        for sl in range(no_slices):
            slices_array[sl, :] = trace_norm[ (slice_length-overlap)*sl:((slice_length-overlap)*sl)+slice_length]


        # If output name is undefined then will be a modified version of input file:
        if out_name == '': 
            out_fname = f"./SLICES_{input_fname}_ch{ch}.txt"
        else:
            out_fname = f"./SLICES_{out_name}_ch{ch}.txt"

        # Output sliced data
        np.savetxt(out_fname, slices_array)
        print('Saved ', out_fname)

if __name__ == "__main__":

    # Get args from input:
    input_fname = sys.argv[1]
    slice_length = int(sys.argv[2])
    overlap = int(sys.argv[3])
    out_name = sys.argv[4]


    gen_mse_slices(input_fname, slice_length, overlap, out_name=out_name)
