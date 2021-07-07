import sys

import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.signal.rotate import rotate2zne
from obspy.signal.filter import bandpass
import numpy as np

def download_mars(output_name, stime, etime, channel='BH?', fmin=0.125, fmax=2, samprate=1000, network='XB', station='ELYSE', location='02', client='IRIS',):

    # S0173a
    t1 = UTCDateTime(stime) # 2019-05-23T02:22:00.000
    t2 = UTCDateTime(etime) # "2019-05-23T02:45:00.000")

    # Fetch waveform from IRIS FDSN web service into a ObsPy stream object
    # and automatically attach correct response
    fdsn_client = Client(client)
    st = fdsn_client.get_waveforms(network=network, station=station, location=location,
                                   channel=channel, starttime=t1, endtime=t2,
                                   attach_response=True)
    # define a filter band to prevent amplifying noise during the deconvolution
    pre_filt = (0.005, 0.006, 30.0, 35.0)
    st.remove_response(output='DISP', pre_filt=pre_filt)

    # Integrate to get dislacement:
    for j in range(3):
        st[j].integrate()

    BHU_data = st[0].data
    BHV_data = st[1].data
    BHW_data = st[2].data

    # Rotate using channel metadata accessed from: https://www.seis-insight.eu/en/science/seis-data/seis-metadata-access
    output_tuple = rotate2zne(data_1=BHU_data, azimuth_1=135.1, dip_1=-29.4,
                              data_2=BHV_data, azimuth_2=15, dip_2=-29.2,
                              data_3=BHW_data, azimuth_3=255, dip_3=-29.7)


    # Assign BP-filtered, rotated data back to traces:
    for i_update_index in range(3):
        # Bandpass the rotated data and re-assign to the stream
        st[i_update_index].data = bandpass(data=output_tuple[i_update_index],
                                           freqmin=fmin,
                                           freqmax=fmax,
                                           df=st[i_update_index].stats.sampling_rate,
                                           zerophase=True)
        # Up sample at 1000 Hz in accordance with parameter space data
        st[i_update_index].interpolate(sampling_rate=samprate, method="linear")


    print('Downloaded...')
    print(st)
    # Use one of these to output the file

    np.savetxt(output_name + '.ascii', [st[0].data, st[1].data, st[2].data]) # 'int_S0173a_0.125_2BP.ascii'
    print('Saved to ', output_name)
    st.write(output_name + '.MSEED')
    print('Written as MSEED: ', output_name + '.MSEED')

if __name__ == "__main__":

    # Get args from input:
    output_name = sys.argv[1]
    stime = sys.argv[2]
    etime = sys.argv[3]

    download_mars(output_name, stime, etime, channel='BH?', fmin=0.125, fmax=2, samprate=1000, network='XB',
                      station='ELYSE', location='02', client='IRIS', )
