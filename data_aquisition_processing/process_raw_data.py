import numpy as np
import os
import obspy.core.trace as TR

lst = os.listdir("../data/p0.2/")

out_dir = "../data/p0.2/processed/"

dec = 2
chls = "RTZ"

for folder in lst:
    if "p0.2" in folder:
        master_dir = f"../data/p0.2/{folder}/"
        print(f"master dir: {master_dir}")

        time = np.loadtxt(fname=f"{master_dir}/data_time.ascii")
        time = time[::dec]

        dt = time[1] - time[0]

        if os.path.isdir(f"{out_dir}/{folder}/")==False:
            os.mkdir(f"{out_dir}/{folder}/")



        time_out_str = f"{out_dir}/{folder}/time_data"
        np.savetxt(fname=time_out_str, X=np.transpose([time]))

        for file in os.listdir(master_dir):
            if "M_0" in file:
                print(f"Opening {file}")
                data = np.loadtxt(fname=f"{master_dir}/{file}")
                for ch in range(3):
                    tr = TR.Trace()
                    tr.stats.delta = dt
                    d = data[:,ch]
                    tr.data = d[::dec]
                    tr.filter(type="lowpass", freq=2.0)

                    out_str = f"{out_dir}/{folder}/{file[:-6]}.{chls[ch]}"
                    np.savetxt(fname=out_str, X=np.transpose([tr.data]))
                    print(f"Processed {out_str}")



