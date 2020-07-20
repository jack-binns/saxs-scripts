import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import seaborn as sns
import math

sns.set()


def last_4chars(x):
    sploot = x.split("/")
    # print(sploot[-1][-8:-4])
    return sploot[-1][-8:-4]


def averageprofiles(file_list):
    # First find the length of the data
    look = np.loadtxt(file_list[0], skiprows=2)
    print(look.shape)
    # Average a profile over the full length
    print("")
    print("Averaging radial profile for " + str(len(files)) + " files...")
    datasum = np.zeros([look.shape[0]], dtype=float)  ## create empty sum
    c = 0
    for i, file in enumerate(file_list):
        data = np.loadtxt(file, skiprows=2)
        datasum = datasum + data[:, 1]
        c += 1
    datasum *= 1.0 / float(c)
    av_pro = np.array([data[:, 0], datasum])
    av_pro = np.transpose(av_pro)
    print("...Updated average radial profile")
    return av_pro


def rfac_calc(profile, model):
    # print(profile.shape)
    # print(model.shape)
    r_it = 0.0
    r_av = 0.0
    for k, pnt in enumerate(profile):
        # print(pnt[1], model[k][1])
        r_it = r_it + (pnt[1] - model[k][1]) ** 2
        r_av = r_av + model[k][1] ** 2
        r_tot = r_it / r_av
    return r_tot


if __name__ == "__main__":
    datapath = "/home/jack/work/Martin_16186a/exp1/reduced/"
    # tifpath = "/media/jack/jack_drive/Martin_14680/Experiment1/images/"
    datatag = "Plate3_A*"
    outpath = "/home/jack/work/Martin_16186a/exp1/analysis/rfac_nm1/"
    # outpath = "/home/jack/python/saxs-scripts/test/"

    datasets = ["Plate4_A*", "Plate5_B*", "Plate3_A*", "Plate2_A*", "Plate5_A*", "Plate1_B*", "plate_al2_30_r1_", "plate_al2_25_r1_", "testplate_al_5_r1_", "testplate_al_4_r1_", "testplate_al_3_r1_", "plate_al2_26_r1_", "plate_al2_27_r1_"]

    for datatag in datasets:
        if not os.path.isdir(outpath):
            os.mkdir(outpath)

        # these values are the limits in q that we're looking at
        qmin = 0.05
        qmax = 0.3

        # this is a value used in the peak hunting algorithm
        prominence = 5

        prefiles = glob.glob(datapath + datatag + "*.dat")
        files = sorted(prefiles, key=last_4chars)
        print(files)
        print("Dataset : ", datatag)
        print("Total number of files in data set : ", len(files))

        # Let's sort the

        # print(average)

        # Now let's compute the R-factor vs average
        rfac_vs_dset_av = []
        average = averageprofiles(files)
        for i, fn in enumerate(files):
            print(fn)
            prof = np.loadtxt(fn, skiprows=2, usecols=(0,1))
            # print(prof[0])
            r = rfac_calc(prof, average)
            rfac_vs_dset_av.append([i, r])
        np.savetxt(outpath + datatag + '_rfac_vs_dset_av.dat', rfac_vs_dset_av)

        # # Rolling R-factor vs previous frm
        # rfac_vs_nm1 = []
        # for i, fn in enumerate(files[1:]):
        #     # print(fn, files[i - 1])
        #     prof = np.loadtxt(fn, skiprows=2, usecols=(0, 1))
        #     prof_nm1 = np.loadtxt(files[i - 1], skiprows=2, usecols=(0, 1))
        #     r = rfac_calc(prof, prof_nm1)
        #     rfac_vs_nm1.append([i, r])
        # np.savetxt(outpath + datatag + '_rfac_vs_nm1.dat', rfac_vs_nm1)
