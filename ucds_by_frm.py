"""
Unit-cell clustering

Jack Binns
"""

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import seaborn as sns
import math

sns.set()


def PeakHunt(profile, prom):
    peak_ind, _ = find_peaks(profile[..., 1], prominence=prom)
    peaks_in_q = []
    for i, datap in enumerate(profile):
        if i in peak_ind:
            peaks_in_q.append(profile[i][0])
    peaks_in_q = np.array(peaks_in_q)
    return peaks_in_q


def QMatrix(files, qmin, qmax):
    # generate matrix over the q range for all profiles
    # in the submitted data set
    qdata_all = []
    for fn in files:
        qdata = []
        data = np.loadtxt(fn)
        for point in data:
            if qmin <= point[0] <= qmax:
                qdata.append(point)
        qdata = np.array(qdata)
        qdata_all.append(qdata)
    np.array(qdata_all)
    return qdata_all


def Qtrim(profile, qmin, qmax):
    # trim any profile between q_min and q_max
    trim_len = 0
    for point in profile:
        if qmin <= point[0] <= qmax:
            trim_len = trim_len + 1
    trim = np.zeros([trim_len, 3], dtype=float)
    k = 0
    for j, point in enumerate(profile):
        if qmin <= point[0] <= qmax:
            trim[k] = profile[j]
            k = k + 1
    return trim


def DiffVec(test, model):
    # calculate a vector and return it's magnitude
    # If this is above some value we'll reject it
    d = 0.0
    if len(test) == 0:
        magd = 0.05
        return magd
    for i, peak in enumerate(model):
        d = d + (test[i] - peak) ** 2
        magd = np.sqrt(d)
        return magd


def RasterMap(log, datatag, outpath):
    xlen = 24
    xrem = len(log) % xlen
    zero_vec = np.zeros((xlen - xrem,))
    b = np.append(log, zero_vec)
    ylen = int(len(b) / xlen)
    c = b.reshape(ylen, xlen)
    np.savetxt(outpath + datatag + "_map.txt", c)


##
## Define the lattice ratios
##
pn3m = np.array([np.sqrt(2), np.sqrt(3), np.sqrt(4), np.sqrt(6), \
                 np.sqrt(8), np.sqrt(9), np.sqrt(10), np.sqrt(11), np.sqrt(12)])
hex = np.array([1, np.sqrt(3), np.sqrt(4)])
im3m = np.array([np.sqrt(2), np.sqrt(4), np.sqrt(6), np.sqrt(8), np.sqrt(10), np.sqrt(12), np.sqrt(14)])
ia3d = np.array(
    [np.sqrt(6), np.sqrt(8), np.sqrt(14), np.sqrt(16), np.sqrt(20), np.sqrt(22), np.sqrt(24), np.sqrt(26), np.sqrt(28)])
lam = np.array([1, 2])

if __name__ == "__main__":

    #############################################################################
    # CONTROL THROUGH THESE INPUTS

    datapath = "/home/jack/work/Martin_16186a/exp1/reduced/"
    # tifpath = "/media/jack/jack_drive/Martin_14680/Experiment1/images/"
    datatag = "Plate5_B1t_"
    outpath = "/home/jack/work/Martin_16186a/exp1/analysis/ucds_by_frm/"
    if os.path.isdir(outpath) == False:
        os.mkdir(outpath)

    # ' on Windows the file paths must include an extra \ : '

    # datapath = "C:\\path\\to\\the\\SAXS\\data\\"

    # don't forget the trailing \\

    # these values are the limits in q that we're looking at
    qmin = 0.05
    qmax = 0.3

    # this is a value used in the peak hunting algorithm
    prominence = 5

    #############################################################################
    ## Heading
    print("###############################################")
    print("Filtering data sets by peak fitting")
    print("###############################################")
    print("")
    print("")

    # Here we create a list 'files' which is a list of strings which
    # are the paths to each profile in the folder you've directed it to
    files = sorted(glob.glob(datapath + datatag + "*.dat"))
    print("Dataset : ", datatag)
    print("Total number of files in data set : ", len(files))

    # Parameters that define the 96 well latest
    x_scan_count = 6
    y_scan_count = 6
    x_step = .35
    y_step = 0.25
    x_big_step = 3.0
    x_big_count = 12

    # Controls the positions of the wells , top and bottom
    x = 0
    y_b = 0
    y_t = 1.8
    x_big = 0
    raw_list = []
    for i, dat in enumerate(files):
        print(dat)
        profile = np.loadtxt(dat, skiprows=2)
        profile = Qtrim(profile, qmin, qmax)
        peaks_in_q = PeakHunt(profile, prominence)
        sploot = dat.split("/")
        splat = sploot[-1].split("_")
        print(splat[-1][:-4])

        frm = int(splat[-1][:-4])
        if len(peaks_in_q) == 0:
            peaks_in_q = [0.0]
        print(frm, peaks_in_q[0])
        raw_list.append([int(frm), float(peaks_in_q[0])])

        # fig = plt.figure(figsize=(7, 5))
        # plt.plot(profile[..., 0], profile[..., 1])
        # for k in np.arange(len(peaks_in_q)):
        #     plt.axvline(x=peaks_in_q[k], ymin=0, ymax=0.2, color=sns.xkcd_rgb["pale red"])
        # plt.xlabel(r'q (nm$^{-1}$)')
        # plt.ylabel(r'Intensity (arbitrary units)')
        # plt.show()

    raw_list = sorted(raw_list)
    # print(raw_list)

    raw_list = np.array(raw_list)
    # output file:
    np.savetxt(outpath + datatag + 'frm_peak0.dat', raw_list)

    #
    #
    #
    #
    # # Display the average profile and the found peaks
    # fig = plt.figure(figsize=(7, 5))
    # plt.plot(average_profile[..., 0], average_profile[..., 1])
    # for i in np.arange(len(peaks_in_q)):
    #     plt.axvline(x=peaks_in_q[i], ymin=0, ymax=0.2, color=sns.xkcd_rgb["pale red"])
    # plt.xlabel(r'q (nm$^{-1}$)')
    # plt.ylabel(r'Intensity (arbitrary units)')
    # plt.show()
    #
    # # Here we manually set the index of which peak we want to use as the first peak for some phase
    # # in your script we will want to automate this
    # index = int(input("Which peak in the average profile should be used as the key? [py list rules]: "))
    # root_peak = float(peaks_in_q[index])
    # # 2% tolerance
    # tol = root_peak * 0.02
    # peak_num = len(peaks_in_q)
    #
    # # Which lattice type will be rooted to root_peak?
    # lattice = model
    #
    # # Find lattice parameter for chosen model
    # if model_name == "hex":
    #     latt_a = 2 / (math.sqrt(3) * root_peak)
    # else:
    #     latt_a = lattice[0] / root_peak
    # print("a = ", str("{0:.5f}".format(latt_a)), "nm")
    #
    # # Generate peaks for this model
    # print("Generating peaks for the model lattice...")
    # model_peaks = []
    # for h in range(len(lattice)):
    #     s_h = float(lattice[h]) / float(latt_a)
    #     model_peaks.append(s_h)
    #
    # model_peaks = np.array(model_peaks)

    # Don't think you'll need this
    # Create a matrix of all profiles in current working list within Q RoI
    # print("")
    # print("Generating a matrix of all ",str(len(files))," profiles...")
    # profile_matrix = QMatrix(files,qmin,qmax)

    # Let's play around with the difference measure
    # diff_list = []
    # for i, profile in enumerate(profile_matrix):
    #     peaks = PeakHunt(profile,prominence)
    #     diff = DiffVec(peaks,model_peaks)
    #     diff_list.append(diff)
    # plt.plot(diff_list)
    # plt.show()
    # np.savetxt(outpath+datatag+"_"+str(len(diff_list))+"_diff.txt",diff_list)
    # RasterMap(diff_list,datatag,outpath)
    #
    # Now to filter
    # print("")
    # print("Filtering profiles...")
    # rejects_f  = open(outpath+datatag+"_"+str(len(files))+"_rej.txt", "w+")
    # manifest_f  = open(outpath+datatag+"_"+str(len(files))+"_mf.txt", "w+")
    # manifest = []
    # rejects = []
    # diff_co = float(input("Difference cut-off :  "))
    # for i,profile in enumerate(profile_matrix):
    #     peaks = PeakHunt(profile,prominence)
    #     diff = DiffVec(peaks,model_peaks)
    #     if len(peaks) == 0:
    #         rejects.append(files[i])
    #         rejects_f.write(str(files[i])+"\n")
    #     elif diff < diff_co   :     #and len(peaks) == peak_num
    #         x = files[i].split("/")
    #         tiffname = tifpath + x[-1][:-4]+'.tif'
    #         manifest.append(files[i])
    #         manifest_f.write(str(tiffname)+"\n")
    #     else:
    #         print("REJECT")
    # rejects.append(files[i])
    # rejects_f.write(str(files[i])+"\n")

    # summarise the stats
    # print("")
    # print(len(rejects)," profiles were rejected ",'{0:.2f}'.format(len(rejects)/len(files)*100)," %" )
    # print(len(manifest)," profiles passed")

    # Generate a new average
    # print("")
    # print("Calculating updated average profile...")
    # ud_average_profile = AverageProfiles(manifest)
    # print("")
    # print("Updated profile written to: ",outpath)
    # np.savetxt(outpath+datatag+"_"+str(len(manifest))+"_average.txt",ud_average_profile)

    # Plot the updated average profile
    # print("Updated average plot: ")
    # ud_average_profile = Qtrim(ud_average_profile,qmin,qmax)
    # fig= plt.figure(figsize=(7,5))
    # plt.plot(ud_average_profile[...,0],ud_average_profile[...,1],color = "g")
    # plt.xlabel( r'q (nm$^{-1}$)' )
    # plt.ylabel( r'Intensity (arbitrary units)' )
    # plt.show()

    # Calculate the average of the reject profile
    # print("")
    # print("Calculating rejected profile...")
    # if len(rejects) > 0 :
    #     reject_profile     = AverageProfiles(rejects)
    #     reject_profile = Qtrim(reject_profile,qmin,qmax)
    # np.savetxt(outpath+datatag+"_"+str(len(manifest))+"_avdiff.txt",reject_profile)
    #
    # print("")
    # print("Reject profile:")
    # fig= plt.figure(figsize=(7,5))
    # plt.plot(reject_profile[...,0],reject_profile[...,1],color = "k")
    # plt.xlabel( r'q (nm$^{-1}$)' )
    # plt.ylabel( r'Intensity (arbitrary units)' )
    # plt.show()
