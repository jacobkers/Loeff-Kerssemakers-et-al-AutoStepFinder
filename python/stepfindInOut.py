import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
import glob


def getWorklist(demo=0):
    """
    get a directory name and one or more filenames to analyze
    """
    # 1) build file list: nr 1 is the single-file option (CHECK: bring to sio)
    if np.floor(demo) == 0:
        workPath = "nopath"
        workList = ["nofile.txt"]
    else:
        # get source (path, file) to 1 or more files
        root = tk.Tk()
        root.withdraw()
        fullpath = filedialog.askopenfilename()
        workPath = os.sep.join(fullpath.split("/")[:-1])
        curPath = os.getcwd()
        os.chdir(Path(workPath))
        workList = [f for f in glob.glob("*.txt")]
        os.chdir(curPath)
        if np.floor(demo) == 1:
            workList = [(fullpath.split("/")[-1])]  # fileList is just clicked one
    return workPath, workList


def LoadData(workPath, fileName):
    """Load single-column data"""
    loadpath = Path(workPath + "\\" + fileName)
    inoutName = loadpath.stem
    with open(loadpath, "r") as inf:
        dataX = np.loadtxt(inf)
    return dataX, inoutName


def SavePlot(workPath, inoutName, dataX, Fits, S_curves, best_shots, steptable, demo):
    """save and plot results"""
    # create header for csv file
    header = [
        "index",
        "level before",
        "level after",
        "step",
        "dwell before",
        "dwell after",
        "predicted error",
        "measured error",
    ]

    # create outpath
    if workPath != "nopath":
        outPath = workPath + str("\\steppy_results\\")
        if not Path(outPath).is_dir():
            Path(outPath).mkdir()
    else:
        outPath = ""

    with open(
        outPath + inoutName + str("_steps") + str(".csv"), "w"
    ) as csv_f:  # will overwrite existing
        # create the csv writer
        writer = csv.writer(csv_f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)

    # save characteristic variables to csv file in outpath
    with open(
        outPath + inoutName + str("_steps") + str(".csv"), "a"
    ) as csv_f:  # will overwrite existing  # will append to file
        # create the csv writer
        writer = csv.writer(csv_f, delimiter="\t", lineterminator="\n")
        # write data to the csv file
        sz = np.shape(steptable)
        for sti in range(sz[0]):
            step_row = steptable[sti]
            writer.writerow(step_row)
            # show final result
    if demo - np.floor(demo) >= 0.1:
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.plot(dataX)
        ax1.plot(np.transpose(Fits))
        ax1.set_title("intermediate & final result")
        ax1.set_xlabel("time, a.u.")
        ax1.set_ylabel("position")
        ax1.set_box_aspect(0.5)
        ax2.plot(0 * S_curves[0])
        ax2.plot(np.transpose(S_curves))
        if len(best_shots) == 1:
            ax2.plot(best_shots[0] - 1, S_curves[best_shots[0] - 1], "ro")
        else:
            for ii in range(len(best_shots)):
                ax2.plot(best_shots[ii] - 1, S_curves[ii, best_shots[ii] - 1], "ro")
        ax2.set_title("S-curve evaluation")
        ax2.set_xlabel("no. of steps in fit")
        ax2.set_ylabel("S-value, a.u.")
        ax2.set_box_aspect(0.7)
        plt.show(block=False)
        png_nm = outPath + inoutName + str("_fit") + str(".png")
        fig.savefig(png_nm, dpi=500)
        if demo < 1.5:  # if not in batch mode
            input("Press Enter to continue...")
        else:
            plt.pause(0.1)
        plt.close()


def SimulateData(N_st, demo=0):
    """make a simple test trace for testing purposes"""
    if N_st == 0:  # minimal, check small segments, no steps etc
        inoutName = "simMinimal"
        # segment = [-0.5, -0.5, 0.5]
        segment = np.random.randn(
            100,
        )
    elif N_st == 1:  # single step
        inoutName = "simsingleStep"
        hs = 10  # stepsize/2
        np.random.randn
        segment_left = np.random.randn(1, 100) - hs
        segment_right = np.random.randn(1, 100) + hs
        segment = np.concatenate((segment_left, segment_right), axis=None)
    elif N_st > 1:  # multi-step
        inoutName = "simstepTrain"
        # random or constant size
        stepstyle = "constant step"
        if stepstyle == "vary step":
            stepsizes = 1 + 4 * np.random.random_sample(
                N_st,
            )
        elif stepstyle == "constant step":
            stepsizes = 4 * np.ones(
                N_st,
            )
        stepsizes = np.rint(stepsizes)
        stepsizes = stepsizes.astype(int)

        # random or constant dwell
        dwellstyle = "constant dwell"
        if dwellstyle == "vary dwell":
            plateaulengths = 40 + 20 * np.random.random_sample(
                N_st,
            )
        elif dwellstyle == "constant dwell":
            plateaulengths = 100 + np.zeros(
                N_st,
            )
        plateaulengths = np.rint(plateaulengths)
        plateaulengths = plateaulengths.astype(int)

        # build the segments:
        segment = np.random.randn(1, plateaulengths[0])
        level = 0
        flipsign = 1
        for ii in range(0, N_st, 1):
            level = level + flipsign * stepsizes[ii]
            Nw = plateaulengths[ii]
            new_segment = (
                np.random.randn(
                    Nw,
                )
                + level
            )
            blockwave = build_blockwave(Nw)
            new_segment = new_segment + blockwave
            segment = np.concatenate((segment, new_segment), axis=None)
            # use to go down:
            flipsign = flipsign

    return segment, inoutName


def build_blockwave(N=200):
    """make an up-and down blockwave of length N"""
    sub_N = 10  # per plateau
    sub_hstp = 1
    reps = int(np.ceil(N / (2 * sub_N)))
    low_seg = (
        np.zeros(
            sub_N,
        )
        - sub_hstp
    )
    hi_seg = (
        np.zeros(
            sub_N,
        )
        + sub_hstp
    )
    segment = np.concatenate((low_seg, hi_seg), axis=None)
    blockwave = np.tile(segment, (reps))
    blockwave = blockwave[0:N]

    return blockwave
