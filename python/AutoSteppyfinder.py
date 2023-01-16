# -*- coding: utf-8 -*-
## Credits
"""
Created 2022
@author: Jacob Kerssemakers
Title: 'Steppyfinder'
Description: simplified python transcript of matlabcode from [1] to detect steps in data traces

## settings
#inputs: 
single-column DataX

## run: multiPass(demo, tresH, N_iter) 

## run options
* demo number sets the different modii: 
# first digit
0 create a simple trace; 
1 (default): load single trace
2 run batch style on a directory

 2nd digit: 
.0: only save; 
.1 (default) show and save end graph 
.2 show intermediate passes
.3 show intermediate trace fit
.4 show intermediate fit segments
note: .2-.4 are suppressed in batch mode

* S_peak_treshold (0.15)
* iteration_range (to run at max=tracelenth/4, set to 0)

## example: 
   multiPass(demo=2.1, tresH=0.15, N_iter=100) 
    will run:
     -batch-style fit on a user-chosen directory with txt files, 
     ... with a fit round acceptance of 0.15
     ... with max. 100 iterations per fit
     ... and save both the step table and a .png fit picture to a 'result' directory 

## Refs:
[1] AutoStepfinder: A fast and automated step detection method for single-molecule analysis
Luuk Loef
f, Jacob W J Kerssemakers, Chirlmin Joo, Cees Dekker
Patterns . 2021 Apr 30;2(5):100256. 
doi: 10.1016/j.patter.2021.100256.
Matlab package on: https://zenodo.org/record/4657659#.Y6LxUXbMI2y

[2] Assembly dynamics of microtubules at molecular resolution
Jacob W J Kerssemakers et al
Nature 2006 Aug 10;442(7103):709-12. 
doi: 10.1038/nature04928.

"""
################################################################################
## Libraries
import matplotlib.pyplot as plt
import numpy as np

# custom
import stepfindCore as core
import stepfindInOut as sio
import stepfindTools as st

"""This section contains the 'core' loop of the stepfinder:
a single full iteration is done and a best fit is determined
@author: jkerssemakers march 2022       
"""


def multiPass(demo=1.0, tresH=0.15, N_iter=100):
    """This is the main, multi-pass loop of the autostepfinder
    @author: jkerssemakers march 2022"""
    # get worklist
    workPath, workList = sio.getWorklist(demo)

    # 2 run through worklist:
    fL = len(workList)
    S_curves = []
    for fi in range(0, fL, 1):
        if workList[0] == "nofile.txt":
            dataX, inoutName = sio.SimulateData(N_st=20)
        else:
            fileName = workList[fi]
            dataX, inoutName = sio.LoadData(workPath, fileName)
        FitX = 0 * dataX
        print("working on:" + inoutName)
        # multipass:
        for ii in range(0, 3, 1):
            # work remaining part of data:
            residuX = dataX - FitX
            newFitX, _, _, S_curve, best_shot = core.stepfindcore(
                residuX, demo, tresH, N_iter
            )
            FitX = st.AppendFitX(newFitX, FitX, dataX)
            # storage for plotting:
            if ii == 0:
                Fits = np.copy(FitX)
                S_curves = np.copy(S_curve)
                best_shots = [best_shot]
            elif best_shot > 0:
                Fits = np.vstack([Fits, FitX])
                S_curves = np.vstack([S_curves, S_curve])
                best_shots = np.hstack([best_shots, best_shot])

        # steps from final fit:
        steptable = st.Fit2Steps(dataX, FitX)
        sio.SavePlot(
            workPath, inoutName, dataX, Fits, S_curves, best_shots, steptable, demo
        )


# run -----------------------------------------
multiPass(demo=0.1, tresH=0.15, N_iter=1000)
""" # run options
* demo number sets the different modii: 
# first digit
0 create a simple trace; 
1 (default): load single trace
2 run batch style on a directory

 2nd digit: 
.0: only save; 
.1 (default) show and save end graph 
.2 show intermediate passes
.3 show intermediate trace fit
.4 show intermediate fit segments
note: .2-.4 are suppressed in batch mode

* S_peak_treshold (0.15)
* iteration_range (to run at max=tracelenth/4, set to 0) """
