# -*- coding: utf-8 -*-
""" 
    Basic tools for stepfinding:
    -split_fast
    Jacob Kers 2021 
    """
import matplotlib.pyplot as plt
import numpy as np


def splitFast(segment, demo=0):
    """
    Determines the best step-fit in a one-dim array 'segment'.
    To save time, use of function 'mean' is avoided in the deep loop
    If the result is invalid, this is passed via 'rank' value=0
    Jacob Kers 2021
    """
    invalidFit = 0
    Nmin = 2
    Ns = len(segment)
    if Ns > 2:
        var_q = np.ones(
            Ns,
        )
        # set up averaging buffers
        avl = segment[0]
        avr = np.sum(segment[1 : Ns + 1]) / (Ns - 1)
        ava = np.sum(segment) / Ns
        for ii in range(1, Ns - 1, 1):
            # current plateau lengths:
            n_L = ii
            n_R = Ns - ii
            # update left average:
            avl = (avl * (n_L - 1) + segment[ii]) / (n_L)
            # update right average:
            avr = (avr * (n_R + 1) - segment[ii]) / (n_R)
            # relative change in averages left and right
            delta_l = avl - ava
            delta_r = avr - ava
            # total variance correction (NEUTRALIZED):
            varcor_L = 1 + 0 * (n_L + 1) / (n_L)
            varcor_R = 1 + 0 * (n_R + 1) / (n_R)
            # total variance:
            delta_q = delta_l**2 * n_L * varcor_L + delta_r**2 * (n_R) * varcor_R
            var_q[ii] = -delta_q
        # wrapping up:
        idx = np.argmin(var_q)
        if (idx < Nmin - 1) | (Ns - idx < Nmin - 1):
            invalidFit = 1
        else:
            avl_fin = np.mean(segment[0:idx])
            avr_fin = np.mean(segment[idx + 1 : Ns + 1])
            rankit = (avr_fin - avl_fin) ** 2 * Ns
            errorcurve = var_q / Ns
            # demo section:
            if demo - np.floor(demo) > 0.35:
                plt.close("All")
                fig, (ax1, ax2) = plt.subplots(2, 1)
                # data
                ax1.plot(segment)
                ax1.set_title("single step")
                ax1.set_xlabel("index (a.u)")
                ax1.set_ylabel("value (a.u)")
                # residual variance of stepfit
                ax2.plot(var_q, "ro-")
                ax2.set_title("residual variance of stepfit")
                ax2.set_xlabel("index (a.u)")
                ax2.set_ylabel("value (a.u)")
                fig.show()
                input("Press Enter to continue...")
                plt.close()
    else:
        invalidFit = 1
    if invalidFit:
        idx = 0
        avl_fin = segment[0]
        avr_fin = segment[0]
        rankit = 0
        errorcurve = segment * 0
    return idx, avl_fin, avr_fin, rankit, errorcurve


def Indices2Fit(dataX, indices, how="mean"):
    """ "This function builds a step fit"""
    ixlo = 0
    LX = len(dataX)
    FitX = 0 * dataX
    indices_ext = np.append(-1, indices)
    indices_ext = np.append(indices_ext, LX)
    Lix = len(indices_ext)
    for ii in range(0, Lix - 1, 1):
        ixlo = indices_ext[ii] + 1
        ixhi = indices_ext[ii + 1] + 1
        if ixhi >= ixlo:
            if how == "mean":
                FitX[ixlo:ixhi] = np.mean(dataX[ixlo:ixhi])
            elif how == "median":
                FitX[ixlo:ixhi] = np.median(dataX[ixlo:ixhi])
        else:
            FitX[ixlo:ixhi] = dataX[ixlo:ixhi]

    return FitX


def Fit2Steps(dataX, FitX):
    """
    Build a table of step properties from a step fit and the raw data
    """
    # get an index ('time') array:
    Lx = len(dataX)
    T = np.arange(Lx)
    # get a noise estimate via the residu:
    globalnoise = np.std(np.diff(dataX - FitX)) / 2**0.5
    # get indices of last points prior to steps and include start and end:
    ixes0 = np.ndarray.flatten(np.argwhere(np.diff(FitX) != 0))
    Lix = len(ixes0)
    ixes = np.hstack([[-1], ixes0, [Lx]])
    steptable = np.zeros(shape=(Lix, 8))
    # get data per step:
    for ii in range(1, Lix, 1):
        ix_pre = ixes[ii - 1]
        ix = ixes[ii]
        ix_aft = ixes[ii + 1]
        lev_pre = FitX[ix]
        lev_aft = FitX[ix + 1]
        step = FitX[ix + 1] - FitX[ix]
        dwell_pre = ix - ix_pre
        dwell_aft = ix_aft - ix
        # error based on global noise and local window,+/-95%:
        error_pred = (
            2
            * (globalnoise**2 / dwell_pre + globalnoise**2 / dwell_aft) ** 0.5
            / 2**0.5
        )
        # error based on local VAR,+/-95%:
        rms_pre = np.std(dataX[ix_pre + 1 : ix])
        rms_aft = np.std(dataX[ix + 1 : ix_aft + 1])
        error_meas = (
            2
            * ((rms_pre**2 / dwell_pre + rms_aft**2 / dwell_aft) ** 0.5)
            / 2**0.5
        )
        new_row_entry = np.array(
            [ix, lev_pre, lev_aft, step, dwell_pre, dwell_aft, error_pred, error_meas]
        )
        steptable[ii - 1] = new_row_entry

    return steptable


def hatCurve(N=100, edz=4):
    """
    Make a soft-edged 'hat' curve.
    To be used to suppress edge values of a profile
    """
    flatpart = np.ones(N - 2 * edz)
    edges = np.hanning(2 * edz)
    left = edges[0:edz]
    right = edges[edz : 2 * edz + 1]
    hat = np.concatenate((left, flatpart, right), axis=None)

    return hat


def AppendFitX(newFitX, FitX, dataX):
    """combine different fit rounds
    includes check for too closely spaced indices from differend rounds
    """
    # get step locations:
    combiFitX = FitX + newFitX
    Lx = len(combiFitX)
    ixes0 = np.ndarray.flatten(np.argwhere(np.diff(combiFitX) != 0))
    if len(ixes0) > 0:
        # pad with start and end
        Lix = len(ixes0)
        ixes = np.hstack([[-1], ixes0, [Lx]])
        # find indices too close together:
        Nmin = 2
        whereblips = np.argwhere(np.diff(ixes) < Nmin)
        # redo a stepfit over this part to pick just one step location:
        for ix in whereblips:
            lo = ix[0] - 1
            ixlo = ixes[lo]
            ixhi = ixes[lo + 3]
            segment = dataX[ixlo + 1 : ixhi]
            idx, avl, avr, rankit, errorcurve = splitFast(segment)
            combiFitX_old = np.copy(combiFitX)
            combiFitX[ixlo + 1 : ixlo + idx + 2] = avl
            combiFitX[ixlo + idx + 2 : ixhi + 1] = avr
    dum = 1
    return combiFitX
