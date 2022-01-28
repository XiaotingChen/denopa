# -*- coding: utf-8 -*-
# @Time    : 18-5-27 下午5:06
# @Author  : Matrix
# @Site    :
# @File    : call_peak.py
# @Software: PyCharm
from __future__ import print_function
import scipy as sp
from scipy import optimize
from scipy.stats import gamma
import pandas as pd
import multiprocessing as mp


def loss_fun(p, q):
    a, l = p
    x = [gamma.cdf(i, a, scale=l) for i in q]
    x = [(i - j)**2 for i, j in zip(x, [0.25, 0.5, 0.75])]
    return sum(x)


def get_para(din):
    q = [sp.percentile(din, i) for i in [25, 50, 75]]
    aInit, _, lInit = gamma.fit(din[din > 0], floc=0)
    y = optimize.minimize(loss_fun,
                          x0=sp.asarray([aInit, lInit]),
                          args=(q, ),
                          method="Nelder-Mead")
    assert y.success, "The estimation does not converge."
    # assert y.fun <= 0.01, "The last value %.2f is too large." % y.fun
    return y.x


def __call_for_a_chrom(para):
    """
    Core function for calling peaks (candidate of searching regions) from the coverage max min track. Parameters are
    unzipped from the para tuple as follows:
    :param chrname: The name of the chromosome called;
    :param data: the max_min track got by adATAC.signal_track_builder.make_max_min_track;
    :param p: the p value used to select the candidate summit points.
    :param pminth: the p-value threshold to determine whether the value of the min point is small enough
    :param merge_dist: when two mearby peaks whose distance less than this parameter, the two would be merged.
    :return: a data frame which mimic the bed format with following columns:
        0: the name of the chromosome;
        1: the start position of the peak;
        2: the stop position of the peak;
        3: the position of the summit (maximum) in the peak;
        4: the signal track value of the sumit point;
        5: the p-value of the peak (summit);
    """
    chrname, data, p, pminth, merge_dist = para
    pmax = get_para(data[data[:, 2] > 0, 1])
    pmin = get_para(data[data[:, 2] < 0, 1])
    maxth = gamma.ppf(1 - p, pmax[0], scale=pmax[1])
    minth = gamma.ppf(pminth, pmin[0], scale=pmin[1])
    major = sp.where((data[:, 2] > 0) & (data[:, 1] > maxth))[0]
    peaks = []
    for k in major:
        try:
            # if data[k, 0] - data[peaks[-1][-1], 0] > merge_dist:
            #     peaks.append([k])
            # else:
            #     peaks[-1].append(k)
            # The policy of merging mearby major maxs is changed to check if the vally between them are deep enough.
            if any([
                    data[i, 2] < 0 and data[i, 1] <= minth
                    for i in range(peaks[-1][-1], k)
            ]):
                peaks.append([k])
            else:
                peaks[-1].append(k)
        except IndexError:
            peaks.append([k])
    peaks.sort(key=lambda u: u[0])
    extend_peaks = []
    for k in peaks:
        ldx = k[0]
        while ldx >= 0:
            if data[ldx, 2] < 0 and data[ldx, 1] <= minth:
                break
            else:
                ldx -= 1
        ldx = max(ldx, 0)
        rdx = k[-1]
        while rdx < data.shape[0]:
            if data[rdx, 2] < 0 and data[rdx, 1] <= minth:
                break
            else:
                rdx += 1
        if rdx > data.shape[0] - 1:
            rdx = data.shape[0] - 1
        extend_peaks.append([ldx, rdx])
    if merge_dist == None:
        merge_dist = sp.mean(
            [data[b, 0] - data[a, 0] for a, b in extend_peaks]) * 4
    print("The merge distance is %f" % merge_dist)
    final_peaks = []
    extend_peaks.sort(key=lambda u: u[0])
    for pk in extend_peaks:
        try:
            if data[pk[0], 0] - data[final_peaks[-1][-1], 0] > merge_dist:
                final_peaks.append(pk)
            else:
                final_peaks[-1][-1] = pk[-1]
        except IndexError:
            final_peaks.append(pk)
    peak_pos = []
    for k1, k2 in final_peaks:
        s = sp.argmax(data[(k1):(k2), 1]) + k1
        m = data[s, 1]
        p = 1 - gamma.cdf(m, pmax[0], scale=pmax[1])
        peak_pos.append([chrname, data[k1, 0], data[k2, 0], data[s, 0], m, p])
    peak_pos = pd.DataFrame(peak_pos)
    peak_pos[1] = sp.asarray(peak_pos[1], dtype=int)
    peak_pos[2] = sp.asarray(peak_pos[2], dtype=int)
    # f = (peak_pos[2] - peak_pos[1]) <= sp.percentile(peak_pos[2] - peak_pos[1], 95)
    f = sp.ones(peak_pos.shape[0], dtype=bool)
    return peak_pos.loc[f, :].copy()


def call_candidate_regions(max_min_track,
                           p_value,
                           p_min=0.5,
                           merge_dist=None,
                           proc=1):
    params = [(k, v, p_value, p_min, merge_dist)
              for k, v in max_min_track.items()]
    pool = mp.Pool(proc)
    y = map(__call_for_a_chrom, params)
    pool.close()
    pool.join()
    y = pd.concat(y, axis=0, join="inner")
    y.index = range(y.shape[0])
    return y
