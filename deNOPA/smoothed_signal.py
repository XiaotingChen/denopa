# -*- coding: utf-8 -*-
# @Time    : 18-5-24 下午10:14
# @Author  : Matrix
# @Site    : 
# @File    : smoothed_signal.py
# @Software: PyCharm

import h5py
import numpy as np
from scipy.stats import norm
import gc
import itertools as it


def get_gauss_kernel0(scale):
    scale = float(scale)
    h = np.arange(-3 * int(scale), 3 * int(scale) + 1, dtype=np.float64)
    y = norm.pdf(h, scale=scale)
    y[0] /= 2
    y[-1] /= 2
    return y


def get_gauss_kernel1(scale):
    scale = float(scale)
    h = np.arange(-3 * int(scale), 3 * int(scale) + 1, dtype=np.float64)
    y = -1 / (scale ** 2) * h * norm.pdf(h, scale=scale)
    y[0] = y[0] / 2
    y[-1] = y[-1] / 2
    return y


def get_gauss_kernel2(scale):
    scale = float(scale)
    h = np.arange(-3 * int(scale), 3 * int(scale) + 1, dtype=np.float64)
    y = 1 / (scale ** 4) * (h ** 2 - scale ** 2) * norm.pdf(h, scale=scale)
    y[-1] = y[-1 / 2]
    y[0] = y[0] / 2
    return y


def build_smooth_track(track_handle, out_put_file, scale):
    k0 = get_gauss_kernel0(scale)
    k1 = get_gauss_kernel1(scale)
    k2 = get_gauss_kernel2(scale)
    with h5py.File(out_put_file, 'w') as fout:
        for key, value in track_handle.items():
            hd = fout.create_dataset(key, shape=(3, len(value)))
            s = np.asarray(value)
            for idx, kn in enumerate([k0, k1, k2]):
                print key, idx
                y = np.convolve(s, kn, mode="full")[(3 * scale):(len(s) + 3 * scale)]
                y[:(3 * scale)] = 0
                y[-(3 * scale):] = 0
                hd[idx, :] = y
                del y
                gc.collect()
    return


def find_max_and_min(handle):
    out = {}
    for key, value in handle.items():
        x = np.asarray(value[1, :])
        raw = np.asarray(value[0, :])
        print "read finish."
        cmax = np.where((x[:-1] >= 0) & (x[1:] <= 0) & (~(x[:-1] == 0)))[0]
        cmin = np.where((x[:-1] <= 0) & (x[1:] >= 0) & (~(x[:-1] == 0)))[0]
        print "Find finish."
        vmax = raw[cmax]
        vmin = raw[cmin]
        print "Raw value finish."
        a = [[k, v, 1] for k, v in it.izip(cmax, vmax)] + [[k, v, -1] for k, v in it.izip(cmin, vmin)]
        a.sort(key=lambda u: u[0])
        print "Sort finish."
        out[key] = np.asarray(a)
        print key
    return out
