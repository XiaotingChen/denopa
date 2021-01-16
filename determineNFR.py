# -*- coding: utf-8 -*-
# @Time    : 2019/10/29 10:57
# @Author  : Matrix
# @Site    : 
# @File    : determineNFR.py
# @Software: PyCharm

from scipy import *
from scipy import stats
from statsmodels.stats.multitest import multipletests
import h5py
import pandas as pd
import signal_track_builder
import call_peak


class NFRDetection(object):
    def __init__(self, numReads, shortTrack, qvalue, trackName="short"):
        self.numReads = numReads.copy().sort_values([0, 1, 5])
        self.qvalue = qvalue
        with h5py.File(shortTrack, 'r') as hdf:
            self.mms = signal_track_builder.MakeMaxMinTrack(
                hdf["%s/0" % trackName],
                hdf["%s/1" % trackName],
                hdf["%s/2" % trackName]
            )()
        self.vmaxs = {}

    def get_pmax(self):
        self.vmaxs = {}
        for k, data in self.mms.items():
            pmax = call_peak.get_para(data[data[:, 2] > 0, 1])
            pvalues = asarray([stats.gamma.sf(data[data[:, 2] > 0, 1], pmax[0], scale=pmax[1])]).T
            self.vmaxs[k] = append(
                data[data[:, 2] > 0, :][:, [0, 1]],
                pvalues,
                axis=1
            )

    def __call__(self, *args, **kwargs):
        if not self.vmaxs:
            self.get_pmax()
        y = []
        for chrom in set(self.numReads[0]):
            a = self.numReads.loc[self.numReads[0] == chrom, :].iloc[:-1, :]
            b = self.numReads.loc[self.numReads[0] == chrom, :].iloc[1:, :]
            a6 = asarray(a[6])
            a5 = asarray(a[5])
            b1 = asarray(b[1])
            b6 = asarray(b[6])
            f = (a6 == b6) & (b1 - a5 > 75)
            a5 = a5[f].copy()
            b1 = b1[f].copy()
            f1 = searchsorted(a5, self.vmaxs[chrom][:, 0], side="left")
            f2 = searchsorted(b1, self.vmaxs[chrom][:, 0], side="right")
            f = (f1 - f2 == 1)
            z = {}
            for i, j, v in zip(a5[f2[f]], b1[f2[f]], self.vmaxs[chrom][f, :]):
                z.setdefault((i, j), []).append(v)
            z = {k: list(max(v, key=lambda t: t[1])) for k, v in z.items()}
            y.extend(
                [
                    [chrom] + list(k) + list(v) for k, v in z.items()
                ]
            )
        y = pd.DataFrame(y)
        z = [
            list(y[0]),
            list(y[1] + 1),
            list(y[2]),
            ["ATAC_%d" % (1 + k) for k in range(y.shape[0])],
            [0 for k in range(y.shape[0])],
            ["." for k in range(y.shape[0])],
            list(y[4]),
            list(y[5]),
            list(multipletests(y[5], self.qvalue, method="fdr_bh")[1]),
            list(y[3] - (y[1] + 1))
        ]
        z = pd.DataFrame(z).T
        z[1] = z[1].astype(int)
        z[2] = z[2].astype(int)
        z[9] = z[9].astype(int)
        return z.loc[z[8] <= self.qvalue, :].copy()
