# -*- coding: utf-8 -*-
# @Time    : 2019/6/10 15:26
# @Author  : Matrix
# @Site    :
# @File    : dbscan_model.py
# @Software: PyCharm
from __future__ import print_function
import sys
import os
from scipy import *
from sklearn.neighbors import KDTree
from sklearn.cluster import DBSCAN
import pandas as pd
from scipy.stats import norm
from scipy.linalg import eigh
from collections import Counter
from statsmodels.stats.multitest import multipletests
import h5py
import pickle as pk
from numpy import *


class FinalModel(object):
    def __init__(
        self, adjusted_p_threshold=0.1, dbscan_only=True, q_eps=0.9, q_mps=0.1
    ):
        self.adjusted_p = adjusted_p_threshold
        self.qeps = q_eps
        self.qmps = q_mps
        self.dbscan_only = dbscan_only

    def loadFeatures(self, data):
        d = data.apply(lambda v: log(abs(v[5] - v[1] - 156) / 157.0 + 1), axis=1)

        def fun(x):
            return (x <= 0) * ((2 * norm(0, 1).pdf(0)) * x) + (x > 0) * 2 * (
                norm(0, 1).cdf(x) - 0.5
            )

        def AT(a, b):
            m = b[0]
            return -(2 * sqrt(a + 0.375) - (2 * sqrt(m + 0.375) - 1 / (4 * sqrt(m))))

        x = asarray([AT(a, b) for a, b in zip(data[10], data[11])])
        x = x - mean(x)
        s = sum(x[x > 0] ** 2) / sum(x > 0)
        x = fun(x / sqrt(s))
        """
            z = asarray([AT(a, b) for a, b in zip(data[14], data[15])])
            z = mean(z) - z
            s = sum(z[z > 0] ** 2) / sum(z > 0)
            z = fun(z / sqrt(s))
        """

        def mu(n, p):
            c = sqrt(n + 0.5)
            x = n * p
            g = (2 * x - n) / (n + 0.75)
            f0 = c * arcsin(g)
            f2 = g / (1 - g**2) ** (1.5) * c
            g1 = 2 / (n + 0.75)
            return f0 + 0.5 * f2 * g1**2 * n * p * (1 - p)

        def BAT(x, n, p):
            return sqrt(n + 0.5) * arcsin((2 * x - n) / (n + 0.75)) - mu(n, p)

        y = asarray([BAT(a, b, c[1]) for a, b, c in zip(data[8], data[7], data[9])])
        y = y - mean(y)
        s = sum(y[y > 0] ** 2) / sum(y > 0)
        y = fun(y / sqrt(s))

        self.rawData = asarray([d, x, y]).T

    def transformFeatures(self):
        V = cov(self.rawData, rowvar=False)
        A, T = eigh(V)
        U = T.dot(diag(A ** (-0.5)))
        self.transData = (self.rawData - self.rawData.mean(axis=0)).dot(U)

    def choosePara(self):
        kdtree = KDTree(self.transData)
        dis, _ = kdtree.query(self.transData, k=2)
        self.eps = percentile(dis[:, 1], self.qeps * 100.0)
        nbs = kdtree.query_radius(self.transData, self.eps, count_only=True) - 1
        self.minpts = int(max(nbs) * self.qmps) + 1

    def detectOutliers(self):
        dbscan = DBSCAN(eps=self.eps, min_samples=self.minpts, algorithm="kd_tree")
        l = dbscan.fit_predict(self.transData)
        c = Counter(l)
        d = max(c.items(), key=lambda u: u[1])[0]
        self.label = l == d

    def detectSig(self, rawData):
        rawData[12] = [
            1 - (1 - p[1]) * (1 - q[2]) for p, q in zip(rawData[11], rawData[9])
        ]
        self.sig = multipletests(rawData[12], method="fdr_bh", alpha=self.adjusted_p)[0]

    def __call__(self, rawData, label=None):
        # rawData = rawData.loc[rawData[14] > 0, :].copy()
        if label == None:
            self.loadFeatures(rawData)
            self.transformFeatures()
            self.choosePara()
            self.detectOutliers()
            self.detectSig(rawData)
            label = self.label | (
                zeros_like(self.label, dtype=bool) if self.dbscan_only else self.sig
            )
        else:
            label = label
        x = []
        for _, v in rawData.loc[label, :].copy().iterrows():
            x.append(
                [
                    v[0],
                    int(v[1]),
                    int(v[5] + 1),
                    int(v[2]),
                    int(v[3]),
                    int(v[4]) + 1,
                    v[7],
                    v[8],
                    v[10],
                    v[9][2],
                    v[11][1],
                    v[12],
                    v[13],
                ]
            )
        x = pd.DataFrame(x)
        x.index = arange(rawData.shape[0])[label]
        x.columns = range(x.shape[1])
        for i in range(1, 7):
            x[i] = x[i].astype(int)
        return x


class CandidatesWithNOC(object):
    def __init__(self, candidates, nocFile, nosFile, trackName="Raw"):
        """
        Add the nucleosome occupancy signals to the candidate nucleosomes and make decision based on the number of fragments covering
        the inner, the number of fragments initiating from the inner and the maximum nucleosome occupancy.
        :param candidates: DataFrame of full nucleosome candidates or pickle file name from which it can be loaded.
        :param nocFile: nucleosome occupancy signal file generated by nucleosomeOccupancyTrack function in fragmentLengthsDist
        :param trackName: Name of the track used in the signal track file.
        """
        if isinstance(candidates, pd.DataFrame):
            self.cand = candidates.copy()
        else:
            with open(candidates) as fin:
                self.cand = pk.load(fin)
        self.nocFile = nocFile
        self.nosFile = nosFile
        self.track = trackName

    def addToCand(self):
        def aPvalue(x, n):
            an = (2 * log(n)) ** (-0.5)
            bn = (2 * log(n)) ** 0.5 - 0.5 * ((2 * log(n)) ** (-0.5)) * (
                log(log(n)) + log(4 * pi)
            )
            t = (x - bn) / an
            return exp(-exp(-t))

        y = []
        with h5py.File(self.nocFile, "r") as hdf:
            for ix, (c, a, b) in enumerate(
                zip(
                    self.cand[0], self.cand[2].astype(int), self.cand[4].astype(int) + 1
                )
            ):
                # u = list(hdf[self.track][c][aa:a]) + list(hdf[self.track][c][b:bb])
                t = hdf[self.track][c][a:b]
                mx = max(t)
                w = where(t == mx)[0]
                w = w[int(len(w) / 2)] + a
                y.append([w, mx, (mx - mean(t)) / std(t), len(t)])
        y = asarray(y)
        flag = (isnan(y)).sum(axis=1) == 0
        y = y[flag, :].copy()
        self.cand = self.cand.loc[flag, :].copy()
        self.cand[3] = y[:, 0].astype(int)
        self.cand[14] = [[i[1], i[2], aPvalue(i[2], i[3]), i[3]] for i in y]
        self.cand[3] = y[:, 0]
        self.cand[12] = [
            1 - (1 - p2[1]) * (1 - p3[2])
            for p1, p2, p3 in zip(self.cand[9], self.cand[11], self.cand[14])
        ]

    def getOutPut(self, alpha=0.1):
        flag = asarray([i[-2] > alpha for i in self.cand[14]])
        print("%d / %d" % (sum(flag), len(flag)))
        y = self.cand.loc[flag, :].copy()
        return y.sort_values([0, 1, 2])


class varCandidatesWithNOC(object):
    def __init__(self, candidates, nocFile, nosFile, trackName="Raw"):
        """
        Add the nucleosome occupancy signals to the candidate nucleosomes and make decision based on the number of fragments covering
        the inner, the number of fragments initiating from the inner and the maximum nucleosome occupancy.
        :param candidates: DataFrame of full nucleosome candidates or pickle file name from which it can be loaded.
        :param nocFile: nucleosome occupancy signal file generated by nucleosomeOccupancyTrack function in fragmentLengthsDist
        :param trackName: Name of the track used in the signal track file.
        """
        if isinstance(candidates, pd.DataFrame):
            self.cand = candidates.copy()
        else:
            with open(candidates) as fin:
                self.cand = pk.load(fin)
        self.nocFile = nocFile
        self.nosFile = nosFile
        self.track = trackName

    def addToCand(self):
        def aPvalue(x, n):
            an = (2 * log(n)) ** (-0.5)
            bn = (2 * log(n)) ** 0.5 - 0.5 * ((2 * log(n)) ** (-0.5)) * (
                log(log(n)) + log(4 * pi)
            )
            t = (x - bn) / an
            return exp(-exp(-t))

        y = []
        flag = []
        with h5py.File(self.nocFile, "r") as hdf, h5py.File(self.nosFile, "r") as sdf:
            for ix, (c, a, b) in enumerate(
                zip(
                    self.cand[0], self.cand[2].astype(int), self.cand[4].astype(int) + 1
                )
            ):
                s0 = asarray(hdf[self.track][c][a:b])
                s1 = asarray(sdf[self.track]["1"][c][a:b])
                s2 = asarray(sdf[self.track]["2"][c][a:b])
                f = where(
                    (s1[:-1] >= 0) & (s1[1:] <= 0) & (s2[:-1] < 0) & (s2[1:] < 0)
                )[0]
                if len(f) > 0:
                    s = s0[f]
                    p = a + f[argmax(s)]
                    q = max(s0)
                    y.append(
                        [p, q, aPvalue((q - mean(s0)) / std(s0), len(s0)), len(s0)]
                    )
                    flag.append(True)
                else:
                    flag.append(False)
        y = asarray(y)
        self.cand = self.cand.loc[flag, :].copy()
        self.cand[3] = y[:, 0].astype(int)
        self.cand[14] = [[i[1], i[2], i[3]] for i in y]
        # self.cand[12] = [1 - (1 - p1[2]) * (1 - p2[1]) * (1 - p3[2]) for p1, p2, p3 in
        #                  zip(self.cand[9], self.cand[11], self.cand[14])]

    def getOutPut(self, alpha=0.1):
        flag = asarray([i[1] > alpha for i in self.cand[14]]) & (self.cand[8] > 0)
        print("%d / %d" % (sum(flag), len(flag)))
        y = self.cand.loc[flag, :].copy()
        return y.sort_values([0, 1, 2])
