# -*- coding: utf-8 -*-
# @Time    : 2019/5/24 16:02
# @Author  : Matrix
# @Site    : 
# @File    : determine_dynamic.py
# @Software: PyCharm
from __future__ import print_function
from scipy import *
import h5py
from scipy.stats import binom


class DetermineDynamics(object):
    def __init__(self, reads, smooth_file, cand_region, pvalue):
        """
        Initialize for determination of dynamics of each peak in a candidate, this will be embedded with ocsvm_model.calc_ov_frags
        :param reads: Reads in the candidate regions. Same as reads in ocsvm_model.calc_ov_frags.
        :param smooth_file: File path of the smoothed tracks.
        :param cand_region: The (chrom,start,stop) of the candidate region where the test would be implemented.
        :param pvalue: Pvalue used to determine whether reads are directed.
        """
        super(DetermineDynamics, self).__init__()
        self.chr, self.start, self.stop = cand_region
        self.pvalue = pvalue
        self.reads = reads - self.start

        # Load used tracks.
        with h5py.File(smooth_file, 'r') as hdf:
            self.s0 = asarray(hdf["sites/0"][self.chr][self.start:(self.stop + 1)])
            self.s1 = asarray(hdf["sites/1"][self.chr][self.start:(self.stop + 1)])
            self.s2 = asarray(hdf["sites/2"][self.chr][self.start:(self.stop + 1)])
            self.right_edges = where(
                (self.s2[:-1] <= 0) & (self.s2[1:] >= 0)
            )[0]
            self.left_edges = where(
                (self.s2[:-1] >= 0) & (self.s2[1:] <= 0)
            )[0]

    def __call__(self, peak):
        """
        Determine if the peak is dynamic.
        :param peak: tuple representing a peak as (chromname, start, local min, stop).
        :return: One of "NoDynam" (not dynamic), "NoSigDynam" (there is a hidden cut peak but reads from the peak are
            not significantly coincided with the nucleosome), "SigDynam" (there is a hidden cut peak but reads from the
            peak are significantly coincided with the nucleosome).
        """
        chrom, start, mid, stop = peak
        start = start - self.start
        mid = mid - self.start
        stop = stop - self.start
        r1 = "Main"
        r2 = "Main"

        # If there is dynamic regions in left side.
        edges = [i for i in self.right_edges if start <= i and i <= mid]
        if len(edges) > 1:
            a, b = edges[:2]
            p = [i for i in self.left_edges if a <= i and i <= b]
            if len(p) > 0:
                p = p[0]
                if stop - p > 105:
                    r1 = "NoSig"
                    p1 = sum((a <= self.reads[:, 0]) & (self.reads[:, 0] <= b))
                    p2 = sum((a <= self.reads[:, 1]) & (self.reads[:, 1] <= b))
                    if binom(p1 + p2, 0.5).sf(p1) < self.pvalue:
                        r1 = "Sig"

        # If there is dynamic regions in right side.
        edges = [i for i in self.left_edges if mid <= i and i <= stop]
        if len(edges) > 1:
            a, b = edges[-2:]
            p = [i for i in self.right_edges if a <= i and i <= b]
            if len(p) > 0:
                p = p[-1]
                if p - start > 105:
                    r2 = "NoSig"
                    p1 = sum((a <= self.reads[:, 0]) & (self.reads[:, 0] <= b))
                    p2 = sum((a <= self.reads[:, 1]) & (self.reads[:, 1] <= b))
                    if binom(p1 + p2, 0.5).sf(p2) < self.pvalue:
                        r2 = "Sig"

        if r1 == "Sig" or r2 == "Sig":
            return "SigDynam"
        if r1 == "NoSig" or r2 == "NoSig":
            return "NoSigDynam"
        return "NoDynam"
