# -*- coding: utf-8 -*-
# @Time    : 19-1-3 上午9:41
# @Author  : Matrix
# @Site    :
# @File    : ocsvm_model.py
# @Software: PyCharm
from __future__ import print_function
import os, sys
from scipy import *
import itertools as it
from scipy.stats import binom, poisson, norm
import pysam
import pandas as pd
import multiprocessing as mp
import time
from collections import Counter
from sklearn.svm import OneClassSVM
import logging
from . import determine_dynamic
from numpy import *

logging.basicConfig(format="%(asctime)s: %(message)s", level=logging.ERROR)


class testing(object):
    def __init__(self, func=None):
        self.times = []
        self.idx = 0
        self.params = []
        self.func = func

    def __call__(self, *args, **kwargs):
        if self.func == None:
            return testing(func=args[0])
        else:
            self.params.append((args, kwargs))
            t = time.time()
            y = self.func(*args, **kwargs)
            self.times.append(time.time() - t)
            self.idx += 1
            sys.stdout.write("\r%d" % self.idx)
            sys.stdout.flush()
            return y

    def clear(self):
        self.params = []
        self.idx = 0
        self.times = []


def __reads_in_peaks(
    sam_list,
    peak_chr,
    peak_start,
    peak_stop,
    frag_len,
    left_shift=+4,
    right_shift=-5,
    return_ends=True,
):
    """
    Read each sam file in sam_list and find fragments overlapped with given region denoted by peak_chr:peak_start-peak_stop,
        then report the positions of these fragments (shifted by left_shift and right_shift) and their shifted ends.
    :return: reads: A list of the shifted position (left most, right most) of each fragment.
             read_ends: Cutting positions (5' end) of overlapped fragments.
    """
    reads = []
    if return_ends:
        read_ends = []
    for sam_file in sam_list:
        with pysam.AlignmentFile(sam_file) as sam:
            genome_size = dict(zip(sam.references, sam.lengths))
            pos_start = max(0, peak_start - frag_len - 1)
            pos_stop = min(genome_size[peak_chr], peak_stop + frag_len + 1)
            rds = {}
            for r2 in sam.fetch(contig=peak_chr, start=pos_start, stop=pos_stop):
                try:
                    r1 = rds.pop(r2.query_name)
                    r1, r2 = r1, r2 if (r2.is_reverse) else (r2, r1)
                    p1 = r1.reference_start + left_shift
                    p2 = r2.reference_end + right_shift - 1
                    if min(p2, peak_stop) + 1 - max(p1, peak_start) >= 0:
                        reads.append((p1, p2))
                        if return_ends:
                            if peak_start <= p1 <= peak_stop:
                                read_ends.append(p1)
                            if peak_start <= p2 <= peak_stop:
                                read_ends.append(p2)
                except KeyError:
                    rds[r2.query_name] = r2
            del rds
    reads = asarray(sorted(reads), dtype=double)
    if return_ends:
        read_ends = asarray(sorted(read_ends))
        return reads, read_ends
    else:
        return reads


# @testing()
def __calc_ov_frags(para):
    """
    Core function calculating the number of reads overlapped for each candidate nucleosome.
    :param para: tuple including the needed parameters as follows:
        :cand_mms: A data frame containing positions of candidate nucleosomes coming from the same peak with the same
            format as the function cand_mm_process.merge_candidate_mms gives.
        :peak: The peak where the candidate nucleosomes come from.
        :sam_list: The sam/bam file list where the reads comes from.
        :frag_len: Maximum fragment length of the sequenced reads.
        :left_shift/right_shift: same as pileup_signals.build_signal_track.
        :sommth_file: file name of smoothing tracks.
        :dynamic_p: pvalue for dynamic determination.
    :return: A copy of cand_mms with following columns added:
        column 7: number of fragments covering the local minimar in the peak;
        column 8: number of fragments covering the whole peak.
        column 9: a list of information about fragment coverage of the peak including:
            [
                the expected probability of fragments with length >= the length of the peak.
                the expected propability of a fragment intersecting the peaks which will cover the peak.
                the p value of the the number of fragments covering the peak is significantly large.
            ]
        column 10: number of reads starting (5' end) between the two edges of the peak.
        column 11: information of the number of reads starting between the two edges of the peak, including:
            [
                the expected number,
                the p value of the observed number of reads significantly simall.
            ]
        column 12: combined p valued of the fragment coverage of the peak and number of read ends from it.
        column 13: if the peak is dynamic.
    """
    logger = logging.getLogger(__name__)
    cand_mms, peak, sam_list, frag_len, left_shift, right_shift, smooth_file, dp = para
    peak_idx = list(set(cand_mms[6]))
    assert len(peak_idx) == 1, "The number of peaks should be exactly 1."
    peak_idx = peak_idx[0]
    peak_chr, peak_start, peak_stop = tuple(peak[:3])
    cand_cp = cand_mms.copy()

    # load read positions in this condidate region.
    reads, read_ends = __reads_in_peaks(
        sam_list,
        peak_chr,
        peak_start,
        peak_stop,
        frag_len,
        left_shift=left_shift,
        right_shift=right_shift,
    )

    # Initialize to determine if the peak is dynamic.
    isDynamic = determine_dynamic.DetermineDynamics(
        reads, smooth_file, (peak_chr, peak_start, peak_stop), dp
    )

    def fun(cand_pos, pos):
        """
        Find rows a from cand_pos satisfying a[0] <= pos[0] and a[1] >= pos[1] at the same time.
        :param cand_pos: n x 2 array sorted by its first column.
        :param pos: list with 2 element, denoting the position being overlapped.
        :return: array of found rows, sorted by its first column.
        """
        cand_pos = cand_pos.copy()
        y = cand_pos[: searchsorted(cand_pos[:, 0], pos[0], side="right"), :]
        y = asarray(sorted([list(i) for i in y], key=lambda a: a[1]))
        z = y[searchsorted(y[:, 1], pos[1], side="left") :, :]
        return asarray(sorted(z, key=lambda a: a[0]))

    # calculate expect values in this condidate regions.
    peak_lengths = pd.Series(Counter(reads[:, 1] - reads[:, 0] + 1))
    peak_length_dis = peak_lengths / double(sum(peak_lengths))
    mean_end_per_base = len(read_ends) / float(peak[2] - peak[1] + 1)

    # Calculate the observed values.
    outs = [[], []]
    for _, record in cand_cp.iterrows():
        try:
            x = fun(reads, [record[3], record[3]])
        except IndexError:
            outs[0].append(0)
            outs[1].append(0)
            continue
        try:
            y = fun(x, [record[1], record[5]])
        except IndexError:
            outs[0].append(x.shape[0])
            outs[1].append(0)
            continue
        outs[0].append(x.shape[0])
        outs[1].append(y.shape[0])
    outs = pd.DataFrame(outs).T
    outs.index = cand_cp.index
    cand_cp = pd.concat([cand_cp, outs], join="inner", axis=1)
    cand_cp.columns = range(cand_cp.shape[1])

    def coverage_pvalue(record):
        idx_a = cand_cp.shape[1] - 2
        idx_c = cand_cp.shape[1] - 1
        L = record[5] - record[1] + 1
        p0 = sum(peak_length_dis[peak_length_dis.index >= L])
        p = sum(
            [
                pl * double(l - L + 1) / double(l)
                for l, pl in peak_length_dis.items()
                if l >= L
            ]
        )
        return [p0, p, binom(record[idx_a], p).sf(record[idx_c])]

    cand_cp[cand_cp.shape[1]] = cand_cp.apply(coverage_pvalue, axis=1)
    cand_cp[cand_cp.shape[1]] = searchsorted(
        read_ends, cand_cp[4], side="right"
    ) - searchsorted(read_ends, cand_cp[2], side="left")

    def end_pvalue(record):
        idx = cand_cp.shape[1] - 1
        lmd = mean_end_per_base * (record[4] - record[2] + 1)
        return [
            lmd,
            poisson(mean_end_per_base * (record[4] - record[2] + 1)).cdf(
                record[idx] - 1
            ),
        ]

    cand_cp[cand_cp.shape[1]] = cand_cp.apply(end_pvalue, axis=1)
    cand_cp[cand_cp.shape[1]] = cand_cp.apply(
        lambda u: 1 - (1 - u[9][2]) * (1 - u[11][1]), axis=1
    )

    # Determine if the peak is dynamic.
    cand_cp[cand_cp.shape[1]] = [
        isDynamic((v[0], v[1], v[3], v[5])) for _, v in cand_cp.iterrows()
    ]
    logger.info("Peak %d processed." % peak_idx)
    return cand_cp


def calc_ov_frags(
    sam_list,
    cand_mm,
    peaks,
    frag_len,
    left_shift,
    right_shift,
    smooth_file,
    dynamic_pvalue=0.05,
    proc=1,
):
    pool = mp.Pool(proc)
    cand_split = {}
    for _, record in cand_mm.iterrows():
        cand_split.setdefault(record[6], []).append(record)
    params = []
    for k, v in cand_split.items():
        params.append(
            (
                pd.DataFrame(v),
                peaks.loc[k, :],
                sam_list,
                frag_len,
                left_shift,
                right_shift,
                smooth_file,
                dynamic_pvalue,
            )
        )
    y = (pool.map if proc > 1 else map)(__calc_ov_frags, params)
    pool.close()
    pool.join()
    y = pd.concat(y, axis=0, join="inner")
    y.index = range(y.shape[0])
    return y


class FinalModeling(object):
    def __init__(self, raw_data, output_prefix, pvalue=0.05, fraction=0.1):
        super(FinalModeling, self).__init__()
        self.raw_data = raw_data
        self.features = None
        self.pos = []
        self.neg = []
        self.pvalue = pvalue
        self.fraction = fraction
        self.outPrefix = output_prefix
        self.output = []

    def getFeatures(self):
        raw = self.raw_data.loc[:, [0, 1, 5]]
        raw.columns = range(3)
        raw[1] = asarray(raw[1] + 1, dtype=int)
        raw[2] = asarray(raw[2] + 1, dtype=int)
        raw.to_csv("%s_raw.bed" % self.outPrefix, header=None, sep="\t", index=None)
        features = []
        index = []
        for key, value in self.raw_data.iterrows():
            if value[8] == 0:
                self.neg.append(key)
            else:
                if value[10] == 0:
                    if value[9][2] < self.pvalue:
                        self.pos.append(key)
                    else:
                        self.neg.append(key)
                else:
                    if value[12] <= self.pvalue:
                        index.append(key)
                        features.append(
                            [
                                log(abs(value[5] - value[1] - 156) / 157.0 + 1),
                                log(float(value[8]) / value[7] / value[9][1]),
                                log(value[10] / value[11][0]),
                            ]
                        )
                    else:
                        self.neg.append(key)
        self.raw_features = pd.DataFrame(features, index=index)
        neg = set(self.neg)
        naive = self.raw_data.loc[
            [i for i in self.raw_data.index if not i in neg], [0, 1, 5]
        ].copy()
        naive[1] = asarray(naive[1] + 1, dtype=int)
        naive[5] = asarray(naive[5] + 1, dtype=int)
        naive.to_csv("%s_naive.bed" % self.outPrefix, header=None, sep="\t", index=None)

    def featureTransform(self):
        transFun = fun = lambda x: (x <= 0) * ((2 * norm(0, 1).pdf(0)) * x) + (
            x > 0
        ) * 2 * (norm(0, 1).cdf(x) - 0.5)
        self.transformed_features = pd.DataFrame(
            index=self.raw_features.index, columns=self.raw_features.columns
        )
        self.transformed_features[0] = self.raw_features[0] / std(self.raw_features[0])
        self.transformed_features[1] = transFun(
            (self.raw_features[1] - mean(self.raw_features[1]))
            / std(self.raw_features[1])
        )
        self.transformed_features[1] = self.transformed_features[1] / std(
            self.transformed_features[1]
        )
        self.transformed_features[2] = transFun(
            -(self.raw_features[2] - mean(self.raw_features[2]))
            / std(self.raw_features[2])
        )
        self.transformed_features[2] = self.transformed_features[2] / std(
            self.transformed_features[2]
        )

    def svmFilter(self):
        svm = OneClassSVM(nu=self.fraction, kernel="rbf").fit(self.transformed_features)
        y = svm.predict(self.transformed_features)
        for a, (b, c) in zip(y, self.transformed_features.iterrows()):
            if a > 0:
                self.pos.append(b)
            else:
                if c[1] > 0 and c[2] > 0:
                    self.pos.append(b)
                else:
                    self.neg.append(b)
        self.pos.sort()
        self.neg.sort()

    def outputToFile(self):
        name = os.path.split(self.outPrefix)[-1]
        pos = self.raw_data.loc[self.pos, [0, 1, 5]].copy()
        pos.columns = range(3)
        pos[1] = asarray(pos[1] + 1, dtype=int)
        pos[2] = asarray(pos[2] + 1, dtype=int)
        pos[3] = map(lambda u: "%s_%d" % (name, u), pos.index)
        pos[4] = asarray(self.raw_data.loc[self.pos, 3] + 1, dtype=int)
        pos[5] = self.raw_data.loc[self.pos, 8]
        pos[6] = self.raw_data.loc[self.pos, 10]
        pos[7] = map(lambda u: u[2], self.raw_data.loc[self.pos, 9])
        pos[8] = map(lambda u: u[1], self.raw_data.loc[self.pos, 11])
        pos[9] = self.raw_data.loc[self.pos, 12]
        pos.to_csv("%s_peaks.bed" % self.outPrefix, header=None, index=None, sep="\t")

    def run(self):
        self.getFeatures()
        self.featureTransform()
        self.svmFilter()
        self.outputToFile()
