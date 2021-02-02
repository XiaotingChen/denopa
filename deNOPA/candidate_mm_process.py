# -*- coding: utf-8 -*-
# @Time    : 18-9-4 下午4:03
# @Author  : Matrix
# @Site    : 
# @File    : candidate_mm_process.py
# @Software: PyCharm

from scipy import *
import h5py
import pandas as pd
import gc
import itertools as it
import sys, os
import pysam
from scipy import sparse
from collections import Counter


def filter_mm_candidates(track, peaks, p_over_ctrl, min_sep=15, max_sep=50):
    """
    For each max point in candidate region, determine whether it's value is significantly larger than values 
    not in candidate regions and its distance to nearest edges are appropriate.   
    :param track: site_max_min track with edge points added. 
    :param peaks: candidate regions with index to the sites max min track added. 
    :param p_over_ctrl: The max point is determined as significant when its value is p-value larger than this cutoff.   
    :return: A list of detected candidate nucleosome ends for each candidate regions, where each item is composed by 
        the following items:
            The index of the max point in the sites max min track; 
            The genome position of this max point;
            The genome position of its left edge (-1 if its distance to the edge not lying between min_sep and max_sep);
            The genome position of its right edge.
    """
    outs = []
    for _, peak in peaks.iterrows():
        out = []
        for idx in range(peak["left"], peak["right"]):
            if track[peak[0]][idx, 2] == 1 and track[peak[0]][idx, 3] > p_over_ctrl:
                o = [idx, track[peak[0]][idx, 0], -1, -1]
                try:
                    assert idx > 0
                    if track[peak[0]][idx - 1, 2] == 2:
                        if min_sep <= track[peak[0]][idx, 0] - track[peak[0]][idx - 1, 0] <= max_sep:
                            o[2] = track[peak[0]][idx - 1, 0]
                        else:
                            o[2] = -track[peak[0]][idx - 1, 0]
                except AssertionError:
                    pass
                try:
                    assert idx < track[peak[0]].shape[0] - 1
                    if track[peak[0]][idx + 1, 2] == -2:
                        if min_sep <= track[peak[0]][idx + 1, 0] - track[peak[0]][idx, 0] <= max_sep:
                            o[3] = track[peak[0]][idx + 1, 0]
                        else:
                            o[3] = -track[peak[0]][idx + 1, 0]
                except AssertionError:
                    pass
                if o[2] > 0 or o[3] > 0:
                    out.append(o)
        outs.append(out)
        if _ % 1000 == 0:
            sys.stdout.write("\r%s" % _)
            sys.stdout.flush()
    return outs


def merge_candidate_mms(candidate_mm, chroms, track, min_sep=100, max_sep=200):
    """
    Merge nearby local maximas to bulid candidate nucleosome positions.
    :param candidate_mm: output of filter_mm_candidates.
    :param chroms: the chromosome names of the candidate regions for each set of candidate_mm.
    :param track: Same as the track argument in filter_mm_candidates.
    :param min_sep: min distance between two neighbour local maximas.
    :param max_sep: max distance between two neighbour local maximas.
    :return: A data frame for each candidate nucleosome, each row contains:
        0, chr: The chromosome names of the nucleosome.
        1, start: The start position (local maxima) of the nucleosome.
        2, edge1: The left edge of the nucleosome.
        3, center: Mimina between the two neighbour local maximas.
        4, edge2: The right edge of the nucleosome.
        5, stop: The stop position (local maxima) of the nucleosome.
        6, chridx: Which of candidate region the nucleosome comes from.
    """
    outs = []
    for chridx, (chrname, record) in enumerate(zip(chroms, candidate_mm)):
        for idx, (m1, m2) in enumerate(zip(record[:-1], record[1:])):
            # Originally, the criterion of distances between local maximas and nearby edge points should be satisfied at both directions.
            # Now, it is only needed from one direction. // 2019-09-19
            if ((m1[3] > 0) or (m2[2] > 0)) and (min_sep <= (m2[1] - m1[1]) <= max_sep):
                if m1[3] == -1 or m2[2] == -1:
                    continue
                minv = +inf
                idxmin = 0
                for gdx in range(m1[0] + 1, m2[0]):
                    if track[chrname][gdx][2] == -1 and track[chrname][gdx][1] < minv:
                        idxmin = track[chrname][gdx][0]
                        minv = track[chrname][gdx][1]
                outs.append([chrname, m1[1], abs(m1[3]), idxmin, abs(m2[2]), m2[1], chridx])
    outs = pd.DataFrame(outs)
    return outs


class fragmentEndsMap(object):
    def __init__(self, sam, chrom, start, stop, genome_size, max_frag_len=2000, shift_5p=+4, shift_3p=-5):
        self.chrom = chrom
        self.start = max(0, start - max_frag_len)
        self.stop = min(genome_size[chrom], stop + max_frag_len)
        reads = {}
        for r in sam.fetch(contig=self.chrom, start=self.start, stop=self.stop):
            reads.setdefault(r.query_name, []).append(r)
        reads = {k: v for k, v in reads.iteritems() if len(v) == 2 and
                 v[0].is_read1 ^ v[1].is_read1}
        counts = []
        for r1, r2 in reads.itervalues():
            if r1.is_proper_pair and r2.is_proper_pair:
                p1 = min(r1.reference_start, r2.reference_start) + shift_5p
                if self.start <= p1 <= self.stop:
                    p1 -= self.start
                else:
                    continue
                p2 = max(r1.reference_end, r2.reference_end) + shift_3p
                if self.start <= p2 <= self.stop:
                    p2 -= self.start
                else:
                    continue
                counts.append((int(p1), int(p2)))
        row_ind = []
        col_ind = []
        data = []
        for (k1, k2), v in Counter(counts).iteritems():
            row_ind.append(k1)
            col_ind.append(k2)
            data.append(v)
        self.data = sparse.csr_matrix((data, (row_ind, col_ind)), dtype=float,
                                      shape=(self.stop - self.start + 1, self.stop - self.start + 1))

    def __getitem__(self, item):
        f = lambda u: None if u is None else u - self.start
        item_out = []
        for i in item:
            if isinstance(i, int):
                item_out.append(i - self.start)
            else:
                item_out.append(
                    slice(f(i.start), f(i.stop), f(i.step))
                )
        return self.data.__getitem__(tuple(item_out))


