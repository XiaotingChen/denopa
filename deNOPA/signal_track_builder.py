# -*- coding: utf-8 -*-
# @Time    : 18-5-12 下午1:47
# @Author  : Matrix
# @Site    :
# @File    : signal_track_builder.py
# @Software: PyCharm
from __future__ import print_function
from scipy import *
from scipy.stats import norm
import h5py
import pandas as pd
import gc
import itertools as it
import multiprocessing as mp
import logging

logging.basicConfig(format="%(asctime)s: %(message)s", level=logging.INFO)


def gauss_convolve(track_in, track_out):
    lock = mp.Lock()
    proc1 = GaussConvolve(track_in, track_out, "coverage", 72, lock=lock)
    proc2 = GaussConvolve(track_in, track_out, "sites", 24, lock=lock)
    proc1.start()
    proc2.start()
    proc1.join()
    proc2.join()
    assert proc1.exitcode == 0, "Coverage smoothed track created unsuccessfully!"
    assert proc2.exitcode == 0, "Site smoothed track created unsuccessfully!"
    returnF


class GaussConvolve(mp.Process):
    def __init__(self, track_in, track_out, track_name, scale, lock=None, third_dev=True):
        super(GaussConvolve, self).__init__()
        self.logger = logging.getLogger(__name__)
        self.track_in = track_in
        self.track_out = track_out
        self.track_name = track_name
        self.scale = int(scale)
        self.third_dev = third_dev
        self.lock = lock

    def get_gauss_kernel(self, d=0):
        scale = float64(self.scale)
        h = arange(-3 * int(scale), 3 * int(scale) + 1, dtype=float64)
        y = {
            0: lambda: norm.pdf(h, scale=scale),
            1: lambda: -1 / (scale ** 2) * h * norm.pdf(h, scale=scale),
            2: lambda: 1 / (scale ** 4) * (h ** 2 - scale ** 2) * norm.pdf(h, scale=scale),
            3: lambda: (3 * h / scale ** 4 - h ** 3 / scale ** 6) * norm.pdf(h, scale=scale)
        }[d]()
        y[0] /= 2
        y[-1] /= 2
        return y

    def run(self):
        if self.track_name == "coverage":
            kernels = [self.get_gauss_kernel(d) for d in range(3)]
        else:
            if self.third_dev:
                kernels = [self.get_gauss_kernel(d) for d in range(4)]
            else:
                kernels = [self.get_gauss_kernel(d) for d in range(3)]
        with h5py.File(self.track_in, 'r') as tin:
            for idx, g in enumerate(kernels):
                for k, v in tin[self.track_name].items():
                    v = asarray(v)
                    y = convolve(v, g, mode="full")
                    y = y[(3 * self.scale):(len(v) + 3 * self.scale)]
                    y[:(6 * self.scale + 1)] = 0
                    y[-(6 * self.scale + 1):] = 0
                    tname = "%s/%d/%s" % (self.track_name, idx, k)
                    self.logger.info("%s track of %s at %d level has been created. " % (
                        self.track_name, k, idx))
                    if self.lock == None:
                        with h5py.File(self.track_out, 'a') as tout:
                            tout.create_dataset(tname, data=y)
                    else:
                        with self.lock:
                            with h5py.File(self.track_out, 'a') as tout:
                                tout.create_dataset(tname, data=y)


class MakeMaxMinTrack(object):
    def __init__(self, track0, track1, track2):
        """
        Make a track of max and min values from smoothed values tracks.
        :param track0: track handle of the values.
        :param track1: track handle of the differential values.
        :param track2: track handle of the second differential values.
        """
        self.tracks = [track0, track1, track2]
        self.values = {k: [] for k in self.tracks[0].keys()}
        self.currentChrom = ""

    def loadChrom(self, chrom):
        """
        Make the chrome for a single chromosome. For each local maxima, local minimas nearest to it are recorded.
        :param chrom: Name of the chromosome.
        """
        self.currentChrom = chrom
        v = [asarray(i[chrom]) for i in self.tracks]
        values = []
        self.__getMax(values, *v)
        self.__getMin(values, *v)
        del v
        gc.collect()
        values.sort()
        a, b = (0, len(values) - 1)
        while (values[a][2] == 1):
            a += 1
        while (values[b][2] == 1):
            b -= 1
        values = values[a:(b + 1)]
        for ix in where([i[2] == 1 for i in values])[0]:
            ia = ix - 1
            ib = ix + 1
            while values[ia][2] == 1:
                print("There are not local minimar between two maximars")
                ia -= 1
            while values[ib][2] == 1:
                print("There are not local minimar between two maximars")
                ib += 1
            self.values[chrom].extend(
                [values[ia], values[ix], values[ib]]
            )
        self.values[chrom].sort()

    def __getMax(self, vout, v0, v1, v2):
        flags = where(
            (v1[:-1] >= 0) & (v1[1:] <= 0) & (
                (v2[:-1] < 0) & (v2[1:] < 0)
            )
        )[0]
        x = [[flags[0]]]
        for i in flags[1:]:
            if i - x[-1][-1] == 1:
                x[-1].append(i)
            else:
                x.append([i])
        for i in x:
            if len(i) == 0:
                idx = i[0] if v0[i[0]] > v0[i[0] + 1] else i[0] + 1
            else:
                idx = max(i, key=lambda u: v0[u])
            vout.append(
                (idx, v0[idx], +1)
            )

    def __getMin(self, vout, v0, v1, v2):
        for i in where(
                (v1[:-1] <= 0) & (v1[1:] >= 0) & (
                    ~((v2[:-1] == 0) & (v2[1:] == 0))
                )
        )[0]:
            idx = i if v0[i] < v0[i + 1] else i + 1
            vout.append(
                (idx, v0[idx], -1)
            )

    def __call__(self):
        """
        Make the tracks for the whole genome.
        :return: A dict whose keys are the chromosome names and values the max_min track for each chromosome as following
        columns:
            0: the position of the point (0 based);
            1: the signal value of the point;
            2: whether the point is a local maximar (+1) or a local minimar (-1);
        """
        self.values = {k: [] for k in self.tracks[0].keys()}
        for k in self.values.keys():
            self.loadChrom(k)
            self.values[k] = sorted(list(set(self.values[k])))
        return {k: asarray(v) for k, v in self.values.items()}


# Delete this if MakeMaxMinTrack is OK.
def make_max_min_track(track, diff_track):
    """
    Find the local maximum and minimum points of a signal track.
    :param track: The raw (smoothed) signal track.
    :param diff_track: The differential of the raw signal tracks.
    :return: A dict whose keys are the chromosome names and values the max_min track for each chromosome as following
        columns:
        0: the position of the point (0 based);
        1: the signal value of the point;
        2: whether the point is a local maximar (+1) or a local minimar (-1);
    """
    df = {}
    for k, v in track.items():
        value = asarray(v)
        diff_value = asarray(diff_track[k])
        target_points = where(diff_value[1:] * diff_value[:-1] <= 0)[0]
        flag = ~((diff_value[target_points] == 0) &
                 (diff_value[target_points + 1] == 0))
        target_points = target_points[flag]
        target_max = target_points[diff_value[target_points] >= 0]
        target_min = target_points[diff_value[target_points + 1] >= 0]
        target_max = [[i, value[i], +1] for i in target_max if
                      diff_value[i] * diff_value[i + 1] < 0 or (diff_value[i] > 0 or diff_value[i + 1] < 0)]
        target_min = [[i, value[i], -1] for i in target_min if
                      diff_value[i] * diff_value[i + 1] < 0 or (diff_value[i] < 0 or diff_value[i + 1] > 0)]
        t = target_min + target_max
        t.sort(key=lambda u: u[0])
        df[k] = asarray(t)
        del value
        del diff_value
        gc.collect()
    return df


def split_max_min_into_peaks(max_min_track, peaks):
    """
    Add the annotation to each peaks decoding form which index of the max_min_track to which index the positions of the
    max_min_tracks lies in it.
    :param max_min_track: a max min track returned by make_max_min_track function.
    :param peaks: a bed like data frame (such as returned by call_peak.call_candidate_regions);
    :return: a copy of the peaks with two columns added:
        left: the index of the max_min_track from which the peaks begins;
        right: the index of the max_min_track to which the peaks ends;
    """
    pks = peaks.copy()
    pk_split = {k: peaks.loc[peaks[0] == k, :] for k in set(peaks[0])}
    out = []
    for k, v in pk_split.items():
        a1 = searchsorted(max_min_track[k][:, 0], list(v[1]), side="left")
        a2 = searchsorted(max_min_track[k][:, 0], list(v[2]), side="right")
        out.extend([[idx, a, b] for idx, a, b in zip(v.index, a1, a2)])
    out.sort(key=lambda u: u[0])
    out = asarray(out)
    pks["left"] = out[:, 1]
    pks["right"] = out[:, 2]
    return pks


def compare_with_max_not_in_peaks(track, peaks):
    """
    Add a column to the max min track which decoding the percentage of the corresponding signal values in distribution
    of the maximum points not lying in any peaks.
    :param track: the max min tracks which is output of the make_max_min_track function.
    :param peaks: peaks output by the split_max_min_into_peaks function.
    :return: a max min track with a new column added, which are the percentages of the corresponding values in the
        distribution of local maximum points which do not lie in any peaks.
    """
    flags = {k: ones(v.shape[0], dtype=bool) for k, v in track.items()}
    for _, i in peaks.iterrows():
        flags[i[0]][i["left"]:i["right"]] = False
    npdis = {k: track[k][f & (track[k][:, 2] > 0), 1]
             for k, f in flags.items()}
    for v in npdis.values():
        v.sort()
    for k, v in track.items():
        z = asarray(
            [searchsorted(npdis[k], v[:, 1], side="right") / float(len(npdis[k]))]).T
        track[k] = hstack((track[k], z))
    return track


def add_second_diff(diff_1, diff_2, track):
    """
    Add the edge set (zero points of LoG) into the max_min track. The 3rd column of these points are decoded as 2 for
    left edge and -2 for right.
    :param diff_1: The first derivation of the smoothed signal.
    :param diff_2: The second derivation of the smoothed signal.
    :param track: The max_min track where the edges be added.
    :return: A max_min track where the edges added, the points are sorted by the coordinates.
    """
    diff_track = make_max_min_track(diff_1, diff_2)
    for k, v in diff_track.items():
        v[:, 2] *= 2
    track_out = {}
    for k in track.keys():
        t = [list(i) for i in track[k]]
        dt = [list(i) + [0] for i in diff_track[k]]
        t = t + dt
        t.sort(key=lambda u: u[0])
        track_out[k] = asarray(t)
    return track_out
