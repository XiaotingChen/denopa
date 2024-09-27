# -*- coding: utf-8 -*-
# @Time    : 18-5-6 下午9:53
# @Author  : Matrix
# @Site    :
# @File    : __init__.py
# @Software: PyCharm
from __future__ import print_function
from . import pileup_signals
from . import smoothed_signal
from . import call_peak  # , aid_scripts
from . import signal_track_builder
from . import candidate_mm_process
from . import ocsvm_model
from . import determine_dynamic
from . import determineNFR
from . import dbscan_model
from . import fragmentLengthsDist
import os
import sys
import h5py
import pickle as pk
import multiprocessing as mp
from numpy import *
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import copy

def makeSignalTracks(
    filesIn, folderOut, pname, bufferSize=1000000, chromSkip="", chromInculde="",nuc_number=0,fl_est_only=False,fragLenCutOff=0
):
    lock = mp.Lock()
    filesIn = [os.path.abspath(i) for i in filesIn]
    pwd = os.getcwd()
    os.chdir(folderOut)
    fl = pileup_signals.build_signal_track(
        filesIn,
        "%s_pileup_signal" % pname,
        chrom_skip=chromSkip,
        chrom_inculde=chromInculde
    )
    with open("%s_frag_len.pkl" % pname, "wb") as fout:
        pk.dump(fl, fout)
    # fl plot
    fl_df = pd.DataFrame.from_dict(fl, orient='index')
    fl_df.sort_index(inplace=True)
    fl_df.reset_index(inplace=True)
    fl_df.rename(axis=1, mapper={'index': 'fl', 0: 'count'}, inplace=True)
    plt.subplots(figsize=(15, 5))
    sns.lineplot(fl_df, x='fl', y='count')
    plt.savefig('fl.pdf', bbox_inches='tight')
    # remove fragment count below cutoff
    fl_threshold=copy.copy(fl)
    for l in range(fragLenCutOff):
        _=fl_threshold.pop(l,None)
    # fl model
    fl = fragmentLengthsDist.fragmentLengthModel(fl_threshold,nuc_number)
    # export model
    if nuc_number!=0:
        with open('nuc_{}_model'.format(nuc_number), 'wb') as f:
            pickle.dump(fl, f)
    else:
        with open('best_AIC_model', 'wb') as f:
            pickle.dump(fl, f)
    with open('deNOPA_param','wb') as f:
        pickle.dump(fragmentLengthsDist.fix_parameters(fl.params),f)
    if fl_est_only==False:
        fl.nucFreeTrack(
            filesIn,
            "%s_pileup_signal.hdf" % pname,
            "%s_smooth.hdf" % pname,
            chrom_inculde=chromInculde,
            chromSkip=chromSkip,
        )
        proc1 = signal_track_builder.GaussConvolve(
            "%s_pileup_signal.hdf" % pname,
            "%s_smooth.hdf" % pname,
            "coverage",
            72,
            lock=lock,
        )
        proc2 = signal_track_builder.GaussConvolve(
            "%s_pileup_signal.hdf" % pname, "%s_smooth.hdf" % pname, "sites", 24, lock=lock
        )
        # proc3 = signal_track_builder.GaussConvolve("pileup_signal.hdf", "smooth.hdf", "short", 72)
        proc1.start()
        proc2.start()
        proc1.join()
        proc2.join()
        # proc3.start()
        # proc3.join()
        os.chdir(pwd)
    return


def candidateNucleosomes(
    samFiles,
    pname,
    outputFolder,
    maxLen=2000,
    leftShift=+4,
    rightShift=-5,
    proc=1,
    candPvalue=0.1,
    nfrQvalue=0.1,
    arerFile="Candidate_peaks_0.1.bed",
):
    pileupFile = "%s_pileup_signal.hdf" % pname
    smoothFile = "%s_smooth.hdf" % pname
    pwd = os.getcwd()
    samFiles = [os.path.abspath(i) for i in samFiles]
    os.chdir(outputFolder)
    nocFile = os.path.abspath(smoothFile)
    pileupFile = os.path.abspath(pileupFile)
    smoothFile = os.path.abspath(smoothFile)
    with h5py.File(pileupFile, "r") as raw, h5py.File(smoothFile, "r") as smooth:
        df = signal_track_builder.MakeMaxMinTrack(
            smooth["coverage/0"], smooth["coverage/1"], smooth["coverage/2"]
        )()
        peaks_denovo = call_peak.call_candidate_regions(
            df, candPvalue, merge_dist=1000, proc=1
        )
        print(peaks_denovo.shape[0])
        peaks_denovo.to_csv(arerFile, header=None, sep="\t", index=None)
        sites_max_min_track = signal_track_builder.make_max_min_track(
            smooth["sites/0"], smooth["sites/1"]
        )
        pks_mm = signal_track_builder.split_max_min_into_peaks(
            sites_max_min_track, peaks_denovo
        )
        sites_max_min_track = signal_track_builder.compare_with_max_not_in_peaks(
            sites_max_min_track, pks_mm
        )
        sites_max_min_track = signal_track_builder.add_second_diff(
            smooth["sites/1"], smooth["sites/2"], sites_max_min_track
        )
        pks_edge = signal_track_builder.split_max_min_into_peaks(
            sites_max_min_track, peaks_denovo
        )
        cand_mm = candidate_mm_process.filter_mm_candidates(
            sites_max_min_track, pks_edge, 0.05, min_sep=15, max_sep=50
        )
        cand_mg = candidate_mm_process.merge_candidate_mms(
            cand_mm, peaks_denovo[0], sites_max_min_track, min_sep=100, max_sep=215
        )
    num_reads = ocsvm_model.calc_ov_frags(
        samFiles,
        cand_mg,
        peaks_denovo,
        maxLen,
        left_shift=leftShift,
        right_shift=rightShift,
        smooth_file=smoothFile,
        proc=proc,
    )
    with open("%s_candidates.pkl" % pname, "wb") as fout:
        pk.dump(num_reads, fout)
    x = dbscan_model.FinalModel()(num_reads)
    x.to_csv("%s_nucleosomes.txt" % pname, header=None, sep="\t", index=None)
    nfr = determineNFR.NFRDetection(num_reads, smoothFile, nfrQvalue)()
    nfr.to_csv("%s_NFR.txt" % pname, header=None, sep="\t", index=None)
    os.chdir(pwd)
    return


def calcNUC(folderIn):
    makeSignalTracks([os.path.join(folderIn, "ddup.bam")], folderIn)
    pwd = os.getcwd()
    os.chdir(folderIn)
    candidateNucleosomes(["ddup.bam"], "pileup_signal.hdf", "smooth.hdf", ".", proc=16)
    os.remove("pileup_signal.hdf")
    os.remove("smooth.hdf")
    os.chdir(pwd)


def shell():
    pwd = os.getcwd()
    y = []
    f = ["0.025", "0.500", "0.250", "0.750", "0.050", "0.075", "0.100"]
    for i in f:
        for j in range(1, 10):
            try:
                calcNUC("%s/%d" % (i, j))
                print("%s/%d" % (i, j))
            except:
                y.append((i, j))
            finally:
                os.chdir(pwd)
