# -*- coding: utf-8 -*-
# @Time    : 18-5-6 下午8:41
# @Author  : Matrix
# @Site    :
# @File    : pre_processing_reads.py
# @Software: PyCharm
from __future__ import print_function
from scipy import *
import h5py
import pysam
from . import aid_scripts
import gc
from numpy import *


class BAMUnsortedError(Exception):
    def __init__(self, file_name):
        super(Exception, self).__init__()
        self.file_name = file_name

    def __str__(self):
        return "SAM/BAM file %s unsorted. " % s


def test_make_singal_track(
    sam_file,
    signal_file,
    chrom_mask=(),
    shift_plus=4,
    shift_minus=-5,
    extend=5,
    chunk_size=10000,
):
    # Open a hdf5 file to store the signal tracts.
    with h5py.File(signal_file, "w") as hdf:
        # Open the sam_file.
        with pysam.AlignmentFile(sam_file) as sam_in:
            if (
                "SO" in sam_in.header["HD"]
                and sam_in.header["HD"]["SO"] == "coordinate"
            ):
                pass
            else:
                raise BAMUnsortedError(sam_file)
            sq = {
                k: sam_in.header.get_reference_length(k)
                for k in sam_in.header.references
                if not k in chrom_mask
            }
            chrom_name = ""
            cov_sig = None
            site_sig = None
            rds = {}
            idx = 0
            for r in sam_in.fetch(until_eof=True):
                if r.reference_name in chrom_mask or not r.is_proper_pair:
                    continue
                if not r.reference_name == chrom_name:
                    if (not (cov_sig is None)) and (not (site_sig is None)):
                        hdf.create_dataset("coverage/raw/%s" % chrom_name, data=cov_sig)
                        hdf.create_dataset("sites/raw/%s" % chrom_name, data=site_sig)
                    del cov_sig
                    del site_sig
                    gc.collect()
                    chrom_name = r.reference_name
                    cov_sig = zeros(sq[r.reference_name], dtype=float64)
                    site_sig = zeros(sq[r.reference_name], dtype=float64)
                rds.setdefault(r.query_name, []).append(r)
                if len(rds[r.query_name]) == 2:
                    r1, r2 = rds[r.query_name]
                    r1, r2 = (r1, r2) if r2.is_reverse else (r2, r1)
                    p1 = r1.reference_start + shift_plus
                    p2 = r2.reference_end + shift_minus
                    cov_sig[max(0, p1) : min(p2 + 1, sq[r.reference_name])] += 1
                    site_sig[
                        max(0, p1 - extend) : min(sq[r.reference_name], p1 + 1 + extend)
                    ] += 1
                    site_sig[
                        max(0, p2 - extend) : min(sq[r.reference_name], p2 + 1 + extend)
                    ] += 1
                    del rds[r.query_name]
                    idx += 1
                    print(idx)
        hdf.create_dataset("coverage/raw/%s" % chrom_name, data=cov_sig)
        hdf.create_dataset("sites/raw/%s" % chrom_name, data=site_sig)
    return
