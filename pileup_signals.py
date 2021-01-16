# -*- coding: utf-8 -*-
# @Time    : 18-5-22 下午10:38
# @Author  : Matrix
# @Site    :
# @File    : pileup_signals.py
# @Software: PyCharm

import pysam
import numpy as np
import h5py
import gc
import logging

logging.basicConfig(format="%(asctime)s: %(message)s", level=logging.INFO)


def build_signal_track(sam_file_list,
                       out_put_prefix,
                       chrom_skip=None,
                       buffer_size=1000000,
                       left_shift=+4,
                       right_shift=-5,
                       extend=10):
    # Open file for out put.
    logger = logging.getLogger(__name__)
    with h5py.File(out_put_prefix + ".hdf", 'w') as fout:
        # Check if all SAM files are coordinately sorted and load the chromosome lengths.
        chrom_size = {}
        for sam_file in sam_file_list:
            with pysam.AlignmentFile(sam_file) as sam:
                assert "SO" in sam.header.as_dict(
                )["HD"] and sam.header.as_dict()["HD"][
                    "SO"] == "coordinate", "Input file %s is not coordinately sorted!" % sam_file
                for k, v in zip(sam.references, sam.lengths):
                    chrom_size.setdefault(k, []).append(v)
        chrom_size = {k: list(set(v)) for k, v in chrom_size.items()}
        assert all(
            [len(v) == 1 for v in chrom_size.values()]
        ), "The chromosome lengths decoded in alignment files are not all the samm. "
        try:
            chrom_size = {
                k: v[0]
                for k, v in chrom_size.items() if not k in chrom_skip
            }
        except TypeError:
            chrom_size = {k: v[0] for k, v in chrom_size.items()}
        # Create groups for signal tracks.
        # Group for coverage of fragments.
        fout.create_group("coverage")
        # Group for coverage of TN5 enzyme sites (5' and 3' most of a fragment, each end is extended both sided by extend).
        fout.create_group("sites")
        # Create data sets
        for k, v in chrom_size.items():
            fout["coverage"].create_dataset(k, shape=(v, ))
            fout["sites"].create_dataset(k, shape=(v, ))
        # Begin loop sam files.
        # Denoting the fragment length distribution.
        fragLenDist = {}
        # Make a counter to denote the total number of reads processed.
        idx = 0
        for sam_file in sam_file_list:
            # Open the sam file.
            with pysam.AlignmentFile(sam_file) as sam:
                # Scan the files and process the reads.
                for r in sam.fetch(until_eof=True):
                    idx += 1
                    if idx % 10000 == 0:
                        logger.info("%d reads has been processed. " % idx)
                    # Discard reads not in interested chromosomes.
                    if not r.reference_name in chrom_size:
                        continue
                    try:
                        # At begin, the variable current_chrom is not defined, making a forced initiation step.
                        if not current_chrom == r.reference_name:
                            # If a new chromosome is reached, all the reads to the current chromosome have been
                            # processed. Write the results to file and re-initiation for the next chromosome.
                            fout["coverage/%s" % current_chrom][
                                buffer_start:buffer_stop] += buffer_cov
                            fout["sites/%s" % current_chrom][
                                buffer_start:buffer_stop] += buffer_ste
                            raise NameError
                    except NameError:
                        # Initiate for a new chromosome.
                        current_chrom = r.reference_name
                        # A buffer (buffer_cov or buffer_ste represent a chromosome location from buffer_start
                        # to buffer_stop).
                        buffer_start = 0
                        # Ensure the buffer is not longer than the chromosome it buffers.
                        buffer_len = min(chrom_size[current_chrom],
                                         buffer_size)
                        buffer_stop = buffer_start + buffer_len
                        # Buffer for coverage track.
                        buffer_cov = np.zeros(buffer_len)
                        # Buffer for sites track.
                        buffer_ste = np.zeros(buffer_len)
                        # Define a dict for read haven't been paired.
                        reads = {}
                    # So far, all the initiations (if needed) has been finished. Begin to process the reads.
                    if r.query_name in reads:
                        # A read pair appears.
                        r1 = reads.pop(r.query_name)
                        # Let the read mapped up stream left and down stream right.
                        rl, rr = (
                            r1, r
                        ) if r1.reference_start <= r.reference_start else (r,
                                                                           r1)
                        # add to the fragment length distribution.
                        fragLen = (rr.reference_end + right_shift) - (
                            rl.reference_start + left_shift)
                        fragLenDist[fragLen] = fragLenDist.setdefault(
                            fragLen, 0) + 1
                        # Left most position.
                        pc1 = rl.reference_start + left_shift
                        # Right most position, note that reference_end is not included in the aligned reads.
                        pc2 = rr.reference_end + right_shift - 1
                        if pc2 + extend + 1 >= buffer_stop:
                            # The reads exceeds the buffered region. Write the current region to file and make a new
                            # one for next region.
                            # Write the current region to file.
                            fout["coverage/%s" % current_chrom][
                                buffer_start:buffer_stop] += buffer_cov
                            fout["sites/%s" % current_chrom][
                                buffer_start:buffer_stop] += buffer_ste
                            # Start position of the new region. Notice that is equals the left most position of reads at
                            # hand minus extend.
                            buffer_start = min(
                                reads.values(),
                                key=lambda u: u.reference_start
                            ).reference_start + left_shift if reads else pc1
                            buffer_start = min(buffer_start, pc1) - extend
                            # Ensure the buffer is not longer than the region it buffers.
                            buffer_len = min(
                                chrom_size[current_chrom] - buffer_start,
                                buffer_size)
                            buffer_stop = buffer_start + buffer_len
                            # Create the buffers.
                            buffer_cov = np.zeros(buffer_len)
                            buffer_ste = np.zeros(buffer_len)
                            # Collect the garbage.
                            gc.collect()
                        # Now the buffers can properly buffers the region the reads represent.
                        # Add to region that the reads covered.
                        # Note pc2 is inclusive.
                        buffer_cov[(pc1 - buffer_start):(pc2 + 1 -
                                                         buffer_start)] += 1
                        try:
                            # Add to the left enzyme site.
                            buffer_ste[(pc1 - extend -
                                        buffer_start):(pc1 + extend -
                                                       buffer_start + 1)] += 1
                            # Add to the right enzyme site.
                            buffer_ste[(pc2 - extend -
                                        buffer_start):(pc2 + extend -
                                                       buffer_start + 1)] += 1
                        except:
                            buffer_ste[max(0, pc1 - extend):(pc1 + extend +
                                                             1)] += 1
                            buffer_ste[(
                                pc2 -
                                extend):min(chrom_size[current_chrom], pc2 +
                                            extend + 1)] += 1
                    else:
                        reads[r.query_name] = r
                # Write the final region to file.
                fout["coverage/%s" %
                     current_chrom][buffer_start:buffer_stop] += buffer_cov
                fout["sites/%s" %
                     current_chrom][buffer_start:buffer_stop] += buffer_ste
                # When a file finished. variable current_chrom is deleted to guarantee the next loop could be
                # initialized.
                del current_chrom
    return fragLenDist
