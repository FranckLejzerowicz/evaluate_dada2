# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import itertools
import subprocess
import numpy as np
import pandas as pd
from os.path import isdir, isfile
from multiprocessing import cpu_count
from qiime2 import Artifact, Metadata
from qiime2.plugins.dada2.methods import denoise_paired
from qiime2.plugins.quality_control.visualizers import evaluate_composition


def load_trimmed_seqs(manifest):
    trimmed_seqs = Artifact.import_data(
        'SampleData[PairedEndSequencesWithQuality]',
        manifest, view_type='PairedEndFastqManifestPhred33V2')
    return trimmed_seqs


def run_evaluation(mock_ref_q2, mock_sam_q2, depth=1):
    evaluation = evaluate_composition(
        expected_features=mock_ref_q2, observed_features=mock_sam_q2,
        depth=depth, palette='Set1', plot_tar=True, plot_r_value=True,
        plot_r_squared=True, plot_jaccard=True, plot_bray_curtis=True,
        plot_observed_features=True, plot_observed_features_ratio=True)
    return evaluation


def run_denoise(combis, trimmed_seqs, denoized_dir):
    res = {}
    for (forward, reverse) in combis:
        fr_dir = '%s/%s-%s' % (denoized_dir, forward, reverse)
        if not isdir(fr_dir):
            os.makedirs(fr_dir)
        tab_fp = '%s/table.qza' % fr_dir
        seq_fp = '%s/sequences.qza' % fr_dir
        sta_fp = '%s/stats.qza' % fr_dir
        if not (isfile(tab_fp) and isfile(seq_fp) and isfile(sta_fp)):
            tab, seq, sta = denoise_paired(
                demultiplexed_seqs=trimmed_seqs,
                trunc_len_f=forward,
                trunc_len_r=reverse,
                max_ee_f=2,
                max_ee_r=2,
                trunc_q=20,
                chimera_method="consensus",
                n_threads=1,
                hashed_feature_ids=True
            )
            tab.save(tab_fp)
            seq.save(seq_fp)
            sta.save(sta_fp)
        res[(forward, reverse)] = (tab_fp, seq_fp, sta_fp)
    return res


def get_results(results_):
    results = {}
    for result in results_:
        for fr, (tab_fp, seq_fp, sta_fp) in result.items():
            results[fr] = (Artifact.load(tab_fp),
                           Artifact.load(seq_fp),
                           Artifact.load(sta_fp))
    return results


def get_combis_split(forwards, reverses):
    """Make splits as per the number of CPUs"""
    combis = [it for it in itertools.product(*[forwards, reverses])]
    n_splits = int(cpu_count() / 4)
    if n_splits > 4:
        n_splits = 4
    combis_split = [x for x in np.array_split(combis, n_splits) if len(x)][:6]
    return combis_split


def spawn_subprocess(cmd):
    # This function will launch the command
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # and retrie£ve the "standard output"
    out, err = p.communicate()
    out = out.decode().strip().split('\n')
    return out


def get_stats_pd(results):
    """Concatenate the stats from the runs"""
    stats_pds = []
    for (f, r), (tab, seq, sta) in results.items():
        stats_pd = sta.view(Metadata).to_dataframe()
        stats_pd['for-rev'] = '%s-%s' % (f, r)
        stats_pd['forward'] = f
        stats_pd['reverse'] = r
        stats_pds.append(stats_pd)
    stats_pd = pd.concat(stats_pds).reset_index()
    return stats_pd

