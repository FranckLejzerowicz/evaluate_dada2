# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
import subprocess
import numpy as np
import pandas as pd
from os.path import isfile
from qiime2 import Artifact, Metadata
from qiime2.plugins.dada2.methods import denoise_single, denoise_paired
from qiime2.plugins.feature_table.methods import filter_samples
from qiime2.plugins.quality_control.visualizers import evaluate_composition


def load_trimmed_seqs(manifest, reverses):
    single, paired = 'Single', ''
    if reverses:
        single, paired = 'Paired', 'PairedEnd'
    trimmed_seqs = Artifact.import_data(
        'SampleData[%sSequencesWithQuality]' % paired,
        manifest, view_type='%sEndFastqManifestPhred33V2' % single)
    return trimmed_seqs


def run_evaluation(ref_q2, sam_q2, depth=1):
    evaluation = evaluate_composition(
        expected_features=ref_q2, observed_features=sam_q2,
        depth=depth, palette='Set1', plot_tar=True, plot_r_value=True,
        plot_r_squared=True, plot_jaccard=True, plot_bray_curtis=True,
        plot_observed_features=True, plot_observed_features_ratio=True)
    return evaluation


def run_denoise(combis, trimmed_seqs, out_files, params):
    for for_rev in combis:
        tab_fp, seq_fp, sta_fp = out_files[tuple(for_rev)]
        if not (isfile(tab_fp) and isfile(seq_fp) and isfile(sta_fp)):
            if len(for_rev) == 2:
                tab, seq, sta = denoise_paired(
                    demultiplexed_seqs=trimmed_seqs,
                    trunc_len_f=for_rev[0],
                    trunc_len_r=for_rev[1],
                    trunc_q=params[0],
                    max_ee_f=params[1],
                    max_ee_r=params[2],
                    n_reads_learn=params[3],
                    chimera_method="consensus",
                    n_threads=1,
                    hashed_feature_ids=True
                )
            else:
                tab, seq, sta = denoise_single(
                    demultiplexed_seqs=trimmed_seqs,
                    trunc_len=for_rev[0],
                    trunc_q=params[0],
                    max_ee=params[1],
                    n_reads_learn=params[3],
                    chimera_method="consensus",
                    n_threads=1,
                    hashed_feature_ids=True
                )
            tab_filt = filter_samples(
                tab,
                min_frequency=1,
                min_features=1,
                filter_empty_features=True
            )
            tab_filt.filtered_table.save(tab_fp)
            seq.save(seq_fp)
            sta.save(sta_fp)


def get_results(out_files):
    dada2 = {}
    for fr, (tab_fp, seq_fp, sta_fp) in out_files.items():
        dada2[fr] = (
            Artifact.load(tab_fp), Artifact.load(seq_fp), Artifact.load(sta_fp)
        )
    return dada2


def get_combis_split(forwards, reverses, n_cores):
    """Make splits as per the number of CPUs"""
    if reverses:
        combis = [it for it in itertools.product(*[forwards, reverses])]
    else:
        combis = [[forward] for forward in forwards]
    combis_split = [x for x in np.array_split(combis, n_cores) if len(x)]
    return combis_split


def spawn_subprocess(cmd):
    # This function will launch the command
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # and retrieve the "standard output"
    out, err = p.communicate()
    out = out.decode().strip().split('\n')
    return out


def get_stats_pd(dada2):
    """Concatenate the stats from the runs"""
    stats_pds = []
    for fr, (tab, seq, sta) in dada2.items():
        stats_pd = sta.view(Metadata).to_dataframe()
        stats_pd['for-rev'] = '-'.join(map(str, fr))
        stats_pd['forward'] = fr[0]
        if len(fr) == 2:
            stats_pd['reverse'] = fr[1]
        else:
            stats_pd['reverse'] = 'None'
        stats_pds.append(stats_pd)
    stats_pd = pd.concat(stats_pds).reset_index()
    return stats_pd

