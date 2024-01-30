# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import shutil
import zipfile
import pandas as pd
from os.path import isdir, isfile


def get_fors_revs(
        mini=150,
        maxi=250,
        step=25,
        values=(),
        f_values=(),
        r_values=()
):
    if f_values and not r_values:
        forwards = [int(x) for x in f_values]
        reverses = []
    elif values:
        forwards = [int(x) for x in values]
        reverses = [int(x) for x in values]
    elif f_values and r_values:
        forwards = [int(x) for x in f_values]
        reverses = [int(x) for x in r_values]
    else:
        forwards = range(mini, (maxi+1), step)
        reverses = range(mini, (maxi+1), step)
    return forwards, reverses


def mk_dirs(to_create):
    for d in to_create:
        if not isdir(d):
            os.makedirs(d)


def define_dirs(base_dir):
    # fastq directory with the samples we are using
    to_create = []
    trimmed_dir = "%s/01_trimmed" % base_dir
    denoized_dir = "%s/02_denoized" % base_dir
    to_create.append(denoized_dir)

    eval_dir = "%s/03_evaluated" % base_dir
    for subdir in ['asv', 'taxo']:
        to_create.append('%s/%s' % (eval_dir, subdir))

    figure_dir = "%s/figures" % base_dir
    to_create.append(figure_dir)
    mk_dirs(to_create)

    pdf_fp = '%s/denoizing_exploration.pdf' % figure_dir
    return trimmed_dir, denoized_dir, eval_dir, pdf_fp


def to_do(out_files):
    for vals in out_files.values():
        s = sum([isfile(x) for x in vals])
        if s == 3:
            continue
        return True


def get_out_files(combis_split, denoized_dir):
    out_files = {}
    for combis in combis_split:
        for for_rev in combis:
            fr_dir = '%s/%s' % (denoized_dir, '-'.join(map(str, for_rev)))
            if not isdir(fr_dir):
                os.makedirs(fr_dir)
            tab_fp = '%s/table.qza' % fr_dir
            seq_fp = '%s/sequences.qza' % fr_dir
            sta_fp = '%s/stats.qza' % fr_dir
            out_files[tuple(for_rev)] = (tab_fp, seq_fp, sta_fp)
    return out_files


def get_metadata(metadata):
    meta = pd.read_table(metadata)
    mock_sams = list(
        meta[meta['control_type'] == 'control positive'].sample_name)
    return meta, mock_sams


def get_fastqs(meta, trimmed_dir):
    fastqs = {}
    for sample_name in meta['sample_name']:
        fastqs[sample_name] = glob.glob(
            '%s/%s_*_R1_*.fastq.gz' % (trimmed_dir, sample_name)
        ) + glob.glob(
            '%s/%s_*_R2_*.fastq.gz' % (trimmed_dir, sample_name))
    return fastqs


def get_trimmed_seqs(fastqs, denoized_dir, reverses):
    manifest = '%s/MANIFEST' % denoized_dir
    if reverses:
        h = 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n'
    else:
        h = 'sample-id\tabsolute-filepath\n'
    with open(manifest, 'w') as o:
        o.write(h)
        for sample_name, rs in sorted(fastqs.items()):
            o.write('%s\t%s\n' % (sample_name, '\t'.join(rs)))
    return manifest


def get_fwd_rev(fr):
    fwd, rev = fr[0], 'None'
    if len(fr) == 2:
        rev = fr[1]
    return fwd, rev


def qzv_unzip(eval_dir, evaluation_fp):
    """
    Adapted from:
    https://forum.qiime2.org/t/how-to-save-the-csv-create-a-table-from-the-barplot-visualisation-using-qiime2-api/17801/5
    """
    exdir = '%s/tmp' % eval_dir
    with zipfile.ZipFile('%s.qzv' % evaluation_fp, 'r') as zip_ref:
        zip_ref.extractall(exdir)
    inf = exdir + '/' + os.listdir(exdir)[-1] + "/data"
    false_neg_fp = inf + '/false_negative_features.tsv'
    misclass_fp = inf + '/misclassifications.tsv'
    underclass_tsv = inf + '/underclassifications.tsv'
    results_fp = inf + '/results.tsv'
    false_neg = pd.read_table(false_neg_fp)
    misclass = pd.read_table(misclass_fp)
    underclass = pd.read_table(underclass_tsv)
    results = pd.read_table(results_fp)
    qzv_outs = {'false_neg': false_neg, 'misclass': misclass,
                'underclass': underclass, 'results': results}
    shutil.rmtree(exdir)
    return qzv_outs
