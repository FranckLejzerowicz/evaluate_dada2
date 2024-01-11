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
from os.path import isdir
from matplotlib.backends.backend_pdf import PdfPages
from evaluate_dada2.mock import get_mock_sams_rep


def get_fors_revs(
        mini=150,
        maxi=250,
        step=25,
        values=(),
        f_values=(),
        r_values=()
):
    if values:
        forwards = [int(x) for x in values]
        reverses = [int(x) for x in values]
    elif f_values and r_values:
        forwards = [int(x) for x in f_values]
        reverses = [int(x) for x in r_values]
    else:
        forwards = range(mini, maxi, step)
        reverses = range(mini, maxi, step)
    return forwards, reverses


def define_dirs(base_dir):
    # fastq directory with the samples we are using
    trimmed_dir = "%s/01_trimmed" % base_dir
    denoized_dir = "%s/02_denoized" % base_dir
    eval_dir = "%s/03_evaluated" % base_dir
    for subdir in ['asv', 'taxo']:
        if not isdir('%s/%s' % (eval_dir, subdir)):
            os.makedirs('%s/%s' % (eval_dir, subdir))

    figure_dir = "%s/figure_dir" % base_dir
    if not isdir(figure_dir):
        os.makedirs(figure_dir)
    pdf_fp = '%s/denoizing_exploration.pdf' % figure_dir
    pdf = PdfPages(pdf_fp)
    return trimmed_dir, denoized_dir, eval_dir, pdf_fp, pdf


def get_metadata(metadata):
    meta = pd.read_table(metadata)
    mock_sams = list(
        meta[meta['control_type'] == 'control positive'].sample_name)
    mock_sams_rep = get_mock_sams_rep(mock_sams)
    return meta, mock_sams, mock_sams_rep


def get_fastqs(meta, trimmed_dir):
    fastqs = {}
    for sample_name in meta['sample_name']:
        fastqs[sample_name] = glob.glob(
            '%s/%s*_R1_*.fastq.gz' % (trimmed_dir, sample_name)
        ) + glob.glob(
            '%s/%s*_R2_*.fastq.gz' % (trimmed_dir, sample_name))
    return fastqs


def get_trimmed_seqs(fastqs, denoized_dir):
    manifest = '%s/MANIFEST' % denoized_dir
    with open(manifest, 'w') as o:
        o.write(
            'sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n')
        for sample_name, (r1, r2) in sorted(fastqs.items()):
            o.write('%s\t%s\t%s\n' % (sample_name, r1, r2))
    return manifest


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
    qza_outs = {'false_neg': pd.read_table(false_neg_fp),
                'misclass': pd.read_table(misclass_fp),
                'underclass': pd.read_table(underclass_tsv),
                'res': pd.read_table(results_fp)}
    shutil.rmtree(exdir)
    return qza_outs
