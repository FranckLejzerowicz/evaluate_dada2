# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd

from qiime2 import Artifact, Metadata
from evaluate_dada2.io import qzv_unzip
from evaluate_dada2.mock import (
    get_asv_mock_sample, get_tax_mock_sample, get_mock_melt)
from evaluate_dada2.q2 import run_evaluation


def eval_asv(eval_dir, t, sam, sam_name, mock_q2s, p, f, r):
    mock_sam = t.drop(columns='perc_ident')
    mock_sam = get_asv_mock_sample(mock_sam, sam_name)
    mock_sam_q2 = Artifact.import_data(
        'FeatureTable[RelativeFrequency]', mock_sam.T)
    evaluation = run_evaluation(mock_q2s[str(p)][0], mock_sam_q2, 1)
    evaluation_fp = '%s/asv/clust-%s_%s_%s-%s' % (eval_dir, p, sam, f, r)
    evaluation.visualization.save(evaluation_fp)
    res = qzv_unzip(eval_dir, evaluation_fp)
    return res, mock_sam


def eval_taxo(eval_dir, mock_sam_pd, sam, mock_q2s, ref_tax_d, ranks, p, f, r):
    mock_tax = get_tax_mock_sample(mock_sam_pd, ref_tax_d, ranks)
    mock_tax_q2 = Artifact.import_data(
        'FeatureTable[RelativeFrequency]', mock_tax)
    evaluation = run_evaluation(mock_q2s[str(p)][1], mock_tax_q2, 7)
    evaluation_fp = '%s/taxo/clust-%s_%s_%s-%s' % (eval_dir, p, sam, f, r)
    evaluation.visualization.save(evaluation_fp)
    res = qzv_unzip(eval_dir, evaluation_fp)
    return res


def collect(out, res, f, r, p, sam, level):
    vs = {'f': f, 'r': r, 'p': p, 'sam': sam, 'type': level}
    for d, dat_ in res.items():
        if not dat_.shape[0]:
            dat = pd.DataFrame({'Taxon': ['None'], 'mock': [np.nan]})
        else:
            dat = dat_.copy()
        for k, v in vs.items():
            dat[k] = v
        out.setdefault(d, []).append(dat)


def get_outs(
        eval_dir, results, mock_sams_rep, hits_pd, mock_q2s, ref_tax_d, ranks):
    out = {}
    for fdx, ((f, r), (tab, _, __)) in enumerate(results.items()):
        print('[%s/%s]' % (fdx + 1, len(results)))
        mock_pd = tab.view(pd.DataFrame).T[sorted(mock_sams_rep)]
        mock_melt = get_mock_melt(mock_pd, hits_pd, list(mock_sams_rep), f, r)
        for sam, sam_name in mock_sams_rep.items():
            for p, t in mock_melt[['perc_ident', 'ref', sam]].groupby(
                    'perc_ident'):
                res, mock_sam_pd = eval_asv(eval_dir, t, sam, sam_name,
                                            mock_q2s, p, f, r)
                collect(out, res, f, r, p, sam, "asv")
                res = eval_taxo(eval_dir, mock_sam_pd, sam, mock_q2s,
                                ref_tax_d, ranks, p, f, r)
                collect(out, res, f, r, p, sam, "taxo")
    outs = {k: pd.concat(v) for k, v in out.items()}
    return outs
