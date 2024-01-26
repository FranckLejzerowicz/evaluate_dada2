# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd

from qiime2 import Artifact
from evaluate_dada2.io import qzv_unzip
from evaluate_dada2.mock import (
    get_asv_mock_sample, get_tax_mock_sample, get_mock_melt)
from evaluate_dada2.q2 import run_evaluation


def eval_asv(eval_dir, mock_q2s, t, m, p, f, r):
    ref_q2 = mock_q2s[str(p)][0]
    sam = get_asv_mock_sample(t)
    sam_q2 = Artifact.import_data('FeatureTable[RelativeFrequency]', sam.T)
    evaluation = run_evaluation(ref_q2, sam_q2, 1)
    evaluation_fp = '%s/asv/clust-%s_%s_%s-%s' % (eval_dir, p, m, f, r)
    evaluation.visualization.save(evaluation_fp)
    res = qzv_unzip(eval_dir, evaluation_fp)
    return res, sam


def eval_tax(eval_dir, mock_q2s, sam, m, p, f, r, refs, ranks):
    ref_q2 = mock_q2s[str(p)][1]
    tax = get_tax_mock_sample(sam, refs, ranks)
    tax_q2 = Artifact.import_data('FeatureTable[RelativeFrequency]', tax)
    evaluation = run_evaluation(ref_q2, tax_q2, 7)
    evaluation_fp = '%s/taxo/clust-%s_%s_%s-%s' % (eval_dir, p, m, f, r)
    evaluation.visualization.save(evaluation_fp)
    res = qzv_unzip(eval_dir, evaluation_fp)
    return res


def collect(out, res, f, r, p, sam, level):
    vs = {'f': f, 'r': r, 'p': p, 'sam': sam, 'type': level}
    for d, dat in res.items():
        if not dat.shape[0]:
            dat = pd.DataFrame({'Taxon': ['None'], 'mock': [np.nan]})
        for k, v in vs.items():
            dat[k] = v
        out[d].append(dat)


def get_outs(dada2, eval_dir, mocks, hits_pd, mock_q2s, refs, ranks):
    out = {'false_neg': [], 'misclass': [], 'underclass': [], 'results': []}
    for fdx, ((f, r), (tab, _, __)) in enumerate(dada2.items()):
        print('[%s/%s] for: %s - rev: %s' % (fdx + 1, len(dada2), f, r))
        mock_pd = tab.view(pd.DataFrame).T[sorted(mocks)]
        mock_melt = get_mock_melt(mock_pd, hits_pd, list(mocks), f, r)
        for m in mocks:
            gb_col = ['perc_ident', 'ref', m]
            for p, t in mock_melt[gb_col].groupby('perc_ident'):
                res, sam = eval_asv(eval_dir, mock_q2s, t, m, p, f, r)
                collect(out, res, f, r, p, m, "asv")
                res = eval_tax(eval_dir, mock_q2s, sam, m, p, f, r, refs, ranks)
                collect(out, res, f, r, p, m, "taxo")
    outs = {k: pd.concat(v) for k, v in out.items()}
    return outs
