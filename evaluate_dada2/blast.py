# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
from qiime2 import Metadata
from evaluate_dada2.q2 import spawn_subprocess
from evaluate_dada2.io import get_fwd_rev


def get_bests_pd(params_pd, max_pident, max_qcovs, or_=0):
    # take the hits with best at both params, or either
    if or_:
        bests_pd = params_pd[
            (params_pd['pident'] == max_pident) |
            (params_pd['qcovs'] == max_qcovs)]
    else:
        bests_pd = params_pd[
            (params_pd['pident'] == max_pident) &
            (params_pd['qcovs'] == max_qcovs)]
    return bests_pd


def get_ref_cause_over80(params_pd, max_pident, max_qcovs, max_bitscore):
    bests_pd = get_bests_pd(params_pd, max_pident, max_qcovs, 1)
    best_score = bests_pd[bests_pd['bitscore'] == max_bitscore]
    if best_score.shape[0] == 1:
        ref = list(best_score['sseqid'])[0]
        cause = 'Min_80_is_best'
    elif best_score.shape[0]:
        ref = "Multiple_hits"
        cause = 'Min_80_is_multiple'
    else:
        ref = "Best_HSP_score_ambiguous"
        cause = 'Min_80_is_multiple'
    return ref, cause


def get_ref(params_pd):
    max_bitscore = params_pd['bitscore'].max()
    max_pident = params_pd['pident'].max()
    max_qcovs = params_pd['qcovs'].max()
    bests_pd = get_bests_pd(params_pd, max_pident, max_qcovs)
    # if more than one perfect hits
    if bests_pd.shape[0] > 1:
        ref = "Multiple_perfect_hits"
        cause = 'Multiple_perfect_hits'
    elif bests_pd.shape[0] == 1:
        ref = list(bests_pd['sseqid'])[0]
        cause = 'One_perfect_hit'
    else:
        if max_pident > 80:
            ref, cause = get_ref_cause_over80(
                params_pd, max_pident, max_qcovs, max_bitscore)
        else:
            ref = "Other"
            cause = 'Less_than_80'
    return ref, cause


def makeblastdb(blast_db):
    cmd = ['makeblastdb', '-dbtype', 'nucl', '-in', blast_db]
    spawn_subprocess(cmd)


def write_seq_to_blast(seq_out, mock_tab, seq, mocks):
    mock_tab = mock_tab[mocks]
    mock_seqs_ids = mock_tab[mock_tab[mocks].sum(1) > 0].index
    mock_seqs = seq.view(Metadata).to_dataframe().loc[mock_seqs_ids]
    with open(seq_out, "w") as o:
        for r, row in mock_seqs.iterrows():
            o.write('>%s\n%s\n' % (r, row["Sequence"]))
    return mock_seqs_ids


def blastn(seq_out, blast_db):
    out_cols = ['qseqid', 'sseqid', 'pident', 'bitscore', 'qcovs']
    cmd = ['blastn', '-query', seq_out, '-db', blast_db,
           '-outfmt', '6 delim=@ %s' % ' '.join(out_cols)]
    blast_out = pd.DataFrame(
        [x.split('@') for x in spawn_subprocess(cmd)],
        columns=out_cols)
    return blast_out


def run_blasts(dada2, eval_dir, mocks, blast_dbs, blast_in, blast_out):
    """Perform the BLAST searches"""
    blast_ins_pds = []
    blast_outs_pds = []
    for fr, (tab, seq, _) in dada2.items():
        fwd, rev = get_fwd_rev(fr)
        mock_tab = tab.view(pd.DataFrame).T
        if set(mocks).difference(set(mock_tab.columns)):
            continue
        seq_out = '%s/%s_toblast.fa' % (eval_dir, '-'.join(map(str, fr)))
        mock_seqs_ids = write_seq_to_blast(seq_out, mock_tab, seq, mocks)
        blast_ins_pds.append([fwd, rev, len(mock_seqs_ids)])
        for p, blast_db in blast_dbs.items():
            blast_out_pd = blastn(seq_out, blast_db)
            blast_out_pd['forward'] = fwd
            blast_out_pd['reverse'] = rev
            blast_out_pd['perc_identity'] = p
            blast_outs_pds.append(blast_out_pd)
        os.remove(seq_out)

    blast_outs_pd = pd.concat(blast_outs_pds)
    blast_outs_pd.to_csv(blast_out, index=False, sep='\t')

    blast_ins_pd = pd.DataFrame(
        blast_ins_pds, columns=['forward', 'reverse', 'nqueries'])
    blast_ins_pd.to_csv(blast_in, index=False, sep='\t')


def get_hits_pd(blast_out):
    """Parse blast results"""
    hits = []
    gb_cols = ['forward', 'reverse', 'perc_identity', 'qseqid']
    for pdx, ((f, r, p, q), params_pd) in enumerate(blast_out.groupby(gb_cols)):
        ref, cause = get_ref(params_pd)
        hits.append([f, r, p, q, ref, cause])
    hits_pd = pd.DataFrame(hits, columns=[
        'forward', 'reverse', 'perc_identity', 'seq', 'ref', 'cause'])
    return hits_pd
