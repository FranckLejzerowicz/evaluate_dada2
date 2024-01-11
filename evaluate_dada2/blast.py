# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
from qiime2 import Artifact, Metadata
from evaluate_dada2.q2 import spawn_subprocess


def get_ref(params_pd):
    max_bitscore = params_pd['bitscore'].max()
    max_pident = params_pd['pident'].max()
    max_qcovs = params_pd['qcovs'].max()
    # take the hits with best at both params
    bests_pd = params_pd[
        (params_pd['pident'] == max_pident) &
        (params_pd['qcovs'] == max_qcovs)]
    # if more than one perfect hits
    if bests_pd.shape[0] > 1:
        ref = "Multiple_perfect_hits"
        cause = 'Multiple_perfect_hits'
    elif bests_pd.shape[0] == 1:
        ref = list(bests_pd['sseqid'])[0]
        cause = 'One_perfect_hit'
    else:
        if max_pident > 80:
            bests_pd = params_pd[
                (params_pd['pident'] == max_pident) |
                (params_pd['qcovs'] == max_qcovs)]
            best_score = bests_pd[
                bests_pd['bitscore'] == max_bitscore]
            if best_score.shape[0] == 1:
                ref = list(best_score['sseqid'])[0]
                cause = 'Min_80_is_best'
            else:
                ref = "Multiple_hits"
                cause = 'Min_80_is_multiple'
        else:
            ref = "Other"
            cause = 'Less_than_80'
    return ref, cause


def makeblastdb(blast_db):
    cmd = ['makeblastdb', '-dbtype', 'nucl', '-in', blast_db]
    spawn_subprocess(cmd)


def write_seq_to_blast(seq_out, tab, seq, mock_sams):
    mock_tab = tab.view(pd.DataFrame).T[mock_sams]
    mock_seqs_ids = mock_tab[mock_tab[mock_sams].sum(1) > 0].index
    mock_seqs = seq.view(Metadata).to_dataframe().loc[mock_seqs_ids]
    with open(seq_out, "w") as o:
        for r, row in mock_seqs.iterrows():
            o.write('>%s\n%s\n' % (r, row["Sequence"]))
    return mock_seqs_ids


def blastn(seq_out, blast_db):
    out_cols = ['qseqid', 'sseqid', 'pident', 'bitscore', 'qcovs']
    cmd = ['/Users/franck/programs/ncbi-blast-2.12.0+/bin/blastn',
           '-query', seq_out, '-db', blast_db,
           '-outfmt', '6 delim=@ %s' % ' '.join(out_cols)]
    blast_out = pd.DataFrame(
        [x.split('@') for x in spawn_subprocess(cmd)],
        columns=out_cols)
    return blast_out


def run_blasts(eval_dir, results, mock_sams, blast_dbs):
    """Perform the BLAST searches"""
    blast_ins = []
    blast_outs = []
    for (fwd, rev), (tab, seq, _) in results.items():
        seq_out = '%s/%s-%s_toblast.fa' % (eval_dir, fwd, rev)
        mock_seqs_ids = write_seq_to_blast(seq_out, tab, seq, mock_sams)
        blast_ins.append([fwd, rev, len(mock_seqs_ids)])
        for p, blast_db in blast_dbs.items():
            blast_out = blastn(seq_out, blast_db)
            blast_out['forward'] = fwd
            blast_out['reverse'] = rev
            blast_out['perc_identity'] = p
            blast_outs.append(blast_out)
        os.remove(seq_out)

    tab_out = '%s/blast_out.tsv' % eval_dir
    blast_out = pd.concat(blast_outs)
    blast_out.to_csv(tab_out, index=False, sep='\t')
    blast_out = pd.read_table(tab_out)

    tab_in = '%s/blast_in.tsv' % eval_dir
    blast_in = pd.DataFrame(
        blast_ins, columns=['forward', 'reverse', 'nqueries'])
    blast_in.to_csv(tab_in, index=False, sep='\t')
    blast_in = pd.read_table(tab_in)
    return blast_out, blast_in


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
