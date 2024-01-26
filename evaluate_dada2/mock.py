# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import itertools

import pandas as pd
from os.path import dirname
from qiime2 import Artifact
from qiime2.plugins.vsearch.pipelines import cluster_features_open_reference
from qiime2.plugins.feature_table.methods import relative_frequency
from evaluate_dada2.blast import makeblastdb


def get_plot_pd(meta_combis, meta, mock_sam, tab_mock, value_name):
    plot_pds = []
    for cdx, meta_combi in enumerate(meta_combis):
        sams_d = meta.groupby('sample_name').apply(lambda x: '_'.join(
            map(str, x[list(meta_combi)].values[0]))).to_dict()
        sams_d[mock_sam] = 'mock (% reads)'
        combi_pd = tab_mock.copy()
        combi_pd = combi_pd.loc[:, combi_pd.sum() > 0]
        combi_pd = combi_pd / combi_pd.sum()
        combi_pd.columns = [sams_d.get(x, x) for x in combi_pd.columns]
        combi_ml = combi_pd.melt(
            id_vars=['mock (% reads)'],
            value_name=value_name,
            ignore_index=False
        ).reset_index()
        combi_ml['comparison'] = ' & '.join(meta_combi)
        plot_pds.append(combi_ml)
    plot_pd = pd.concat(plot_pds)
    return plot_pd


def get_clusters(ref_seqs, seq, tab):
    clusters = {}
    for p, (_, __, ref_seq) in ref_seqs.items():
        print('Clustering vs DB version p="%s"' % p)
        open_table, open_seqs, _ = cluster_features_open_reference(
            sequences=seq, table=tab, reference_sequences=ref_seq,
            perc_identity=float(p), threads=1)
        print('Turning to relative frequencies')
        clusters[p] = relative_frequency(open_table)
    return clusters


def get_mock_refs(ref_seqs, refs, ranks):
    """Make BLAST databases from the mock references"""
    blast_dbs = {}
    mock_q2s = {}
    for p, (ref_seq_fp, ref_table_fp, ref_seq) in ref_seqs.items():
        makeblastdb(ref_seq_fp)
        blast_dbs[p] = ref_seq_fp
        mock_q2s[p] = get_db_q2(ref_table_fp, refs, ranks)
    return blast_dbs, mock_q2s


def get_tab_mock(tab_clust, mock_sam, mock_sams):
    tab_mock = tab_clust.loc[tab_clust[mock_sam] > 0]
    other_sams = [x for x in mock_sams if x != mock_sam]
    if other_sams:
        tab_mock = tab_mock.drop(columns=other_sams)
    return tab_mock


def get_lmplots(plots_pds, relab, mocks, meta, meta_cols, f, r, p):
    val = 'sample (% reads)'
    tab_clust = relab.relative_frequency_table.view(pd.DataFrame).T
    meta_combis = [it for n in range(len(meta_cols))
                   for it in itertools.combinations(meta_cols, n + 1)]
    for mock in mocks:
        tab_mock = get_tab_mock(tab_clust, mock, mocks)
        empties = 100 * ((tab_mock.sum() == 0).sum() / tab_mock.shape[1])
        plots_pd = get_plot_pd(meta_combis, meta, mock, tab_mock, val)
        plots_pd['mock_sample'] = mock
        plots_pd['forward'] = f
        plots_pd['reverse'] = r
        plots_pd['perc_identity'] = p
        plots_pd['perc_empty_samples'] = empties
        plots_pds.append(plots_pd)


def open_ref(dada2, ref_seqs, mocks, meta, meta_cols):
    """
    Perform open-reference clustering of the mock sample ASVs onto the reference mock sequences
    For the three different `perc_identity` at which the reference mock sequences do cluster
    """
    plots_pds = []
    for (f, r), (tab, seq, _) in dada2.items():
        clusters = get_clusters(ref_seqs, seq, tab)
        for p, relab in clusters.items():
            get_lmplots(plots_pds, relab, mocks, meta, meta_cols, f, r, p)
    plots_pd = pd.concat(plots_pds)
    return plots_pd


def get_db_q2(table_fp, refs, ranks):
    tab = pd.read_table(table_fp)
    asv_q2 = Artifact.import_data(
        'FeatureTable[RelativeFrequency]',
        tab.set_index('featureid').T)
    tab['featureid'] = [refs.get(x, 'd__Eukaryota') for x in tab['featureid']]
    tab['featureid'] = ['%s; %s' % (x, '; '.join([
        '%s__' % r for r in ranks[ranks.index(x.split('; ')[-1][0]) + 1:]]
    )) for x in tab['featureid']]
    tab = tab.groupby('featureid').sum()
    tax_q2 = Artifact.import_data('FeatureTable[RelativeFrequency]', tab.T)
    return asv_q2, tax_q2


def get_mock_melt(mock_sams_pd, hits_pd, sams, f, r):
    fr_pd = hits_pd[
        (hits_pd['forward'] == f) &
        (hits_pd['reverse'] == r)]
    fr_us = fr_pd.drop(
        columns=['forward', 'reverse', 'cause']
    ).set_index(
        ['perc_identity', 'seq']
    ).unstack().T
    fr_us.index = fr_us.index.droplevel()
    fr_mock_pd = pd.concat([fr_us, mock_sams_pd], axis=1)
    fr_mock_pd = fr_mock_pd.apply(lambda x: x.fillna(x.index.to_series()))
    mock_melt = fr_mock_pd.melt(
        id_vars=sams, value_name='ref', var_name='perc_ident'
    ).groupby(
        ['perc_ident', 'ref']
    ).sum().reset_index()
    return mock_melt


def get_asv_mock_sample(t):
    sam = t.drop(columns='perc_ident')
    sam = sam.set_index('ref')
    sam.index.name = 'featureid'
    sam.columns = ['mock']
    sam = sam / sam.sum()
    return sam


def get_tax_mock_sample(sam, refs, ranks):
    sam.index = [refs.get(x, 'd__') for x in sam.index]
    sam.index = ['%s; %s' % (x, '; '.join(
        ['%s__' % r for r in ranks[ranks.index(x.split('; ')[-1][0]) + 1:]]
    )) for x in sam.index]
    tax = sam.groupby(level=0).sum().T
    return tax


def get_ref_seqs(mock_ref_dir):
    """Get the reference mock community sequences and taxonomy"""
    ref_seqs = {}
    ref_clust_fps = glob.glob('%s/clustering/*/sequences.fasta' % mock_ref_dir)
    for ref_clust_fp in ref_clust_fps:
        ref_table_fp = "%s/relative_abundances.tsv" % dirname(ref_clust_fp)
        p = ref_clust_fp.split('/')[-2]
        ref_seqs[p] = (
            ref_clust_fp, ref_table_fp,
            Artifact.import_data('FeatureData[Sequence]', ref_clust_fp))
    return ref_seqs


def get_refs(mock_ref_dir, ref_tax_file):
    ref_tax_fp = '%s/%s' % (mock_ref_dir, ref_tax_file)
    ref_tax = pd.read_table(ref_tax_fp, index_col=0)
    refs = ref_tax.to_dict()['Taxon']
    return refs
