# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
import multiprocessing
from matplotlib.backends.backend_pdf import PdfPages

from evaluate_dada2.q2 import (
    load_trimmed_seqs, get_combis_split, run_denoise, get_results, get_stats_pd)
from evaluate_dada2.io import (
    get_fors_revs, define_dirs, get_metadata, get_fastqs,
    get_trimmed_seqs, get_out_files, to_do)
from evaluate_dada2.plots import (
    plot_regressions, get_txts, make_heatmap_classifs, make_heatmap_stats,
    make_heatmap_outputs, make_heatmap_blast_asv)
from evaluate_dada2.mock import get_ref_seqs, get_refs, open_ref, get_mock_refs
from evaluate_dada2.blast import run_blasts, get_hits_pd
from evaluate_dada2.eval import get_outs


def run_dada2(
        base_dir,
        metadata,
        mock_ref_dir,
        ref_tax_file,
        meta_cols,
        ranks,
        trim_range,
        trim_lengths,
        f_trim_lengths,
        r_trim_lengths,
        n_cores,
        sample_regressions,
        trunc_q,
        max_er,
        max_er_rev
):
    mini, maxi, step = trim_range
    params = [trunc_q, max_er, max_er_rev]
    forwards, reverses = get_fors_revs(mini, maxi, step, trim_lengths,
                                       f_trim_lengths, r_trim_lengths)
    combis_split = get_combis_split(forwards, reverses, n_cores)
    print("Will trim forward reads to", ' nt, '.join(
        map(str, list(forwards))), 'nt')
    if reverses:
        print("Will trim reverse reads to", ' nt, '.join(
            map(str, list(reverses))), 'nt')
    print("That is %s combinations" % len([y for x in combis_split for y in x]))

    print("Metadata file:", metadata)
    print("Mock community files in:", mock_ref_dir)
    print("Mock community taxonomy:", ref_tax_file)
    print("Metadata variables to check:", '; '.join(sorted(meta_cols)))

    print("Getting output folders")
    trimmed_dir, denoized_dir, eval_dir, pdf_fp = define_dirs(base_dir)
    pdf = PdfPages(pdf_fp)
    out_files = get_out_files(combis_split, denoized_dir)
    lmplot_fp = '%s/lmplot_data.tsv' % eval_dir

    # metadata things
    print("Loading metadata")
    meta, mocks = get_metadata(metadata)
    if 'control_type' not in meta_cols:
        meta_cols = ['control_type'] + sorted(meta_cols)

    # mock things
    print("Loading mock community reference(s)")
    ref_seqs = get_ref_seqs(mock_ref_dir)
    refs = get_refs(mock_ref_dir, ref_tax_file)

    # DADA2 things
    if to_do(out_files):
        print("Running DADA2")
        fastqs = get_fastqs(meta, trimmed_dir)
        print("Fastq files in", base_dir, "[%s samples detected]" % len(fastqs))
        manifest = get_trimmed_seqs(fastqs, denoized_dir, reverses)
        trimmed = load_trimmed_seqs(manifest, reverses)
        pool = multiprocessing.Pool(n_cores)
        pool_params = [(x, trimmed, out_files, params) for x in combis_split]
        pool.starmap(run_denoise, pool_params)
        pool.close()
        pool.join()

    print("Reading DADA2 results")
    dada2 = get_results(out_files)
    stats_pd = get_stats_pd(dada2)

    print("Making heatmaps from DADA2 stat results")
    make_heatmap_outputs(meta, stats_pd, pdf)
    if sample_regressions:
        if not os.path.isfile(lmplot_fp):
            print("Open-reference clustering on the mock references")
            plots_pd = open_ref(dada2, ref_seqs, mocks, meta, meta_cols)
            plots_pd.to_csv(lmplot_fp, index=False, sep='\t')
        else:
            plots_pd = pd.read_table(lmplot_fp)
        print("Making regressions for relative abundances of samples/mock ASVs")
        plot_regressions(plots_pd, pdf)

    blast_in = '%s/blast_in.tsv' % eval_dir
    blast_out = '%s/blast_out.tsv' % eval_dir
    print("Loading reference mock into qiime2 and for BLASTn")
    blast_dbs, mock_q2s = get_mock_refs(ref_seqs, refs, ranks)
    if not (os.path.isfile(blast_in) and os.path.isfile(blast_out)):
        print("Running BLASTn for ASVs vs mock references")
        run_blasts(dada2, eval_dir, mocks, blast_dbs, blast_in, blast_out)
    blast_out_pd = pd.read_table(blast_out)
    blast_in_pd = pd.read_table(blast_in)
    print("Making heatmap of the BLASTed ASVs numbers")
    make_heatmap_blast_asv(blast_in_pd, pdf)
    print("Parsing the BLASTn hits")
    hits_pd = get_hits_pd(blast_out_pd)
    print("Running Qiime2's evaluate-composition for samples' mocks features")
    outs = get_outs(dada2, eval_dir, mocks, hits_pd, mock_q2s, refs, ranks)
    print("Making heatmap from the Qiime2's evaluate-composition results")
    txts = get_txts()
    make_heatmap_classifs(outs, txts, pdf)
    make_heatmap_stats(outs, txts, pdf)
    pdf.close()
    print('--> Written:', pdf_fp)
