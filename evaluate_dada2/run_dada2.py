# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import multiprocessing

from evaluate_dada2.io import (
    get_fors_revs, define_dirs, get_metadata, get_fastqs, get_trimmed_seqs)
from evaluate_dada2.q2 import (
    load_trimmed_seqs, get_combis_split, run_denoise, get_results, get_stats_pd)
from evaluate_dada2.mock import (
    get_ref_seqs, get_ref_tax_d, perform_open_ref, get_mock_refs)
from evaluate_dada2.plots import (
    get_txts, make_heatmap_classifs, make_heatmap_stats,
    make_heatmap_outputs, make_heatmap_blast_asv)
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
        n_cores
):
    mini, maxi, step = trim_range
    forwards, reverses = get_fors_revs(mini, maxi, step, trim_lengths,
                                       f_trim_lengths, r_trim_lengths)
    print("Will trim forward reads to", ' nt, '.join(
        map(str, list(forwards))), 'nt')
    print("Will trim reverse reads to", ' nt, '.join(
        map(str, list(reverses))), 'nt')

    print("Metadata file:", metadata)
    print("Mock community files in:", mock_ref_dir)
    print("Mock community taxonomy:", ref_tax_file)
    print("Metadata variables to check:", '; '.join(sorted(meta_cols)))

    trimmed_dir, denoized_dir, eval_dir, pdf_fp, pdf = define_dirs(base_dir)
    meta, mock_sams, mock_sams_rep = get_metadata(metadata)
    fastqs = get_fastqs(meta, trimmed_dir)
    print("Fastq files in:", base_dir, "[%s samples detected]" % len(fastqs))

    manifest = get_trimmed_seqs(fastqs, denoized_dir)
    trimmed_seqs = load_trimmed_seqs(manifest)
    if 'control_type' not in meta_cols:
        meta_cols = ['control_type'] + sorted(meta_cols)

    ref_seqs = get_ref_seqs(mock_ref_dir)
    ref_tax_d = get_ref_tax_d(mock_ref_dir, ref_tax_file)

    # Run the run_dada2.R script through qiime2
    pool = multiprocessing.Pool(n_cores)
    combis_split = get_combis_split(forwards, reverses, n_cores)
    combis_split = [(x, trimmed_seqs, denoized_dir) for x in combis_split]
    results_ = pool.starmap(run_denoise, combis_split)
    pool.close()
    pool.join()

    results = get_results(results_)
    stats_pd = get_stats_pd(results)

    make_heatmap_outputs(meta, stats_pd, pdf)
    perform_open_ref(results, ref_seqs, mock_sams, meta, meta_cols, pdf)
    blast_dbs, mock_ref_q2s = get_mock_refs(ref_seqs, ref_tax_d, ranks)
    blast_out, blast_in = run_blasts(eval_dir, results, mock_sams, blast_dbs)
    make_heatmap_blast_asv(blast_in, pdf)
    hits_pd = get_hits_pd(blast_out)
    outs = get_outs(eval_dir, results, mock_sams_rep, hits_pd,
                    mock_ref_q2s, ref_tax_d, ranks)
    txts = get_txts()
    make_heatmap_classifs(outs, txts, pdf)
    make_heatmap_stats(outs, txts, pdf)
    pdf.close()
    print('--> Written:', pdf_fp)
