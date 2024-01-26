# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from evaluate_dada2 import __version__
from evaluate_dada2.run_dada2 import run_dada2

@click.command()
@click.option(
    "-i", "--i-fastq-dir", default=None, nargs=1,
    help="Folder containing the fastq files")
@click.option(
    "-m", "--i-metadata", default=None, nargs=1,
    help="Metadata file (tab-separated) column `sample_name` match fastq files")
@click.option(
    "-mi", "--i-mock-dir", default=None, nargs=1,
    help="Folder containing the mocks sequnces and clusters")
@click.option(
    "-mt", "--i-mock-tax", default=None, nargs=1,
    help="Name of the taxonomy .tsv file in the `--i-mock-dir` folder")
@click.option(
    "-r", "--p-ranks", default=['d', 'p', 'c', 'o', 'f', 'g', 's'],
    multiple=True, type=str, help="Taxonomy ranks for the mock evaluation")
@click.option(
    "-c", "--p-meta-vars",  multiple=True, type=str, show_default=False,
    default=[], help="Metadata variables to use for checking samples vs mock")
@click.option(
    "-r", "--p-trim-range", nargs=3, show_default=False,
    default=(150, 250, 25), help="Min, Max, Step for the trimming length")
@click.option(
    "-l", "--p-trim-lengths", multiple=True, type=int, show_default=False,
    default=[], help="Trimming lengths")
@click.option(
    "-lf", "--p-f-trim-lengths", multiple=True, type=int, show_default=False,
    default=[], help="Trimming lengths (forward reads)")
@click.option(
    "-lr", "--p-r-trim-lengths", multiple=True, type=int, show_default=False,
    default=[], help="Trimming lengths (reverse reads)")
@click.option(
    "-n", "--p-n-cores", type=int, nargs=1, show_default=False,
    default=4, help="Number of cores for multiprocessing")
@click.option(
    "--sample-regressions/--no-sample-regressions",
    default=False, show_default=True,
    help="Make regression for mock features in actual samples")
@click.version_option(__version__, prog_name="evaluate_dada2")


def standalone_dada2(
        i_fastq_dir,
        i_metadata,
        i_mock_dir,
        i_mock_tax,
        p_meta_vars,
        p_ranks,
        p_trim_range,
        p_trim_lengths,
        p_f_trim_lengths,
        p_r_trim_lengths,
        p_n_cores,
        sample_regressions
):

    run_dada2(
        base_dir=i_fastq_dir,
        metadata=i_metadata,
        mock_ref_dir=i_mock_dir,
        ref_tax_file=i_mock_tax,
        meta_cols=p_meta_vars,
        ranks=p_ranks,
        trim_range=p_trim_range,
        trim_lengths=p_trim_lengths,
        f_trim_lengths=p_f_trim_lengths,
        r_trim_lengths=p_r_trim_lengths,
        n_cores=p_n_cores,
        sample_regressions=sample_regressions
    )


if __name__ == "__main__":
    standalone_dada2()
