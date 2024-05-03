# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import tempfile

import qiime2
import skbio.io
import skbio.sequence

import pandas as pd
from qiime2 import Artifact

from q2_assembly._utils import run_commands_with_pipe
from q2_moshpit._utils import run_command
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import (
    MultiBowtie2IndexDirFmt, CasavaOneEightSingleLanePerSampleDirFmt, BAMDirFmt
)


def get_mag_length(fasta_file):
    sequences = skbio.io.read(fasta_file, format='fasta', into=skbio.DNA)
    return sum(len(seq) for seq in sequences)


def estimate_abundance(
        ctx,
        reads,
        mags,
        map,
):
    # TODO: implement single- and paired-end case
    maps = map.view(BAMDirFmt)
    r = reads.view(CasavaOneEightSingleLanePerSampleDirFmt)
    m = mags.view(MAGSequencesDirFmt)

    # get read counts per sample
    _tabulate_counts = ctx.get_action("demux", "tabulate_read_counts")
    counts, = _tabulate_counts([reads])
    counts_df = counts.view(qiime2.Metadata).to_dataframe()

    # calculate MAG lengths
    # TODO: replace by FeatureData[SequenceCharacteristics % length]
    lengths = {}
    for mag_id, mag_fp in m.feature_dict().items():
        lengths[mag_id] = get_mag_length(mag_fp)
    lengths_df = pd.DataFrame.from_dict(lengths, orient="index", columns=["length"])

    # get sample IDs from reads and BAMs
    sample_ids = r.manifest.index.to_list()
    sample_ids_bam = {
        os.path.basename(x).split("_alignment")[0]: x for x
        in glob.glob(os.path.join(str(maps), "*.bam"))
    }
    if set(sample_ids) != set(sample_ids_bam.keys()):
        # TODO: add a better error message
        raise ValueError("Sample IDs in reads and BAMs do not match.")

    with tempfile.TemporaryDirectory() as temp_dir:
        report_dfs = []
        for sample_id, sample_fp in sample_ids_bam.items():
            sample_dir = os.path.join(temp_dir, sample_id)
            os.makedirs(sample_dir)

            output_fp = os.path.join(str(temp_dir), f"{sample_id}.bam")
            coverage_fp = os.path.join(str(temp_dir), f"{sample_id}.coverage.tsv")

            cmd1 = [
                "samtools", "sort", "-o", output_fp, sample_fp
            ]
            run_command(cmd1, verbose=True)

            cmd3 = [
                "samtools", "coverage", "-o", coverage_fp, output_fp
            ]
            run_command(cmd3, verbose=True)

            df = pd.read_csv(coverage_fp, sep="\t", index_col=0)
            df["sample_id"] = sample_id
            df["mag_id"] = df.index.map(lambda x: x.split("_", maxsplit=1)[0])
            report_dfs.append(df)

        coverage_df = pd.concat(report_dfs)
        coverage_summed = coverage_df.groupby(["sample_id", "mag_id"]).sum().reset_index(drop=False)
        coverage_summed = coverage_summed.merge(lengths_df, left_on="mag_id", right_index=True)
        coverage_summed = coverage_summed.merge(counts_df, left_on="sample_id", right_index=True)

        coverage_summed["abundance"] = coverage_summed["numreads"] * 10**6 / (coverage_summed["Demultiplexed sequence count"] * coverage_summed["length"])

    ft = coverage_summed.pivot(index='sample_id', columns='mag_id', values='abundance')
    ft.index.name = "sample-id"

    return Artifact.import_data("FeatureTable[Frequency]", ft)
