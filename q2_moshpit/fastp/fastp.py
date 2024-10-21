# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt
from qiime2.plugin import Int, Float, Bool

def run_fastp(input_sequences: CasavaOneEightSingleLanePerSampleDirFmt, trim_front1: int = 0, trim_tail1: int = 0, cut_window_size: int = 4, cut_mean_quality: int = 20, n_base_limit: int = 5, length_required: int = 15, qualified_quality_phred: int = 15, unqualified_percent_limit: float = 40.0, compression: int = 2, thread: int = 3) -> CasavaOneEightSingleLanePerSampleDirFmt:
    output_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = input_sequences.manifest.view(pd.DataFrame)
    for sample_id, row in manifest.iterrows():
        input_fp = row['absolute-filepath']
        output_fp = os.path.join(output_sequences.path, sample_id + '.fastq.gz')
        cmd = [
            'fastp',
            '--in1', input_fp,
            '--out1', output_fp,
            '--trim_front1', str(trim_front1),
            '--trim_tail1', str(trim_tail1),
            '--cut_window_size', str(cut_window_size),
            '--cut_mean_quality', str(cut_mean_quality),
            '--n_base_limit', str(n_base_limit),
            '--length_required', str(length_required),
            '--qualified_quality_phred', str(qualified_quality_phred),
            '--unqualified_percent_limit', str(unqualified_percent_limit),
            '--compression', str(compression),
            '--thread', str(thread)
        ]
        subprocess.run(cmd, check=True)
    return output_sequences
