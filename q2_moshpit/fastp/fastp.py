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

def run_fastp(
    input_sequences: CasavaOneEightSingleLanePerSampleDirFmt, 
    trim_front1: int = 0, 
    trim_tail1: int = 0, 
    cut_window_size: int = 4, 
    cut_mean_quality: int = 20, 
    n_base_limit: int = 5, 
    length_required: int = 15, 
    qualified_quality_phred: int = 15, 
    unqualified_percent_limit: int = 40, 
    compression: int = 2, 
    thread: int = 3, 
    trim_front2: int = 0, 
    trim_tail2: int = 0, 
    adapter_sequence: str = '', 
    adapter_sequence_r2: str = '', 
    poly_g_min_len: int = 10, 
    poly_x_min_len: int = 10, 
    overlap_len_require: int = 30, 
    overlap_diff_limit: int = 5, 
    overlap_diff_percent_limit: int = 20, 
    correction: bool = False
    ) -> CasavaOneEightSingleLanePerSampleDirFmt:
    output_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = input_sequences.manifest
    for sample_id, row in manifest.iterrows():
        input_fp = row['forward']
        output_fp = os.path.join(output_sequences.path, os.path.basename(row['forward']))
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
            '--thread', str(thread),
            '--adapter_sequence', adapter_sequence,
            '--poly_g_min_len', str(poly_g_min_len),
            '--poly_x_min_len', str(poly_x_min_len),
            '--overlap_len_require', str(overlap_len_require),
            '--overlap_diff_limit', str(overlap_diff_limit),
            '--overlap_diff_percent_limit', str(overlap_diff_percent_limit)
        ]
        if correction:
            cmd.append('--correction')
        if 'reverse' in row and row['reverse'] is not None:
            input_fp2 = row['reverse']
            output_fp2 = os.path.join(output_sequences.path, os.path.basename(row['reverse']))
            cmd.extend([
                '--in2', input_fp2,
                '--out2', output_fp2,
                '--trim_front2', str(trim_front2),
                '--trim_tail2', str(trim_tail2),
                '--adapter_sequence_r2', adapter_sequence_r2
            ])
        subprocess.run(cmd, check=True)
    return output_sequences
