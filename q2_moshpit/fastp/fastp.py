# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile

from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt
from qiime2.plugin import Int, Float, Bool
from q2_moshpit._utils import run_command
from q2_moshpit.fastp.aggregate import aggregate_fastp_reports


def fastp(
        sequences: CasavaOneEightSingleLanePerSampleDirFmt,
        trim_front1: int = 0,
        trim_tail1: int = 0,
        max_len1: int = 0,
        trim_front2: int = 0,
        trim_tail2: int = 0,
        max_len2: int = 0,
        disable_quality_filtering: bool = False,
        n_base_limit: int = 5,
        qualified_quality_phred: int = 15,
        unqualified_percent_limit: int = 40,
        length_required: int = 15,
        compression: int = 2,
        thread: int = 1,
        adapter_sequence: str = '',
        adapter_sequence_r2: str = '',
        poly_g_min_len: int = 10,
        poly_x_min_len: int = 10,
        overlap_len_require: int = 30,
        overlap_diff_limit: int = 5,
        overlap_diff_percent_limit: int = 20,
        correction: bool = False,
        cut_window_size: int = 4,
        cut_mean_quality: int = 20,
        cut_front: bool = False,
        cut_tail: bool = False,
        cut_right: bool = False
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    output_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    with tempfile.TemporaryDirectory as tmp:
        report_dir = os.path.join(tmp.name, 'reports')
        os.makedirs(report_dir, exist_ok=True)
        for sample_id, row in sequences.manifest.iterrows():
            input_fp = row['forward']
            output_fp = os.path.join(output_sequences.path, os.path.basename(row['forward']))
            report_fp = os.path.join(report_dir, f'{sample_id}.html')
            cmd = [
                'fastp',
                '--in1', input_fp,
                '--out1', output_fp,
                '--html', report_fp,
                '--trim_front1', str(trim_front1),
                '--trim_tail1', str(trim_tail1),
                '--max_len1', str(max_len1),
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
                '--overlap_diff_percent_limit', str(overlap_diff_percent_limit),
                '--disable_quality_filtering', str(disable_quality_filtering)
            ]
            if correction:
                cmd.append('--correction')
            if cut_front or cut_tail or cut_right:
                cmd.extend([
                    '--cut_window_size', str(cut_window_size),
                    '--cut_mean_quality', str(cut_mean_quality)
                ])
                if cut_front:
                    cmd.append('--cut_front')
                if cut_tail:
                    cmd.append('--cut_tail')
                if cut_right:
                    cmd.append('--cut_right')
            if 'reverse' in row and row['reverse'] is not None:
                input_fp2 = row['reverse']
                output_fp2 = os.path.join(output_sequences.path, os.path.basename(row['reverse']))
                cmd.extend([
                    '--in2', input_fp2,
                    '--out2', output_fp2,
                    '--trim_front2', str(trim_front2),
                    '--trim_tail2', str(trim_tail2),
                    '--max_len2', str(max_len2),
                    '--adapter_sequence_r2', adapter_sequence_r2
                ])
            run_command(cmd)
        # aggregate_fastp_reports(report_dir, output_sequences.path)
    return output_sequences
