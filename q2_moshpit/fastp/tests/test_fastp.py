# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
import unittest
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt
from q2_moshpit.fastp import run_fastp
from q2_moshpit.fastp.aggregate import aggregate_fastp_reports

class TestFastp(unittest.TestCase):

    def setUp(self):
        self.input_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
        self.output_sequences = CasavaOneEightSingleLanePerSampleDirFmt()

        # Create dummy input files
        for i in range(1, 4):
            input_fp = os.path.join(self.input_sequences.path, f'sample{i}.fastq.gz')
            with open(input_fp, 'w') as f:
                f.write('@SEQ_ID\nGATTTGGGGTTTCCCAGTTG\n+\nIIIIIIIIIIIIIIIIIIII\n')

    def test_run_fastp(self):
        output_sequences = run_fastp(self.input_sequences, trim_front1=5, trim_tail1=5, cut_window_size=4, cut_mean_quality=20, n_base_limit=5, length_required=15, qualified_quality_phred=15, unqualified_percent_limit=40.0, compression=2, thread=3)

        # Check if output files are created
        for i in range(1, 4):
            output_fp = os.path.join(output_sequences.path, f'sample{i}.fastq.gz')
            self.assertTrue(os.path.exists(output_fp))

        # Check if the output files are not empty
        for i in range(1, 4):
            output_fp = os.path.join(output_sequences.path, f'sample{i}.fastq.gz')
            self.assertGreater(os.path.getsize(output_fp), 0)

    def test_aggregate_fastp_reports(self):
        report_dir = 'test_reports'
        output_dir = 'test_output'
        os.makedirs(report_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)

        # Create dummy HTML reports
        for i in range(1, 4):
            report_fp = os.path.join(report_dir, f'sample{i}.html')
            with open(report_fp, 'w') as f:
                f.write(f'<html><body>Report for sample{i}</body></html>')

        aggregate_fastp_reports(report_dir, output_dir)

        # Check if the aggregated report is created
        aggregated_report_fp = os.path.join(output_dir, 'index.html')
        self.assertTrue(os.path.exists(aggregated_report_fp))

        # Check if the aggregated report contains the individual reports
        with open(aggregated_report_fp, 'r') as f:
            aggregated_content = f.read()
            for i in range(1, 4):
                self.assertIn(f'Report for sample{i}', aggregated_content)

if __name__ == '__main__':
    unittest.main()
