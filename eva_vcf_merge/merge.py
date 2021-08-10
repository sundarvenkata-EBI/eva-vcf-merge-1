# Copyright 2021 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os

from ebi_eva_common_pyutils.nextflow import NextFlowPipeline, NextFlowProcess

from eva_vcf_merge.utils import write_files_to_list


class VCFMerger:
    def __init__(self, bgzip_binary, bcftools_binary, nextflow_binary, nextflow_config, output_dir):
        self.bgzip_binary = bgzip_binary
        self.bcftools_binary = bcftools_binary
        self.nextflow_binary = nextflow_binary
        self.nextflow_config = nextflow_config
        self.output_dir = output_dir

    def horizontal_merge(self, vcf_groups):
        """
        Merge groups of vcfs horizontally, i.e. by sample, using bcftools.

        :param vcf_groups: dict mapping a string (e.g. an analysis alias) to a group of vcf files to be merged
        :returns: dict of merged filenames
        """
        pipeline, merged_filenames = self.generate_horizontal_merge_pipeline(vcf_groups)
        workflow_file = os.path.join(self.output_dir, 'merge_workflow.nf')
        pipeline.run_pipeline(
            workflow_file_path=workflow_file,
            working_dir=self.output_dir,
            nextflow_binary_path=self.nextflow_binary,
            nextflow_config_path=self.nextflow_config
        )
        return merged_filenames

    def vertical_concat(self, vcf_groups):
        """
        Merge groups of vcfs vertically, i.e. concatenation.

        :param vcf_groups: dict mapping a string (e.g. an analysis alias) to a group of vcf files to be merged
        """
        raise NotImplementedError('Vertical concatenation not yet implemented.')

    def generate_horizontal_merge_pipeline(self, vcf_groups):
        """
        Generate horizontal merge pipeline, including compressing and indexing VCFs.

        :param vcf_groups: dict mapping a string to a group of vcf files to be merged
        :return: complete NextflowPipeline and dict of merged filenames
        """
        dependencies = {}
        merged_filenames = {}
        for i, (alias, vcfs) in enumerate(vcf_groups.items()):
            index_processes = []
            compressed_vcfs = []
            for j, vcf in enumerate(vcfs):
                # compress vcf if needed
                compress_process = None
                if not vcf.endswith('gz'):
                    compress_process = NextFlowProcess(
                        process_name=f'compress_{i}_{j}',
                        command_to_run=f'{self.bgzip_binary} -c {vcf} > {vcf}.gz'
                    )
                    vcf = f'{vcf}.gz'

                compressed_vcfs.append(vcf)
                index_process = NextFlowProcess(
                    process_name=f'index_{i}_{j}',
                    command_to_run=f'{self.bcftools_binary} index -f -c {vcf}'
                )
                index_processes.append(index_process)
                # each file's index depends only on compress (if present)
                dependencies[index_process] = [compress_process] if compress_process else []

            list_filename = write_files_to_list(compressed_vcfs, alias, self.output_dir)
            merged_filename = os.path.join(self.output_dir, f'{alias}_merged.vcf.gz')
            merge_process = NextFlowProcess(
                process_name=f'merge_{i}',
                command_to_run=f'{self.bcftools_binary} merge --merge all --file-list {list_filename} '
                               f'--threads 3 -O z -o {merged_filename}'
            )
            # each alias's merge process depends on all index processes
            dependencies[merge_process] = index_processes
            merged_filenames[alias] = merged_filename

        return NextFlowPipeline(dependencies), merged_filenames
