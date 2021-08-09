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
import shutil

from ebi_eva_common_pyutils.nextflow import NextFlowPipeline, NextFlowProcess

from eva_vcf_merge.multistage import get_multistage_vertical_concat_pipeline
from eva_vcf_merge.utils import write_files_to_list


class VCFMerger:
    def __init__(self, bgzip_binary, bcftools_binary, nextflow_binary, nextflow_config, output_dir):
        self.bgzip_binary = bgzip_binary
        self.bcftools_binary = bcftools_binary
        self.nextflow_binary = nextflow_binary
        self.nextflow_config = nextflow_config
        self.output_dir = output_dir

    def horizontal_merge(self, vcf_groups, resume=True):
        """
        Merge groups of vcfs horizontally, i.e. by sample, using bcftools.

        :param vcf_groups: dict mapping a string (e.g. an analysis alias) to a group of vcf files to be merged
        :param resume: whether to resume pipeline (default true)
        :returns: dict of merged filenames
        """
        pipeline, merged_filenames = self.generate_horizontal_merge_pipeline(vcf_groups)
        workflow_file = os.path.join(self.output_dir, 'merge_workflow.nf')
        pipeline.run_pipeline(
            workflow_file_path=workflow_file,
            working_dir=self.output_dir,
            nextflow_binary_path=self.nextflow_binary,
            nextflow_config_path=self.nextflow_config,
            resume=resume
        )
        return merged_filenames

    def vertical_merge(self, vcf_groups, chunk_size=500, resume=True):
        """
        Merge groups of vcfs vertically, i.e. concatenation.

        :param vcf_groups: dict mapping a string (e.g. an analysis alias) to a group of vcf files to be merged
        :param chunk_size: number of vcfs to merge at once (default 500)
        :param resume: whether to resume pipeline (default true)
        :returns: dict of merged filenames
        """
        pipeline, merged_filenames = self.generate_vertical_merge_pipeline(vcf_groups, chunk_size)
        workflow_file = os.path.join(self.output_dir, "vertical_concat.nf")
        # TODO logging?
        pipeline.run_pipeline(
            workflow_file_path=workflow_file,
            working_dir=self.output_dir,
            nextflow_binary_path=self.nextflow_binary,
            nextflow_config_path=self.nextflow_config,
            resume=resume
        )
        # move merged files to output directory and rename
        for alias in merged_filenames:
            target_filename = os.path.join(self.output_dir, f'{alias}_merged.vcf.gz')
            shutil.move(merged_filenames[alias], target_filename)
            merged_filenames[alias] = target_filename
        return merged_filenames

    def generate_horizontal_merge_pipeline(self, vcf_groups):
        """
        Generate horizontal merge pipeline, including compressing and indexing VCFs.

        :param vcf_groups: dict mapping a string to a group of vcf files to be merged
        :return: complete NextflowPipeline and dict of merged filenames
        """
        dependencies = {}
        merged_filenames = {}
        for i, (alias, vcfs) in enumerate(vcf_groups.items()):
            deps, index_processes, compressed_vcfs = self.compress_and_index(i, vcfs)
            dependencies.update(deps)

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

    def generate_vertical_merge_pipeline(self, vcf_groups, chunk_size):
        """
        Generate vertical merge (concatenation) pipeline.

        :param vcf_groups: dict mapping a string to a group of vcf files to be merged
        :param chunk_size: number of vcfs to merge at once
        :return: complete NextFlowPipeline and dict of merged filenames
        """
        full_pipeline = NextFlowPipeline()
        merged_filenames = {}
        for i, (alias, vcfs) in enumerate(vcf_groups.items()):
            deps, index_processes, compressed_vcfs = self.compress_and_index(i, vcfs)
            compress_pipeline = NextFlowPipeline(deps)
            # TODO use alias in the merged filename somehow?
            concat_pipeline, merged_filename = get_multistage_vertical_concat_pipeline(
                vcf_files=compressed_vcfs,
                concat_chunk_size=chunk_size,
                concat_processing_dir=self.output_dir,
                bcftools_binary=self.bcftools_binary
            )
            pipeline = NextFlowPipeline.join_pipelines(compress_pipeline, concat_pipeline)
            full_pipeline = NextFlowPipeline.join_pipelines(full_pipeline, pipeline)
            merged_filenames[alias] = merged_filename
        return full_pipeline, merged_filenames

    def compress_and_index(self, alias_index, vcfs):
        """
        Bgzip-compress and CSI-index VCFs.

        :param vcfs:
        :return:
        """
        dependencies = {}
        index_processes = []
        compressed_vcfs = []
        for i, vcf in enumerate(vcfs):
            compress_process = None
            if not vcf.endswith('gz'):
                compress_process = NextFlowProcess(
                    process_name=f'compress_{alias_index}_{i}',
                    command_to_run=f'{self.bgzip_binary} -c {vcf} > {vcf}.gz'
                )
                vcf = f'{vcf}.gz'
            compressed_vcfs.append(vcf)
            index_process = NextFlowProcess(
                process_name=f'index_{alias_index}_{i}',
                command_to_run=f'{self.bcftools_binary} index -f -c {vcf}'
            )
            index_processes.append(index_process)
            # each file's index depends only on compress (if present)
            dependencies[index_process] = [compress_process] if compress_process else []
        # TODO preferably return a NextFlowPipeline rather than dependencies & final processes
        return dependencies, index_processes, compressed_vcfs
