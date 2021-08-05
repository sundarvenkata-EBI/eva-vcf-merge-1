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

from eva_vcf_merge.config import MergeConfig
from eva_vcf_merge.utils import write_files_to_list


# TODO deal with invalid process names
def generate_pipeline(vcf_groups, bgzip_binary, bcftools_binary, output_dir):
    dependencies = {}
    for alias, vcfs in vcf_groups.items():
        index_processes = []
        compressed_vcfs = []
        for i, vcf in enumerate(vcfs):
            index_process = NextFlowProcess(
                process_name=f'index_{alias}_{i}',
                command_to_run=f'{bcftools_binary} index -c {vcf}.gz'
            )
            index_processes.append(index_process)
            if vcf.endswith('gz'):
                compressed_vcfs.append(vcf)
            else:
                compress_process = NextFlowProcess(
                    process_name=f'compress_{alias}_{i}',
                    command_to_run=f'{bgzip_binary} -c {vcf} > {vcf}.gz'
                )
                # each file's index depends only on compress (if present)
                dependencies[index_process] = [compress_process]
                compressed_vcfs.append(f'{vcf}.gz')

        list_filename = write_files_to_list(compressed_vcfs, alias, output_dir)
        merged_filename = os.path.join(output_dir, f'{alias}_merged.vcf.gz')
        merge_process = NextFlowProcess(
            process_name=f'merge_{alias}',
            command_to_run=f'{bcftools_binary} merge --merge all --file-list {list_filename} '
                           f'--threads 3 -O z -o {merged_filename}'
        )
        # each alias's merge process depends on all index processes
        dependencies[merge_process] = index_processes
    return NextFlowPipeline(dependencies)


def horizontal_merge(vcf_groups, cfg=MergeConfig()):
    """
    Merge groups of vcfs horizontally, i.e. by sample, using bcftools.

    :param vcf_groups: dict mapping a string (e.g. an analysis alias) to a group of vcf files to be merged
    :param cfg: MergeConfig with necessary configuration
    """
    pipeline = generate_pipeline(vcf_groups, cfg.bgzip_binary, cfg.bcftools_binary, cfg.output_dir)
    workflow_file = os.path.join(cfg.output_dir, 'merge_workflow.nf')
    pipeline.run_pipeline(
        workflow_file_path=workflow_file,
        working_dir=cfg.output_dir,
        nextflow_binary_path=cfg.nextflow_binary,
        nextflow_config_path=cfg.nextflow_config
    )


def vertical_concat(vcf_groups, cfg=MergeConfig()):
    """
    Merge groups of vcfs vertically, i.e. concatenation.

    :param vcf_groups: dict mapping a string (e.g. an analysis alias) to a group of vcf files to be merged
    :param cfg: MergeConfig with necessary configuration
    """
    raise NotImplementedError('Vertical concatenation not yet implemented.')
