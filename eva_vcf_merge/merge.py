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


def generate_pipeline(vcf_groups, bcftools_binary, output_dir):
    dependencies = {}
    for label, vcfs in vcf_groups.items():
        merge_process = NextFlowProcess(
            process_name=f'merge_{label}',
            command_to_run=f'{bcftools_binary} merge --merge all --file-list all_files_{label}.list '
                           f'--threads 3 -O z -o {os.path.join(output_dir, f"{label}_merged.vcf.gz")}'
        )
        qc_process = NextFlowProcess(
            process_name=f'qc_{label}',
            command_to_run=''  # TODO
        )
        # QC depends on merge, but otherwise the groups are independent.
        dependencies[qc_process] = [merge_process]
    return NextFlowPipeline(dependencies)


def horizontal_merge(vcf_groups, cfg):
    """
    Merge groups of vcfs horizontally, i.e. by sample, using bcftools.

    :param vcf_groups: dict mapping a string (e.g. an analysis alias) to a group of vcf files to be merged
    :param cfg: MergeConfig with necessary configuration
    """
    pipeline = generate_pipeline(vcf_groups, cfg.bcftools_binary, cfg.output_dir)
    workflow_file = os.path.join(cfg.output_dir, 'merge_workflow.nf')
    work_dir = os.path.join(cfg.output_dir, 'work')
    pipeline.run_pipeline(
        workflow_file_path=workflow_file,
        working_dir=work_dir,
        nextflow_binary_path=cfg.nextflow_binary,
        nextflow_config_path=cfg.nextflow_config
    )
    # clean up work dir
    shutil.rmtree(work_dir)
