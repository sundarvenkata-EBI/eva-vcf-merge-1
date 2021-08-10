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

import itertools
from enum import Enum

import pysam

from eva_vcf_merge.utils import are_all_elements_unique


class SampleSetType(Enum):
    # Same sample names in the same order, e.g. [(A, B), (A, B)]
    SINGLE_SET = 'single set'
    # Same sample names but in different orders, e.g. [(A, B), (B, A)]
    UNSORTED_SINGLE_SET = 'unsorted single set'
    # Distinct non-overlapping sample names, e.g. [(A, B), (C, D, E)]
    UNIQUE_SETS = 'unique sets'
    # Overlapping sample names, e.g. [(A, B), (C, D, A)]
    OVERLAPPING_SETS = 'overlapping sets'


class MergeType(Enum):
    # Horizontally mergeable VCFs have non-overlapping sample sets and can be combined to create
    # a single multi-sample file.
    HORIZONTAL = 'horizontal'
    # Vertically mergeable VCFs have identical sample sets and can be combined to create a single
    # multi-chromosome or multi-variant file.
    VERTICAL = 'vertical'


def detect_merge_type(vcf_files):
    """
    Detect what type of merge should be used for a list of VCFs.

    :param vcf_files: list of vcf filepaths
    :return: MergeType or None if it could not be detected
    """
    file_to_sample_names = {}
    # retrieve all the sample_names from the VCF files
    for file_path in vcf_files:
        file_to_sample_names[file_path] = get_samples_from_vcf(file_path)
    # Check that all the samples are the same and in the same order to enable horizontal merging
    sample_info = compare_sample_sets(file_to_sample_names.values())
    if sample_info == SampleSetType.UNIQUE_SETS:
        return MergeType.HORIZONTAL
    elif sample_info == SampleSetType.SINGLE_SET:
        return MergeType.VERTICAL
    return None


def get_samples_from_vcf(vcf_file):
    """Get the list of samples present in a single VCF file."""
    with pysam.VariantFile(vcf_file, 'r') as vcf_in:
        samples = list(vcf_in.header.samples)
    return samples


def compare_sample_sets(list_of_sample_names):
    """
    Compare sample names from different VCFs to determine whether they are identical, distinct, or overlapping.

    :param list_of_sample_names: list of list of sample names in each VCF
    :return: SampleSetType or None if it could not be detected
    """
    set_of_sample_names = set(tuple(s) for s in list_of_sample_names)
    set_of_sample_names_sorted = set([tuple(sorted(list(s))) for s in set_of_sample_names])
    if len(set_of_sample_names) == 1:
        return SampleSetType.SINGLE_SET
    elif len(set_of_sample_names_sorted) == 1:
        return SampleSetType.UNSORTED_SINGLE_SET
    elif are_all_elements_unique(itertools.chain(*list_of_sample_names)):
        return SampleSetType.UNIQUE_SETS
    elif len(set_of_sample_names_sorted) == len(list_of_sample_names):
        return SampleSetType.OVERLAPPING_SETS
    return None
