import os
import shutil

import pytest

from eva_vcf_merge.config import MergeConfig

tests_dir = os.path.dirname(__file__)
resources_dir = os.path.join(tests_dir, 'resources')


@pytest.fixture
def merge_config():
    # creates and cleans up output directory
    output_dir = os.path.join(tests_dir, 'output')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    yield MergeConfig(output_dir=output_dir)
    shutil.rmtree(output_dir)


@pytest.fixture
def horizontal_merge_vcfs():
    return [
        os.path.join(resources_dir, 'chr1_samplesAB.vcf'),
        os.path.join(resources_dir, 'chr1_samplesCD.vcf')
    ]


@pytest.fixture
def vertical_merge_vcfs():
    return [
        os.path.join(resources_dir, 'chr1_samplesAB.vcf'),
        os.path.join(resources_dir, 'chr2_samplesAB.vcf')
    ]
