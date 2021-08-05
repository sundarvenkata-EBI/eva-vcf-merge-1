import glob
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


def clean_up_resources():
    # cleans up any new files created by tests in resources directory
    for ext in ('*.gz', '*.tbi', '*.csi'):
        for fn in glob.glob(os.path.join(resources_dir, ext)):
            os.remove(fn)


@pytest.fixture
def unique_samples_vcfs():
    yield [
        os.path.join(resources_dir, 'chr1_samplesAB.vcf'),
        os.path.join(resources_dir, 'chr1_samplesCD.vcf')
    ]
    clean_up_resources()


@pytest.fixture
def unique_samples_vcfs_2():
    yield [
        os.path.join(resources_dir, 'chr2_samplesAB.vcf'),
        os.path.join(resources_dir, 'chr2_samplesCDE.vcf')
    ]
    clean_up_resources()


@pytest.fixture
def same_samples_vcfs():
    yield [
        os.path.join(resources_dir, 'chr1_samplesAB.vcf'),
        os.path.join(resources_dir, 'chr2_samplesAB.vcf')
    ]
    clean_up_resources()


@pytest.fixture
def overlapping_samples_vcfs():
    yield [
        os.path.join(resources_dir, 'chr1_samplesCD.vcf'),
        os.path.join(resources_dir, 'chr2_samplesCDE.vcf')
    ]
    clean_up_resources()
