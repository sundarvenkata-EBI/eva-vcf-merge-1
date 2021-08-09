import glob
import os
import shutil

import pytest

from eva_vcf_merge.merge import VCFMerger

tests_dir = os.path.dirname(__file__)
resources_dir = os.path.join(tests_dir, 'resources')


@pytest.fixture
def vcf_merger():
    # create output directory
    output_dir = os.path.join(tests_dir, 'output')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # use standard executables for testing
    yield VCFMerger(
        bgzip_binary='bgzip',
        bcftools_binary='bcftools',
        nextflow_binary='nextflow',
        nextflow_config=None,
        output_dir=output_dir
    )
    # clean up after merge tests
    shutil.rmtree(output_dir)
    for ext in ('*.gz', '*.tbi', '*.csi'):
        for fn in glob.glob(os.path.join(resources_dir, ext)):
            os.remove(fn)


@pytest.fixture
def unique_samples_vcfs():
    return [
        os.path.join(resources_dir, 'chr1_samplesAB.vcf'),
        os.path.join(resources_dir, 'chr1_samplesCD.vcf')
    ]


@pytest.fixture
def unique_samples_vcfs_2():
    return [
        os.path.join(resources_dir, 'chr2_samplesAB.vcf'),
        os.path.join(resources_dir, 'chr2_samplesCDE.vcf')
    ]


@pytest.fixture
def same_samples_vcfs():
    return [
        os.path.join(resources_dir, 'chr1_samplesAB.vcf'),
        os.path.join(resources_dir, 'chr2_samplesAB.vcf')
    ]


@pytest.fixture
def overlapping_samples_vcfs():
    return [
        os.path.join(resources_dir, 'chr1_samplesCD.vcf'),
        os.path.join(resources_dir, 'chr2_samplesCDE.vcf')
    ]


@pytest.fixture
def many_vcfs_to_concat():
    return [
        os.path.join(resources_dir, 'concat', 's0.vcf.gz'),
        os.path.join(resources_dir, 'concat', 's1.vcf.gz'),
        os.path.join(resources_dir, 'concat', 's2.vcf.gz'),
        os.path.join(resources_dir, 'concat', 's3.vcf.gz'),
        os.path.join(resources_dir, 'concat', 's4.vcf.gz')
    ]
