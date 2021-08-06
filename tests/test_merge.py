import os

from eva_vcf_merge.merge import horizontal_merge


def test_horizontal_merge(unique_samples_vcfs, merge_config):
    vcfs = {'alias': unique_samples_vcfs}
    horizontal_merge(vcfs, merge_config)
    assert os.path.exists(os.path.join(merge_config.output_dir, 'alias_merged.vcf.gz'))


def test_horizontal_merge_multiple_groups(unique_samples_vcfs, unique_samples_vcfs_2, merge_config):
    vcfs = {'1': unique_samples_vcfs, '2': unique_samples_vcfs_2}
    horizontal_merge(vcfs, merge_config)
    assert os.path.exists(os.path.join(merge_config.output_dir, '1_merged.vcf.gz'))
    assert os.path.exists(os.path.join(merge_config.output_dir, '2_merged.vcf.gz'))
