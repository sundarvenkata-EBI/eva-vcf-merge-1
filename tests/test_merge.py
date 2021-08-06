import os


def test_horizontal_merge(vcf_merger, unique_samples_vcfs):
    vcfs = {'alias': unique_samples_vcfs}
    vcf_merger.horizontal_merge(vcfs)
    assert os.path.exists(os.path.join(vcf_merger.output_dir, 'alias_merged.vcf.gz'))


def test_horizontal_merge_multiple_groups(vcf_merger, unique_samples_vcfs, unique_samples_vcfs_2):
    vcfs = {'1': unique_samples_vcfs, '2': unique_samples_vcfs_2}
    vcf_merger.horizontal_merge(vcfs)
    assert os.path.exists(os.path.join(vcf_merger.output_dir, '1_merged.vcf.gz'))
    assert os.path.exists(os.path.join(vcf_merger.output_dir, '2_merged.vcf.gz'))
