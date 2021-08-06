import os


def assert_all_files_present(filenames):
    assert all(os.path.exists(f) for f in filenames)


def test_horizontal_merge(vcf_merger, unique_samples_vcfs):
    vcfs = {'alias': unique_samples_vcfs}
    filenames = vcf_merger.horizontal_merge(vcfs)
    assert filenames == {'alias': os.path.join(vcf_merger.output_dir, 'alias_merged.vcf.gz')}
    assert_all_files_present(filenames.values())


def test_horizontal_merge_multiple_groups(vcf_merger, unique_samples_vcfs, unique_samples_vcfs_2):
    vcfs = {'1': unique_samples_vcfs, '2': unique_samples_vcfs_2}
    filenames = vcf_merger.horizontal_merge(vcfs)
    assert filenames == {
        '1': os.path.join(vcf_merger.output_dir, '1_merged.vcf.gz'),
        '2': os.path.join(vcf_merger.output_dir, '2_merged.vcf.gz')
    }
    assert_all_files_present(filenames.values())
