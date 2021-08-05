from eva_vcf_merge.merge import horizontal_merge


def test_horizontal_merge(unique_samples_vcfs, merge_config):
    vcfs = {'name': unique_samples_vcfs}
    horizontal_merge(vcfs, merge_config)
    # TODO assert...


def test_horizontal_merge_multiple_groups(unique_samples_vcfs, merge_config):
    pass
