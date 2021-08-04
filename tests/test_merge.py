from eva_vcf_merge.merge import horizontal_merge


def test_horizontal_merge(horizontal_merge_vcfs, merge_config):
    vcfs = {'name': horizontal_merge_vcfs}
    horizontal_merge(vcfs, merge_config)
