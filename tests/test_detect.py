from eva_vcf_merge.detect import compare_sample_sets, SampleSetType, MergeType, detect_merge_type


def test_detect_merge_type(horizontal_merge_vcfs, vertical_merge_vcfs):
    assert detect_merge_type(horizontal_merge_vcfs) == MergeType.HORIZONTAL
    assert detect_merge_type(vertical_merge_vcfs) == MergeType.VERTICAL


def test_compare_sample_sets():
    assert compare_sample_sets([['S1', 'S2', 'S3'], ['S1', 'S2', 'S3']]) == SampleSetType.SINGLE_SET
    assert compare_sample_sets([['S1', 'S2', 'S3'], ['S1', 'S3', 'S2']]) == SampleSetType.UNSORTED_SINGLE_SET
    assert compare_sample_sets([['S1'], ['S2'], ['S3']]) == SampleSetType.UNIQUE_SETS
    assert compare_sample_sets([['S1'], ['S2'], ['S3', 'S4']]) == SampleSetType.UNIQUE_SETS
    assert compare_sample_sets([['S1'], ['S2'], ['S3', 'S1']]) == SampleSetType.OVERLAPPING_SETS
    assert compare_sample_sets([['S1'], ['S2'], ['S3', 'S3']]) == SampleSetType.OVERLAPPING_SETS
