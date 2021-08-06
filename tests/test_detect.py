from eva_vcf_merge.detect import compare_sample_sets, SampleSetType, MergeType, detect_merge_type


def test_detect_merge_type(unique_samples_vcfs, same_samples_vcfs, overlapping_samples_vcfs):
    assert detect_merge_type(unique_samples_vcfs) == MergeType.HORIZONTAL
    assert detect_merge_type(same_samples_vcfs) == MergeType.VERTICAL
    assert detect_merge_type(overlapping_samples_vcfs) == None


def test_compare_sample_sets():
    assert compare_sample_sets([['S1', 'S2', 'S3'], ['S1', 'S2', 'S3']]) == SampleSetType.SINGLE_SET
    assert compare_sample_sets([['S1', 'S2', 'S3'], ['S1', 'S3', 'S2']]) == SampleSetType.UNSORTED_SINGLE_SET
    assert compare_sample_sets([['S1'], ['S2'], ['S3']]) == SampleSetType.UNIQUE_SETS
    assert compare_sample_sets([['S1'], ['S2'], ['S3', 'S4']]) == SampleSetType.UNIQUE_SETS
    assert compare_sample_sets([['S1'], ['S2'], ['S3', 'S1']]) == SampleSetType.OVERLAPPING_SETS
    assert compare_sample_sets([['S1'], ['S2'], ['S3', 'S3']]) == SampleSetType.OVERLAPPING_SETS
