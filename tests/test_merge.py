import glob
import os

from ebi_eva_common_pyutils.command_utils import run_command_with_output


def assert_all_files_present(filenames):
    assert all(os.path.exists(f) for f in filenames)


def test_horizontal_merge(vcf_merger, unique_samples_vcfs):
    vcfs = {'alias': unique_samples_vcfs}
    filenames = vcf_merger.horizontal_merge(vcfs, resume=False)
    assert filenames == {'alias': os.path.join(vcf_merger.output_dir, 'alias_merged.vcf.gz')}
    assert_all_files_present(filenames.values())


def test_horizontal_merge_multiple_groups(vcf_merger, unique_samples_vcfs, unique_samples_vcfs_2):
    vcfs = {'1': unique_samples_vcfs, '2': unique_samples_vcfs_2}
    filenames = vcf_merger.horizontal_merge(vcfs, resume=False)
    assert filenames == {
        '1': os.path.join(vcf_merger.output_dir, '1_merged.vcf.gz'),
        '2': os.path.join(vcf_merger.output_dir, '2_merged.vcf.gz')
    }
    assert_all_files_present(filenames.values())


def test_vertical_merge(vcf_merger, same_samples_vcfs):
    vcfs = {'alias': same_samples_vcfs}
    filenames = vcf_merger.vertical_merge(vcfs, resume=False)
    assert filenames == {'alias': os.path.join(vcf_merger.output_dir, 'alias_merged.vcf.gz')}
    assert_all_files_present(filenames.values())


# TODO fails when run as suite but not alone... two of 0/0 concat + index but nothing else....
def test_concat_uninterrupted(vcf_merger, many_vcfs_to_concat):
    #   s0.vcf.gz   s1.vcf.gz   s2.vcf.gz   s3.vcf.gz   s4.vcf.gz
    #       \           /           \           /
    #        s01.vcf.gz               s23.vcf.gz        s4.vcf.gz       <-------- Stage 0
    #               \                      /
    #                \                   /
    #                 \                /
    #                   s0123.vcf.gz                    s4.vcf.gz       <-------- Stage 1
    #                           \                       /
    #                            \                    /
    #                             \                 /
    #                                 Final_merged                      <-------- Stage 2
    vcfs = {'many': many_vcfs_to_concat}
    vcf_merger.vertical_merge(vcfs, chunk_size=2, resume=False)
    stage_dirs = glob.glob(f'{vcf_merger.output_dir}/vertical_concat/stage*')
    assert len(stage_dirs) == 3
    output_vcf_from_multi_stage_concat = os.path.join(vcf_merger.output_dir, 'many_merged.vcf.gz')
    output_vcf_from_single_stage_concat = f'{vcf_merger.output_dir}/single_stage_concat_result.vcf.gz'
    run_command_with_output('Concatenate VCFs with single stage...', f"{vcf_merger.bcftools_binary} concat {' '.join(many_vcfs_to_concat)} "
                                                                     f"--allow-overlaps --remove-duplicates "
                                                                     f"-O z "
                                                                     f"-o {output_vcf_from_single_stage_concat}"
                            )
    diffs = run_command_with_output('Compare outputs from single and multi-stage concat processes...',
                                    f'bash -c "diff '
                                    f'<(zcat {output_vcf_from_single_stage_concat} | grep -v ^#) '
                                    f'<(zcat {output_vcf_from_multi_stage_concat} | grep -v ^#)"',
                                    return_process_output=True)
    assert diffs.strip() == ''
