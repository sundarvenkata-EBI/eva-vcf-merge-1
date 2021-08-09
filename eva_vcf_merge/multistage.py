import math
import os

from ebi_eva_common_pyutils.nextflow import NextFlowPipeline, NextFlowProcess

from eva_vcf_merge.utils import write_files_to_list


def get_multistage_vertical_concat_pipeline(
        vcf_files,
        concat_processing_dir,
        concat_chunk_size,
        bcftools_binary,
        stage=0,
        prev_stage_processes=[],
        pipeline=NextFlowPipeline()
):
    """
    # Generate Nextflow pipeline for multi-stage VCF concatenation of 5 VCF files with 2-VCFs concatenated at a time (CONCAT_CHUNK_SIZE=2)
    # For illustration purposes only. Usually the CONCAT_CHUNK_SIZE is much higher (ex: 500).
    #
    #		    vcf1		    vcf2		vcf3		    vcf4		vcf5
    #               \		     /		       \		      /
    # Stage0:	     \		   /		        \		    /
    # -------	      \	     /		             \	      /
    #	    	vcf1_2=concat(vcf1,vcf2)	vcf3_4=concat(vcf3,vcf4)	vcf5    <---- 3 batches of concat in stage 0
    #		 		    \ 		                /
    # Stage1:	  		 \		              /
    # -------	   		  \	                /
    #				vcf1_2_3_4=concat(vcf1_2,vcf3_4)		            vcf5    <---- 2 batches of concat in stage 1
    #	 					      \ 		                            /
    # Stage2:	  		 		   \		 	                      /
    # -------	   		  			\	                            /
    #						      vcf1_2_3_4_5=concat(vcf1_2_3_4,vcf5)          <----- Final result
    """
    # If we are left with only one file, this means we have reached the last concat stage
    if len(vcf_files) == 1:
        return pipeline, vcf_files[0]

    num_batches_in_stage = math.ceil(len(vcf_files) / concat_chunk_size)
    curr_stage_processes = []
    output_vcf_files_from_stage = []
    for batch in range(0, num_batches_in_stage):
        # split files in the current stage into chunks based on concat_chunk_size
        files_in_batch = vcf_files[(concat_chunk_size * batch):(concat_chunk_size * (batch + 1))]
        files_to_concat_list = write_files_to_concat_list(files_in_batch, stage, batch, concat_processing_dir)
        output_vcf_file = get_output_vcf_file_name(stage, batch, concat_processing_dir)

        # TODO log file?  log_file_name = os.path.join(concat_processing_dir, f"{concat_stage_batch_name}.log")
        concat_process = NextFlowProcess(
            process_name=f"concat_stage{stage}_batch{batch}",
            command_to_run=f"{bcftools_binary} concat --allow-overlaps --remove-duplicates "
                           f"--file-list {files_to_concat_list} -o {output_vcf_file} -O z"
        )
        index_process = NextFlowProcess(
            process_name=f"index_stage{stage}_batch{batch}",
            command_to_run=f"{bcftools_binary} index --csi {output_vcf_file}"
        )
        # index depends only on this batch's concat
        pipeline.add_dependencies({index_process: [concat_process]})
        # next stage requires indexing to be complete from this stage
        curr_stage_processes.append(index_process)
        output_vcf_files_from_stage.append(output_vcf_file)

        # Concatenation batch in a given stage will have to wait until the completion of
        # n batches in the previous stage where n = concat_chunk_size
        # Ex: In the illustration above stage 1/batch 0 depends on completion of stage 0/batch 0 and stage 0/batch 1
        # While output of any n batches from the previous stage can be worked on as they become available,
        # having a predictable formula simplifies pipeline generation and troubleshooting
        prev_stage_dependencies = prev_stage_processes[(concat_chunk_size * batch):(concat_chunk_size * (batch + 1))]
        pipeline.add_dependencies({concat_process: prev_stage_dependencies})
    prev_stage_processes = curr_stage_processes

    return get_multistage_vertical_concat_pipeline(
        output_vcf_files_from_stage,
        concat_processing_dir, concat_chunk_size,
        bcftools_binary,
        stage=stage+1,
        prev_stage_processes=prev_stage_processes,
        pipeline=pipeline
    )


def write_files_to_concat_list(files_to_concat, concat_stage, concat_batch, concat_processing_dir):
    """
    Write the list of files to be concatenated for a given stage and batch
    """
    output_dir = get_concat_output_dir(concat_stage, concat_processing_dir)
    alias = f'batch{concat_batch}'
    return write_files_to_list(files_to_concat, alias, output_dir)


def get_concat_output_dir(concat_stage_index: int, concat_processing_dir: str):
    """
    Get the file name with the list of files to be concatenated for a given stage and batch in the concatenation process
    """
    return os.path.join(concat_processing_dir, "vertical_concat", f"stage_{concat_stage_index}")


def get_output_vcf_file_name(concat_stage_index: int, concat_batch_index: int, concat_processing_dir: str):
    return os.path.join(get_concat_output_dir(concat_stage_index, concat_processing_dir),
                        f"concat_output_stage{concat_stage_index}_batch{concat_batch_index}.vcf.gz")
