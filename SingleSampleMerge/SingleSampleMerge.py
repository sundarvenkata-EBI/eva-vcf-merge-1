# Copyright 2017 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import SSMergeCommonUtils
import datetime
import gzip
import os
import sys
import traceback

from cassandra.cluster import Cluster
from cassandra.query import BatchStatement, BatchType

from pyspark import SparkConf, SparkContext, SQLContext


class SingleSampleMerge:
    def __init__(self, *args):
        self.bcf_tools_dir, self.study_name, self.study_vcf_file_name, self.cassandraNodeIPs, self.keyspaceName, \
            self.headerTableName, self.variantTableName, self.sampleInsertLogTableName, self.uniquePosTableName = args

        self.cassandraCluster = Cluster(self.cassandraNodeIPs)
        self.cassandraSession = self.cassandraCluster.connect(self.keyspaceName)
        self.sampleName = SSMergeCommonUtils.getSampleName(self.bcf_tools_dir, self.study_vcf_file_name)
        self.filterVCFFileName = os.path.dirname(self.study_vcf_file_name) + os.path.sep + "{0}_filtervcf.gz"\
            .format(self.sampleName)
        self.sample_proc_status = self.get_sample_processed_status()

        self.variantInsertStmt = self.cassandraSession.prepare(
            "insert into {0}.{1} "
            "(chrom,chunk,start_pos,ref,alt,qual,filter,info,"
            "sampleinfoformat,sampleinfo,var_id,var_uniq_id,sampleName) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)".format(self.keyspaceName, self.variantTableName))
        self.headerInsertStmt = self.cassandraSession.prepare("insert into {0}.{1} (samplename, header) VALUES (?,?)"
                                                              .format(self.keyspaceName, self.headerTableName))
        self.uniquePosInsertStmt = self.cassandraSession.prepare(
            "insert into {0}.{1} (chrom, chunk, in_src, start_pos) VALUES (?,?,?,?)"
            .format(self.keyspaceName, self.uniquePosTableName))

    def process_study_file(self):
        """
        Process the study file that this specific SingleSampleMerge object has been given
        
        :return: Tuple - Total variants in sample, Error messages if any  
        """
        err_file_contents, return_err_msg, cassandra_cluster, cassandra_session = None, None, None, None
        total_variants_in_sample = 0

        if self.sample_proc_status == "variants_inserted":
            # return 0,0,None
            return (-1, "Already processed sample: {0} for study: {1}"
                    .format(self.sampleName, self.study_name))

        try:
            # Get the sample prefix and the base directory of the file that we are operating on
            base_dir = os.path.dirname(self.study_vcf_file_name)

            filter_err_file_name = base_dir + os.path.sep + "{0}_filtering_err.txt".format(self.sampleName)

            # region Run bcftools filter command to filter out monomorphic references
            if self.sample_proc_status != "variants_filtered":

                if os.path.isfile(self.filterVCFFileName):
                    os.remove(self.filterVCFFileName)
                os.system("""({0}/bin/bcftools filter -e ALT=\\'.\\' {1} | gzip) 1>{2} 2>{3}""".
                          format(self.bcf_tools_dir, self.study_vcf_file_name, self.filterVCFFileName,
                                 filter_err_file_name))
                err_file_contents = SSMergeCommonUtils.getErrFileContents(filter_err_file_name)
                if err_file_contents:
                    filter_command_result = -1
                else:
                    filter_command_result = 0
                    self.update_processing_status("variants_filtered", 0)
            else:
                filter_command_result = 0
            # endregion

            # region Insert the variants into Cassandra from the filtered VCF files obtained from above
            if filter_command_result != 0:
                return_err_msg = "Failed to process {0} due to error: {1}".\
                    format(self.study_vcf_file_name, err_file_contents)
            else:
                total_variants_in_sample = self.cassandra_insert()
                os.system("echo {0} > {1}_filtered_variant_count.txt".format(str(total_variants_in_sample),
                                                                             self.sampleName))
            # endregion

            # region Update status to indicate that all variants have been inserted into Cassandra
            if not return_err_msg and total_variants_in_sample != 0:
                self.update_processing_status("variants_inserted", total_variants_in_sample)
                # endregion

        except Exception, ex:
            return_err_msg = "Error in processing file:{0}".format(self.study_vcf_file_name) + os.linesep + ex.message \
                           + os.linesep + traceback.format_exc()
        finally:
            if cassandra_cluster is not None and cassandra_session is not None:
                try:
                    cassandra_session.shutdown()
                    cassandra_cluster.shutdown()
                except Exception, ex:
                    print(ex.message)
            if return_err_msg:
                return -1, return_err_msg
            return total_variants_in_sample, None

    def write_variants_to_cassandra(self, lines_to_write):
        """
        Write a set of variant lines to Cassandra
        
        :param lines_to_write: Set of variant lines to write from a given sample
        :type lines_to_write: list[str]        
        :return: None
        """
        batch = BatchStatement(BatchType.UNLOGGED)

        for line in lines_to_write:
            line_comps = [x.strip() for x in line.split("\t")]
            chromosome, position, rs_id, ref, alt, qual, qual_filter, info = line_comps[:8]
            sample_info_format, sample_info = None, None
            if len(line_comps) > 8:
                sample_info_format, sample_info = line_comps[8:]

            # Chromosome + Chunk combination determines which node a variant is written to
            position = long(float(position))
            chunk = int(position / SSMergeCommonUtils.CHR_POS_CHUNKSIZE)

            # Insert the chromosome position that was scanned
            batch.add(self.uniquePosInsertStmt.bind([chromosome, chunk, 1, position]))
            # Also, insert the normalized chromosome positions for each alternate allele
            for altAllele in alt.split(","):
                revised_start_pos = SSMergeCommonUtils.getNormalizedStartPos(chromosome, position, ref, altAllele)
                if revised_start_pos != position:
                    revised_chunk = int(revised_start_pos / SSMergeCommonUtils.CHR_POS_CHUNKSIZE)
                    batch.add(self.uniquePosInsertStmt.bind([chromosome, revised_chunk, 0, revised_start_pos]))

            # Unique Variant ID is chromosome_12-digit zero padded start_MD5(REF_ALT)_20-character samplename
            variant_id = SSMergeCommonUtils.getUniqueVariantID(self.sampleName, chromosome, position, ref, alt)
            # Insert the variant information
            bound_stmt = self.variantInsertStmt.bind([chromosome, chunk, position, ref, alt,
                                                      qual, qual_filter, info, sample_info_format, sample_info,
                                                      rs_id, variant_id, self.sampleName])
            batch.add(bound_stmt)

        # Batch execute the insert statements added in the loop above
        self.cassandraSession.execute(batch, timeout=SSMergeCommonUtils.BATCH_WRITE_TIMEOUT_IN_SECS)

    def write_header_to_cassandra(self, header_lines):
        """
        Write a set of header lines to Cassandra
    
        :param header_lines: Set of header lines to write from a given sample
        :type header_lines: str        
        :return: None
        """
        self.cassandraSession.execute(self.headerInsertStmt.bind([self.sampleName, header_lines]),
                                      timeout=SSMergeCommonUtils.BATCH_WRITE_TIMEOUT_IN_SECS)

    def cassandra_insert(self):
        """
        Insert the contents of a VCF file: variants and headers into relevant Cassandra tables        
        
        :return: Total number of variants that were scanned from the VCF file
        :rtype: long
        """        
        tot_vars_in_sample = 0
        header_lines, sample_name = "", ""
        lines_to_write = []

        with gzip.open(self.filterVCFFileName) as vcfFileHandle:
            for line in vcfFileHandle:
                line = line.strip()
                if line.startswith("#"):
                    header_lines += (line + os.linesep)
                    if line.startswith("#CHROM"):
                        self.write_header_to_cassandra(header_lines.strip())
                        break
            line_batch_index = 0
            for line in vcfFileHandle:
                tot_vars_in_sample += 1
                line = line.strip()
                lines_to_write.append(line)
                line_batch_index += 1
                if line_batch_index == SSMergeCommonUtils.LINE_BATCH_SIZE:
                    self.write_variants_to_cassandra(lines_to_write)
                    line_batch_index = 0
                    lines_to_write = []
            if lines_to_write:
                self.write_variants_to_cassandra(lines_to_write)

            return tot_vars_in_sample

    def get_sample_processed_status(self):
        """
        Check if the variants from a specific sample have already been inserted into the requisite Cassandra tables       
        
        :return: True if variants were inserted for the given sample, False otherwise
        :rtype: bool
        """
        allrows = self.cassandraSession.execute(
            "select proc_status from {0}.{1} where samplename = '{2}' allow filtering;"
            .format(self.keyspaceName, self.sampleInsertLogTableName, self.sampleName),
            timeout=SSMergeCommonUtils.LARGE_QUERY_TIMEOUT_IN_SECS)
        if allrows:
            allrows = iter(allrows)
            first_row = allrows.next()
        else:
            return False
        if first_row.proc_status:
            return first_row.proc_status
        return False

    def update_processing_status(self, proc_status, num_variants):
        """
        Helper function - Update processing status as the sample files are being processed.        
        
        :param proc_status: Processing status ("variants_filtered" or "variants_inserted")
        :type proc_status: str
        :param num_variants: Number of variants inserted for a given sample
        :type num_variants: int
        :return: 
        """
        self.cassandraSession.execute(
            "insert into {0}.{1} (samplename, proc_status, num_variants) values ('{2}', '{3}',{4})"
            .format(self.keyspaceName, self.sampleInsertLogTableName, self.sampleName, proc_status, num_variants),
            timeout=SSMergeCommonUtils.BATCH_WRITE_TIMEOUT_IN_SECS)


class ProcessVariantMatch:
    def __init__(self, *args):

        self.bcf_tools_dir, self.study_name, self.study_vcf_file_name, self.study_files_input_dir, \
            self.variantPositionFileName, self.numTotVariants, self.defaultGenotype, self.missingGenotype, \
            self.cassandraNodeIPs, self.keyspaceName, self.sampleDefaultsTableName, self.variantTableName, \
            self.sampleInsertLogTableName = args

        self.sampleName = SSMergeCommonUtils.getSampleName(self.bcf_tools_dir, self.study_vcf_file_name)
        self.cassandraCluster = Cluster(self.cassandraNodeIPs)
        self.cassandraSession = self.cassandraCluster.connect(self.keyspaceName)
        self.defaultGenotypeInsertPrepStmt = self.cassandraSession.prepare(
            "insert into {0}.{1} (samplename, default_genotype) values (?,?)".format(self.keyspaceName,
                                                                                     self.sampleDefaultsTableName))
        self.variantInsertPrepStmt = self.cassandraSession.prepare(
            "insert into {0}.{1} (chrom,chunk,start_pos,var_uniq_id,ref,alt,samplename, sampleinfoformat, sampleinfo) "
            "values (?,?,?,?,?,?,?,?,?)".format(self.keyspaceName, self.variantTableName))

    def insert_freq_genotype_to_cassandra(self, assumed_frequent_genotype):
        """
        Insert the most frequent genotype, assumed for a given sample, into Cassandra.            
        """
        self.cassandraSession.execute(self.defaultGenotypeInsertPrepStmt.bind
                                      ([self.sampleName, assumed_frequent_genotype]))

    @staticmethod
    def get_variant_info(match_line, sample_name):
        """
        Helper function to get variant information from the match file generated by this program 

        :param match_line: A line from the match file generated by this program
        :type match_line: str
        :param sample_name: Name of the sample from which we are reading the match line
        :type sample_name: str
        """
        chrom_match, var_pos_match, ref, alt, genotype, format_match, sampleinfo_match = match_line.strip().split("\t")
        var_pos_match = long(float(var_pos_match))
        chunk = int(var_pos_match / SSMergeCommonUtils.CHR_POS_CHUNKSIZE)
        variant_id = SSMergeCommonUtils.getUniqueVariantID(sample_name, chrom_match, var_pos_match, ref, alt)
        return chrom_match, var_pos_match, ref, alt, genotype, format_match, sampleinfo_match, chunk, variant_id

    def process_variant_match_file(self, match_output_filename, assumed_frequent_genotype):
        """
        Process a variant match file that is generated for every sample. The variant match file contains only entries 
        at the variant positions (determined from the first pass. see SingleSampleMerge.py) for a given sample.

        :param match_output_filename: Output file that contains the result of the sample file matched against 
        the variant position file
        :type match_output_filename: str
        :param assumed_frequent_genotype: Dominant genotype assumed for the sample
        :type assumed_frequent_genotype: str
        """
        # Insert default genotype assumed for the sample into Cassandra
        self.insert_freq_genotype_to_cassandra(assumed_frequent_genotype)

        # Compare the universal variant list against the matches found
        # at variant positions in the individual sample files
        with gzip.open(self.variantPositionFileName) as varListFile:
            with gzip.open(match_output_filename) as varMatchFile:
                for match_line in varMatchFile:
                    chrom_match, var_pos_match, ref, alt, genotype, \
                        format_match, sampleinfo_match, chunk, variant_id = self.get_variant_info(match_line, 
                                                                                                  self.sampleName)

                    for var_pos_line in varListFile:
                        var_pos_line = var_pos_line.strip()
                        chrom_to_find, var_pos_to_find = var_pos_line.split("\t")
                        var_pos_to_find = long(float(var_pos_to_find))

                        # region Handle an edge case THAT SHOULD NOT HAPPEN
                        # If a variant position from the universal variant position list is farther than
                        # the position in the match file. This can happen if the match file has multiple entries for
                        # the same position. In such cases, forward to a position in the match file which is greater
                        # than or equal to the variant position from the universal variant position list.
                        while var_pos_to_find > var_pos_match:
                            match_line = varMatchFile.readline()
                            if not match_line:
                                return
                            chrom_match, var_pos_match, ref, alt, genotype, \
                                format_match, sampleinfo_match, chunk, variant_id \
                                = self.get_variant_info(match_line, self.sampleName)
                        # endregion

                        # If we assumed the default was a "missing" genotype (ex: ./.) for a sample, insert to Cassandra
                        # if we have evidence otherwise i.e.,a match was found at a variant position in the sample.
                        #
                        # If we assumed the default was a "non-missing" genotype (ex: 0/0) for a sample,
                        # insert to Cassandra if we have evidence otherwise
                        # i.e.,no variant was found at a position in the sample (or)
                        # the variant found at the position had a different genotype than the default genotype (ex: 0/1)
                        if (chrom_to_find, var_pos_to_find) == (chrom_match, var_pos_match):
                            if genotype == assumed_frequent_genotype or alt != '.':
                                break
                            self.cassandraSession.execute(self.variantInsertPrepStmt.bind(
                                [chrom_to_find, chunk, var_pos_to_find, variant_id, ref, alt,
                                 self.sampleName, format_match, sampleinfo_match]))
                            break
                        else:
                            if assumed_frequent_genotype != "./." and assumed_frequent_genotype != ".|.":
                                self.cassandraSession.execute(self.variantInsertPrepStmt.bind(
                                    [chrom_to_find, chunk, var_pos_to_find, variant_id, ref, alt,
                                     self.sampleName, "GT", self.missingGenotype]))

    def get_num_variants_in_sample(self):
        """
        Get number of variant positions recorded for a given sample during the first pass (see SingleSampleMerge.py)

        :return: Number of variant positions recorded for the sample
        :rtype: int
        """
        rows = self.cassandraSession.execute("select num_variants from {0}.{1} where samplename = '{2}';"
                                             .format(self.keyspaceName, self.sampleInsertLogTableName,
                                                     self.sampleName), timeout=None)
        if rows:
            rows = iter(rows)
            return rows.next().num_variants
        else:
            raise Exception("Could not find number of matched variants for the sample: {0}".format(self.sampleName))

    def match_variant_pos_in_study_files(self):
        return_err_msg = None
        chosen_file_to_process = self.study_vcf_file_name
        try:
            sample_name = SSMergeCommonUtils.getSampleName(self.bcf_tools_dir, self.study_vcf_file_name)
            match_output_file_name = sample_name + "_variantmatch.gz"
            err_file_name = sample_name + "_variantmatch.err"
            bcf_variant_match_cmd = "({0}/bin/bcftools query " \
                "-f'%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\t%LINE' -T {1} {3} " \
                "| cut -f1,2,3,4,5,14,15 " \
                "| gzip) " \
                "1> {2} 2> {4}".format(self.bcf_tools_dir, self.variantPositionFileName,
                                       match_output_file_name, chosen_file_to_process, err_file_name)

            os.system(bcf_variant_match_cmd)
            errlines = SSMergeCommonUtils.getErrFileContents(err_file_name)
            if not errlines:
                num_variants_in_sample = self.get_num_variants_in_sample()
                # Default to the missing genotype if there were less than 50% of the total variants in the sample
                if (self.numTotVariants * 1.0 / num_variants_in_sample) > 2:
                    self.process_variant_match_file(match_output_file_name, self.missingGenotype)
                else:
                    self.process_variant_match_file(match_output_file_name, self.defaultGenotype)
            else:
                return_err_msg = "Error in processing file:{0}".format(
                    chosen_file_to_process) + os.linesep + os.linesep.join(errlines)
        except Exception, e:
            return_err_msg = "Error in processing file:{0}".format(
                chosen_file_to_process) + os.linesep + e.message + os.linesep + traceback.format_exc()
        finally:
            self.cassandraCluster.shutdown()
            self.cassandraSession.shutdown()
            if return_err_msg:
                return return_err_msg
            return None


def match_variant_pos_in_study_files(bcf_tools_dir, study_name, study_vcf_file_name, study_files_input_dir,
                                     variant_position_file_name, num_total_variants, default_genotype, missing_genotype,
                                     cassandra_node_ips, keyspace_name, sample_defaults_table_name,
                                     variant_table_name, sample_insert_log_table_name):
    """
    Match the universal set of unique variant positions from across all samples 
    (obtained from first pass, see SingleSampleMerge.py)

    :param bcf_tools_dir: Directory containing the bcftools static binaries
    :type bcf_tools_dir: str
    :param study_name: Study Name
    :type study_name: str
    :param study_vcf_file_name: File Name for the specific sample being matched for variant positions
    :type study_vcf_file_name: str
    :param study_files_input_dir: Input directory where the single sample files reside
    :type study_files_input_dir: str
    :param variant_position_file_name: Full path to file that has all the unique variant positions 
    obtained form first pass
    :type variant_position_file_name: str
    :param num_total_variants: Total number of variants that were inserted into Cassandra during the first pass
    :type num_total_variants: long
    :param default_genotype: Default genotype to use for the sample
    :type default_genotype: str
    :param missing_genotype: Missing genotype to use for the sample
    :type missing_genotype: str
    :param cassandra_node_ips: Set of IP addresses to connect to the Cassandra cluster
    :type cassandra_node_ips: list[str]
    :param keyspace_name: Cassandra key space where the variant table resides
    :type keyspace_name: str
    :param sample_defaults_table_name: Name of the Cassandra table that has the default genotypes per sample
    :type sample_defaults_table_name: str
    :param variant_table_name: Cassandra variant table
    :type variant_table_name: str
    :param sample_insert_log_table_name: Cassandra sample insert log table to determine the 
    number of variants in the sample
    :type sample_insert_log_table_name: str
    :return: Error messages, if any, during the function execution. None otherwise.
    :rtype: str
    """
    proc_var_match = ProcessVariantMatch(bcf_tools_dir, study_name, study_vcf_file_name, study_files_input_dir,
                                         variant_position_file_name, num_total_variants, default_genotype,
                                         missing_genotype, cassandra_node_ips, keyspace_name,
                                         sample_defaults_table_name, variant_table_name, sample_insert_log_table_name)
    return proc_var_match.match_variant_pos_in_study_files()


def process_study_files(bcf_tools_dir, study_name, study_vcf_filename, cassandra_node_ips, keyspace_name,
                        header_table_name, variant_table_name, sample_insert_log_table_name, unique_pos_table_name):
    """
    Process a study file (runs on each core on each Apache Spark node)

    :param study_name: Name of the study (ex: PRJEB21300)
    :type study_name: str
    :param study_vcf_filename: Full path to the study file (ex: /mnt/shared_storage/PRJEB21300/StudySample1.vcf.gz)
    :type study_vcf_filename: str
    :param cassandra_node_ips: Set of "seed" IPs of Cassandra nodes to initialize connection
    :type cassandra_node_ips: list[str]
    :param keyspace_name: Cassandra keyspace (analogous to a database) where variant information should be written to
    :type keyspace_name: str
    :param header_table_name: Name of the Cassandra table that stores variant header information
    :type header_table_name: str
    :param variant_table_name: Name of the Cassandra table that stores variant information
    :type variant_table_name: str
    :param sample_insert_log_table_name: Cassandra table that stores the status of variant insertion from a sample
    :type sample_insert_log_table_name: str
    :param unique_pos_table_name: Name of the Cassandra table that stores unique chromosome + start_pos information
    :type unique_pos_table_name: str
    :param bcf_tools_dir: Directory containing the bcftools static binaries (ex: /mnt/shared_storage/bcftools)
    :type bcf_tools_dir: str
    :return: Tuple - Total variants in sample, Error messages if any 
    """
    ss_merge_obj = SingleSampleMerge(bcf_tools_dir, study_name, study_vcf_filename, cassandra_node_ips, keyspace_name,
                                     header_table_name, variant_table_name,
                                     sample_insert_log_table_name, unique_pos_table_name)
    return ss_merge_obj.process_study_file()


def get_chrom_chunk_set(cassandra_session, keyspace_name, unique_pos_table_name):
    """
    Helper function - Get the set of chromosome, chunk combinations for the study
    
    :param cassandra_session: Cassandra session object
    :param keyspace_name: Cassandra keyspace that holds the table with the unique variant positions 
    :param unique_pos_table_name: Name of the Cassandra table that holds the unique variant positions
    :return: Set of unique chromosome, chunk combinations
    """
    chrom_chunk_set = set()
    # Due to Cassandra's limitations, we cannot filter for in_src = 1. Instead we filter this manually as shown below.
    results = cassandra_session.execute("select distinct chrom, chunk, in_src from {0}.{1};"
                                        .format(keyspace_name, unique_pos_table_name),
                                        timeout=SSMergeCommonUtils.LARGE_QUERY_TIMEOUT_IN_SECS)
    for result in results:
        if result.in_src == 1:
            chrom_chunk_set.add((result.chrom, result.chunk))
    return chrom_chunk_set


def get_cassandra_total_variants_count(cassandra_node_ips, keyspace_name, variant_table_name, chrom, chunk):
    """
    Helper function - Get total variants count from Cassandra

    :param cassandra_node_ips: IP addresses of a few Cassandra nodes 
    :param keyspace_name: Cassandra keyspace that holds the table with the variants 
    :param variant_table_name: Name of the Cassandra variant table
    :param chrom: Chromosome to look for
    :param chunk: Chunk to look for
    :return: Number of records in Cassandra variant table for the given chromosome and chunk 
    """
    cluster = Cluster(cassandra_node_ips)
    session = cluster.connect()
    total_variants_count = None
    results = session.execute("select count(*) as num_variants from {0}.{1} where chrom = '{2}' and chunk = {3};"
                              .format(keyspace_name, variant_table_name, chrom, chunk),
                              timeout=SSMergeCommonUtils.LARGE_QUERY_TIMEOUT_IN_SECS)
    if results:
        results = iter(results)
        total_variants_count = results.next().num_variants
    session.shutdown()
    cluster.shutdown()
    if total_variants_count is not None:
        return total_variants_count
    return -1


def get_cassandra_distinct_var_pos_count(cassandra_node_ips, keyspace_name, unique_pos_table_name, chrom, chunk):
    """
    Helper function - Get total number of distinct variant positions from Cassandra
    
    :param cassandra_node_ips: IP addresses of a few Cassandra nodes
    :param keyspace_name: Cassandra keyspace that holds the table with the unique variant positions
    :param unique_pos_table_name: Name of the Cassandra table that holds the unique variant positions
    :param chrom: Chromosome to look for
    :param chunk: Chunk to look for
    :return: Number of distinct variant positions stored in Cassandra for the given chromosome and chunk
    """
    cluster = Cluster(cassandra_node_ips)
    session = cluster.connect()
    tot_variants_count = None
    results = session.execute("select count(*) as num_variants from {0}.{1} where chrom = '{2}' "
                              "and chunk = {3} and in_src = 1;"
                              .format(keyspace_name, unique_pos_table_name, chrom, chunk),
                              timeout=SSMergeCommonUtils.LARGE_QUERY_TIMEOUT_IN_SECS)
    if results:
        results = iter(results)
        tot_variants_count = results.next().num_variants
    session.shutdown()
    cluster.shutdown()
    if tot_variants_count is not None:
        return tot_variants_count
    return -1


def validate_record_count(cassandra_node_ips, keyspace_name, local_session_var, results, spark_context,
                          study_name, variant_table_name, study_info_table_name, unique_pos_table_name):
    atleast_one_error = False
    total_variants_in_sample_files = 0
    for (totVarsInSample, errMsg) in results:
        if errMsg:
            print(errMsg)
            atleast_one_error = True
        else:
            # Count the number of variants in the sample VCF files
            total_variants_in_sample_files += totVarsInSample
    if not atleast_one_error:
        print("Variant insertion into Cassandra complete!")
        print("Beginning record count validation.........")

        # Count the number of variants from the variants table
        chrom_chunk_set = get_chrom_chunk_set(local_session_var, keyspace_name, unique_pos_table_name)
        chrom_chunk_set_rdd = spark_context.parallelize(chrom_chunk_set)
        total_variants_per_chrom_chunk = chrom_chunk_set_rdd.map(lambda entry:
                                                                 get_cassandra_total_variants_count
                                                                 (cassandra_node_ips, keyspace_name,
                                                                  variant_table_name, entry[0], entry[1])).collect()
        total_variants_in_cassandra = sum(total_variants_per_chrom_chunk)
        total_distinct_variant_positions_per_chrom_chunk = chrom_chunk_set_rdd.map(lambda entry:
                                                                                   get_cassandra_distinct_var_pos_count
                                                                                   (cassandra_node_ips, keyspace_name,
                                                                                    unique_pos_table_name, entry[0],
                                                                                    entry[1])).collect()
        total_distinct_variant_positions_in_cassandra = sum(total_distinct_variant_positions_per_chrom_chunk)
        # Register counts of variants
        local_session_var.execute("insert into {0}.{1} (study_name, tot_num_variants, distinct_num_variant_pos) "
                                  "values ('{2}',{3}, {4})"
                                  .format(keyspace_name, study_info_table_name, study_name, total_variants_in_cassandra,
                                          total_distinct_variant_positions_in_cassandra))

        print("******************************************************************************************************")
        print("Number of variants as counted from source files: {0}".format(total_variants_in_sample_files))
        print("Number of variant inserts recorded as retrieved from Cassandra: {0}".format(total_variants_in_cassandra))
        print("Number of distinct variant positions as retrieved from Cassandra: {0}"
              .format(total_distinct_variant_positions_in_cassandra))
        print("******************************************************************************************************")

        return total_variants_in_cassandra, total_variants_in_sample_files == total_variants_in_cassandra


def init_spark_env(cassandra_node_ips):
    conf = SparkConf().setMaster("spark://{0}:7077".format(SSMergeCommonUtils.get_ip_address())) \
        .setAppName("SingleSampleVCFMerge").set("spark.cassandra.connection.host", cassandra_node_ips[0]) \
        .set("spark.cassandra.read.timeout_ms", SSMergeCommonUtils.SPARK_CASSANDRA_READ_TIMEOUT_IN_MS) \
        .set("spark.cassandra.connection.timeout_ms", SSMergeCommonUtils.SPARK_CASSANDRA_READ_TIMEOUT_IN_MS)
    sc = SparkContext(conf=conf)
    sc.setLogLevel("WARN")
    return conf, sc


def init_cassandra_env(cassandra_node_ips, study_name):
    keyspace_name = "variant_ksp"
    variant_table_name = "variants_{0}".format(study_name.lower())
    header_table_name = "headers_{0}".format(study_name.lower())
    sample_insert_log_table_name = "sample_insert_log_{0}".format(study_name.lower())
    study_info_table_name = "study_info_{0}".format(study_name.lower())
    unique_pos_table_name = "uniq_pos_{0}".format(study_name.lower())
    sample_defaults_table_name = "sample_defaults_{0}".format(study_name.lower())
    local_cluster_var = Cluster(cassandra_node_ips)
    local_session_var = local_cluster_var.connect()
    local_session_var.execute("create keyspace if not exists {0}".format(keyspace_name) +
                              " with replication = {'class': 'SimpleStrategy', 'replication_factor' : 1};",
                              timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(samplename varchar, header varchar, "
                              "primary key(samplename));"
                              .format(keyspace_name, header_table_name),
                              timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(chrom varchar, chunk int, start_pos bigint, ref varchar, alt varchar, qual varchar, "
                              "filter varchar, info varchar, sampleinfoformat varchar, "
                              "sampleinfo  varchar, var_id varchar, var_uniq_id varchar, sampleName varchar,  "
                              "primary key((chrom, chunk), start_pos, var_uniq_id, samplename));"
                              .format(keyspace_name, variant_table_name),
                              timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(samplename varchar, proc_status varchar, num_variants bigint, primary key(samplename));"
                              .format(keyspace_name, sample_insert_log_table_name),
                              timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(study_name varchar, tot_num_variants bigint, distinct_num_variant_pos bigint, "
                              "primary key(study_name));"
                              .format(keyspace_name, study_info_table_name),
                              timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(chrom varchar, chunk int, in_src int, start_pos bigint, "
                              "primary key((chrom,chunk,in_src), start_pos));"
                              .format(keyspace_name, unique_pos_table_name),
                              timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} (samplename varchar, default_genotype varchar, "
                              "primary key(samplename));"
                              .format(keyspace_name, sample_defaults_table_name),
                              timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    return keyspace_name, local_cluster_var, local_session_var, study_info_table_name, variant_table_name, \
        header_table_name, sample_insert_log_table_name, unique_pos_table_name, sample_defaults_table_name


# Helper function to write the unique variant positions from the first pass into a single file
def gen_var_pos_file(spark_context, study_files_input_dir, keyspace_name, unique_pos_table_name):
    variant_position_file_name = study_files_input_dir + os.path.sep + "unique_variant_positions.gz"
    sql = SQLContext(spark_context)
    variants = sql.read.format("org.apache.spark.sql.cassandra"). \
        load(keyspace=keyspace_name, table=unique_pos_table_name)
    variants.registerTempTable("variantsTable")
    result_df = sql.sql("select chrom,start_pos from variantsTable order by 1,2")
    iterator = result_df.toLocalIterator()
    with gzip.open(variant_position_file_name, "wb") as variantPositionFileHandle:
        for result in iterator:
            # For some reason, if you don't store the output of the write, every write result is echoed to the screen!!!
            discard_output = variantPositionFileHandle.write(result["chrom"] + "\t" + str(result["start_pos"]) +
                                                             os.linesep)
    return variant_position_file_name


def main_func():
    if len(sys.argv) != 8:
        print("Usage: SingleSampleMerge.py <Study PRJ ID> <Default Genotype> <Missing Genotype> "
              "<Full Path to study files> <Cassandra node IP1> <Cassandra node IP2> <BCF Tools Directory>")
        sys.exit(1)

    start_time = datetime.datetime.now()
    # Parse arguments
    study_name, default_genotype, missing_genotype, study_files_input_dir = sys.argv[1:5]
    cassandra_node_ips, bcf_tools_dir = [sys.argv[5], sys.argv[6]], sys.argv[7]

    # Get the list of study files
    study_vcf_file_names = SSMergeCommonUtils.getVCFFileNames(study_files_input_dir)

    # region Create keyspace and the required tables for variants, headers, sample insertion log, study info and
    # unique chromosome position tables in Cassandra
    keyspace_name, local_cluster_var, local_session_var, study_info_table_name, variant_table_name, header_table_name, \
        sample_insert_log_table_name, unique_pos_table_name, sample_defaults_table_name = \
        init_cassandra_env(cassandra_node_ips, study_name)
    # endregion

    # region Initialize Apache Spark
    spark_conf, spark_context = init_spark_env(cassandra_node_ips)
    # endregion

    # region Process study files in parallel
    num_partitions = len(study_vcf_file_names)
    study_indiv_rdd = spark_context.parallelize(study_vcf_file_names, num_partitions)
    results = study_indiv_rdd.map(lambda entry: process_study_files(bcf_tools_dir, study_name, entry,
                                                                    cassandra_node_ips, keyspace_name,
                                                                    header_table_name, variant_table_name,
                                                                    sample_insert_log_table_name,
                                                                    unique_pos_table_name)).collect()
    # endregion

    # region Validate Cassandra record counts against record counts in the source file
    total_variants_in_cassandra, record_count_matched = validate_record_count(cassandra_node_ips, keyspace_name,
                                                                              local_session_var, results, spark_context,
                                                                              study_name, variant_table_name,
                                                                              study_info_table_name,
                                                                              unique_pos_table_name)
    if not record_count_matched:
        print("ERROR: Record counts did not match between the source study files and the Cassandra variant table")
    else:
        # Generate a position file with the unique set of variant positions determined from the first pass
        variant_position_file_name = gen_var_pos_file(spark_context, study_files_input_dir, keyspace_name,
                                                      unique_pos_table_name)

        # Scan each sample file against the unique set of variant positions
        # to insert sample genotypes at these positions
        num_partitions = len(study_vcf_file_names)
        study_indiv_rdd = spark_context.parallelize(study_vcf_file_names, num_partitions)
        process_results = study_indiv_rdd.map(lambda study_vcf_file_name:
                                              match_variant_pos_in_study_files(bcf_tools_dir, study_name,
                                                                               study_vcf_file_name,
                                                                               study_files_input_dir,
                                                                               variant_position_file_name,
                                                                               total_variants_in_cassandra,
                                                                               default_genotype, missing_genotype,
                                                                               cassandra_node_ips, keyspace_name,
                                                                               sample_defaults_table_name,
                                                                               variant_table_name,
                                                                               sample_insert_log_table_name)).collect()
        for result in process_results:
            if result:
                print(result)
    # endregion

    # region Print elapsed time
    end_time = datetime.datetime.now()
    print("Time taken to process the study {0}: {1}".format(study_name, end_time - start_time))
    # endregion

    # region Shutdown Cassandra and Apache Spark
    local_session_var.shutdown()
    local_cluster_var.shutdown()
    spark_context.stop()
    # endregion


if __name__ == '__main__':
    main_func()
