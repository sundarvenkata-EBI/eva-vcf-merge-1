# Copyright 2014-2017 EMBL - European Bioinformatics Institute
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

import os, sys, gzip, SSMergeCommonUtils, traceback, datetime
from cassandra.cluster import Cluster
from cassandra.query import BatchStatement, BatchType
from pyspark import SparkConf, SparkContext, SQLContext

class SingleSampleMerge:
    def __init__(self, *args):
        self.bcfToolsDir, self.studyName, self.studyVCFFileName, self.cassandraNodeIPs, self.keyspaceName, \
        self.headerTableName, self.variantTableName, self.sampleInsertLogTableName, self.uniquePosTableName = args

        self.cassandraCluster = Cluster(self.cassandraNodeIPs)
        self.cassandraSession = self.cassandraCluster.connect(self.keyspaceName)
        self.sampleName = SSMergeCommonUtils.getSampleName(self.bcfToolsDir, self.studyVCFFileName)
        self.filterVCFFileName = os.path.dirname(self.studyVCFFileName) + os.path.sep \
                                 + "{0}_filtervcf.gz".format(self.sampleName)
        self.sample_proc_status = self.getSampleProcessedStatus()

        self.variantInsertStmt = self.cassandraSession.prepare(
            "insert into {0}.{1} "
            "(chrom,chunk,start_pos,ref,alt,qual,filter,info,sampleinfoformat,sampleinfo,var_id,var_uniq_id,sampleName) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)".format(self.keyspaceName, self.variantTableName))
        self.headerInsertStmt = self.cassandraSession.prepare("insert into {0}.{1} (samplename, header) VALUES (?,?)"
                                                              .format(self.keyspaceName, self.headerTableName))
        self.uniquePosInsertStmt = self.cassandraSession.prepare(
            "insert into {0}.{1} (chrom, chunk, in_src, start_pos) VALUES (?,?,?,?)"
                .format(self.keyspaceName, self.uniquePosTableName))


    def processStudyFile(self):
        """
        Process the study file that this specific SingleSampleMerge object has been given
        
        :return: Tuple - Total variants in sample, Error messages if any  
        """
        errFileContents, returnErrMsg, cassandraCluster, cassandraSession = None, None, None, None
        totVarsInSample = 0

        if self.sample_proc_status == "variants_inserted":
            # return 0,0,None
            return (-1, "Already processed sample: {0} for study: {1}"
                                                                   .format(self.sampleName, self.studyName))

        try:
            # Get the sample prefix and the base directory of the file that we are operating on
            baseDir = os.path.dirname(self.studyVCFFileName)

            filterErrFileName = baseDir + os.path.sep + "{0}_filtering_err.txt".format(self.sampleName)

            # region Run bcftools filter command to filter out monomorphic references
            if self.sample_proc_status != "variants_filtered":

                if os.path.isfile(self.filterVCFFileName): os.remove(self.filterVCFFileName)
                os.system("""({0}/bin/bcftools filter -e ALT=\\'.\\' {1} | gzip) 1>{2} 2>{3}"""
                          .format(self.bcfToolsDir, self.studyVCFFileName, self.filterVCFFileName, filterErrFileName))
                errFileContents = SSMergeCommonUtils.getErrFileContents(filterErrFileName)
                if errFileContents:
                    filterCommandResult = -1
                else:
                    filterCommandResult = 0
                    self.update_processing_status("variants_filtered", 0)
            else:
                filterCommandResult = 0
            # endregion

            # region Insert the variants into Cassandra from the filtered VCF files obtained from above
            if filterCommandResult != 0:
                returnErrMsg = "Failed to process {0} due to error: {1}".format(self.studyVCFFileName, errFileContents)
            else:
                totVarsInSample = self.cassandraInsert()
                os.system("echo {0} > {1}_filtered_variant_count.txt".format(str(totVarsInSample), self.sampleName))
            # endregion

            # region Update status to indicate that all variants have been inserted into Cassandra
            if not returnErrMsg and totVarsInSample != 0:
                self.update_processing_status("variants_inserted", totVarsInSample)
            # endregion

        except Exception, ex:
            returnErrMsg = "Error in processing file:{0}".format(self.studyVCFFileName) + os.linesep + ex.message \
                           + os.linesep + traceback.format_exc()
        finally:
            if cassandraCluster is not None and cassandraSession is not None:
                try:
                    cassandraSession.shutdown()
                    cassandraCluster.shutdown()
                except Exception, ex:
                    print(ex.message)
            if returnErrMsg: return -1, returnErrMsg
            return totVarsInSample, None

    def writeVariantToCassandra(self, linesToWrite):
        """
        Write a set of variant lines to Cassandra
        
        :param linesToWrite: Set of variant lines to write from a given sample
        :type linesToWrite: list[str]        
        :return: None
        """
        batch = BatchStatement(BatchType.UNLOGGED)

        for line in linesToWrite:
            lineComps = [x.strip() for x in line.split("\t")]
            chromosome, position, rsID, ref, alt, qual, qualFilter, info = lineComps[:8]
            sampleInfoFormat, sampleInfo = None, None
            if len(lineComps) > 8: sampleInfoFormat, sampleInfo = lineComps[8:]

            # Chromosome + Chunk combination determines which node a variant is written to
            position = long(float(position))
            chunk = int(position / SSMergeCommonUtils.CHR_POS_CHUNKSIZE)

            # Insert the chromosome position that was scanned
            batch.add(self.uniquePosInsertStmt.bind([chromosome, chunk, 1, position]))
            # Also, insert the normalized chromosome positions for each alternate allele
            for altAllele in alt.split(","):
                revised_start_pos = SSMergeCommonUtils.getNormalizedStartPos(chromosome, position, ref, altAllele)
                if revised_start_pos != position:
                    revisedChunk = int(revised_start_pos / SSMergeCommonUtils.CHR_POS_CHUNKSIZE)
                    batch.add(self.uniquePosInsertStmt.bind([chromosome, revisedChunk, 0, revised_start_pos]))

            # Unique Variant ID is chromosome_12-digit zero padded start_MD5(REF_ALT)_20-character samplename
            variantID = SSMergeCommonUtils.getUniqueVariantID(self.sampleName, chromosome, position, ref, alt)
            # Insert the variant information
            boundStmt = self.variantInsertStmt.bind([chromosome, chunk, position, ref, alt, qual, qualFilter, info,
                                                     sampleInfoFormat, sampleInfo, rsID, variantID, self.sampleName])
            batch.add(boundStmt)

        # Batch execute the insert statements added in the loop above
        self.cassandraSession.execute(batch, timeout=SSMergeCommonUtils.BATCH_WRITE_TIMEOUT_IN_SECS)



    def writeHeaderToCassandra(self, headerLines):
        """
        Write a set of header lines to Cassandra
    
        :param headerLines: Set of header lines to write from a given sample
        :type headerLines: str        
        :return: None
        """
        self.cassandraSession.execute(self.headerInsertStmt.bind([self.sampleName, headerLines]),
                                      timeout=SSMergeCommonUtils.BATCH_WRITE_TIMEOUT_IN_SECS)

    def cassandraInsert(self):
        """
        Insert the contents of a VCF file: variants and headers into relevant Cassandra tables        
        
        :return: Total number of variants that were scanned from the VCF file
        :rtype: long
        """
        lineBatchSize = 50  # Write variants in batches
        totVarsInSample = 0
        headerLines, sampleName = "", ""
        linesToWrite = []

        with gzip.open(self.filterVCFFileName, 'rb') as vcfFileHandle:
            for line in vcfFileHandle:
                line = line.strip()
                if line.startswith("#"):
                    headerLines += (line + os.linesep)
                    if line.startswith("#CHROM"):
                        self.writeHeaderToCassandra(headerLines.strip())
                        break
            lineBatchIndex = 0
            for line in vcfFileHandle:
                totVarsInSample += 1
                line = line.strip()
                linesToWrite.append(line)
                lineBatchIndex += 1
                if lineBatchIndex == lineBatchSize:
                    self.writeVariantToCassandra(linesToWrite)
                    lineBatchIndex = 0
                    linesToWrite = []
            if linesToWrite:
                self.writeVariantToCassandra(linesToWrite)

            return totVarsInSample

    def getSampleProcessedStatus(self):
        """
        Check if the variants from a specific sample have already been inserted into the requisite Cassandra tables       
        
        :return: True if variants were inserted for the given sample, False otherwise
        :rtype: bool
        """
        allrows = self.cassandraSession.execute("select proc_status from {0}.{1} where samplename = '{2}' allow filtering;"
                                                .format(self.keyspaceName, self.sampleInsertLogTableName,
                                                        self.sampleName),
                                                timeout=SSMergeCommonUtils.LARGE_QUERY_TIMEOUT_IN_SECS)
        if allrows:
            allrows = iter(allrows)
            firstRow = allrows.next()
        else:
            return False
        if firstRow.proc_status: return firstRow.proc_status
        return ""

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
                .format(self.keyspaceName, self.sampleInsertLogTableName, self.sampleName,proc_status, num_variants),
            timeout=SSMergeCommonUtils.BATCH_WRITE_TIMEOUT_IN_SECS)


class ProcessVariantMatch:
    def __init__(self, *args):

        self.bcfToolsDir, self.studyName, self.studyVCFFileName, self.studyFilesInputDir, self.variantPositionFileName, \
        self.numTotVariants, self.defaultGenotype, self.missingGenotype, self.cassandraNodeIPs, self.keyspaceName, \
        self.sampleDefaultsTableName, self.variantTableName, self.sampleInsertLogTableName = args

        self.sampleName = SSMergeCommonUtils.getSampleName(self.bcfToolsDir, self.studyVCFFileName)
        self.cassandraCluster = Cluster(self.cassandraNodeIPs)
        self.cassandraSession = self.cassandraCluster.connect(self.keyspaceName)
        self.defaultGenotypeInsertPrepStmt = self.cassandraSession.prepare(
            "insert into {0}.{1} (samplename, default_genotype) values (?,?)".format(self.keyspaceName,
                                                                                     self.sampleDefaultsTableName))
        self.variantInsertPrepStmt = self.cassandraSession.prepare(
            "insert into {0}.{1} (chrom,chunk,start_pos,var_uniq_id,ref,alt,samplename, sampleinfoformat, sampleinfo) "
            "values (?,?,?,?,?,?,?,?,?)".format(self.keyspaceName, self.variantTableName))

    def insertFreqGenotypeToCassandra(self, assumedFrequentGenotype):
        """
        Insert the most frequent genotype, assumed for a given sample, into Cassandra.            
        """
        self.cassandraSession.execute(self.defaultGenotypeInsertPrepStmt.bind
                                      ([self.sampleName, assumedFrequentGenotype]))

    @staticmethod
    def getVariantInfo(matchLine, sampleName):
        """
        Helper function to get variant information from the match file generated by this program 

        :param matchLine: A line from the match file generated by this program
        :type matchLine: str
        :param sampleName: Name of the sample from which we are reading the match line
        :type sampleName: str
        """
        chromMatch, varPosMatch, ref, alt, genotype, formatMatch, sampleinfoMatch = matchLine.strip().split("\t")
        varPosMatch = long(float(varPosMatch))
        chunk = int(varPosMatch / SSMergeCommonUtils.CHR_POS_CHUNKSIZE)
        variantID = SSMergeCommonUtils.getUniqueVariantID(sampleName, chromMatch, varPosMatch, ref, alt)
        return chromMatch, varPosMatch, ref, alt, genotype, formatMatch, sampleinfoMatch, chunk, variantID

    def processVariantMatchFile(self, matchOutputFileName, assumedFrequentGenotype):
        """
        Process a variant match file that is generated for every sample. The variant match file contains only entries 
        at the variant positions (determined from the first pass. see SingleSampleMerge.py) for a given sample.

        :param matchOutputFileName: Output file that contains the result of the sample file matched against the variant
                                    position file
        :type matchOutputFileName: str
        :param assumedFrequentGenotype: Dominant genotype assumed for the sample
        :type assumedFrequentGenotype: str
        """
        # Insert default genotype assumed for the sample into Cassandra
        self.insertFreqGenotypeToCassandra(assumedFrequentGenotype)

        # Compare the universal variant list against the matches found at variant positions in the individual sample files
        with gzip.open(self.variantPositionFileName, "rb") as varListFile:
            with gzip.open(matchOutputFileName, "rb") as varMatchFile:
                for matchLine in varMatchFile:
                    chromMatch, varPosMatch, ref, alt, genotype, \
                    formatMatch, sampleinfoMatch, chunk, variantID = self.getVariantInfo(matchLine, self.sampleName)

                    for varPosLine in varListFile:
                        varPosLine = varPosLine.strip()
                        chromToFind, varPosToFind = varPosLine.split("\t")
                        varPosToFind = long(float(varPosToFind))

                        # region Handle an edge case THAT SHOULD NOT HAPPEN
                        # If a variant position from the universal variant position list is farther than
                        # the position in the match file. This can happen if the match file has multiple entries for
                        # the same position. In such cases, forward to a position in the match file which is greater than
                        # or equal to the variant position from the universal variant position list.
                        while varPosToFind > varPosMatch:
                            matchLine = varMatchFile.readline()
                            if not matchLine: return
                            chromMatch, varPosMatch, ref, alt, genotype, \
                            formatMatch, sampleinfoMatch, chunk, variantID = self.getVariantInfo(matchLine,
                                                                                                 self.sampleName)
                        # endregion

                        # If we assumed the default was a "missing" genotype (ex: ./.) for a sample, insert to Cassandra
                        # if we have evidence otherwise i.e.,a match was found at a variant position in the sample.
                        #
                        # If we assumed the default was a "non-missing" genotype (ex: 0/0) for a sample, insert to Cassandra
                        # if we have evidence otherwise i.e.,no variant was found at a position in the sample (or)
                        # the variant found at the position had a different genotype than the default genotype (ex: 0/1)
                        if (chromToFind, varPosToFind) == (chromMatch, varPosMatch):
                            if genotype == assumedFrequentGenotype or alt != '.': break
                            self.cassandraSession.execute(self.variantInsertPrepStmt.bind(
                                [chromToFind, chunk, varPosToFind, variantID, ref, alt,
                                 self.sampleName, formatMatch, sampleinfoMatch]))
                            break
                        else:
                            if assumedFrequentGenotype != "./." and assumedFrequentGenotype != ".|.":
                                self.cassandraSession.execute(self.variantInsertPrepStmt.bind(
                                    [chromToFind, chunk, varPosToFind, variantID, ref, alt,
                                     self.sampleName, "GT", self.missingGenotype]))

    def getNumVariantsInSample(self):
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

    def matchVariantPosInStudyFiles(self):
        returnErrMsg = None
        chosenFileToProcess = self.studyVCFFileName
        try:
            sampleName = SSMergeCommonUtils.getSampleName(self.bcfToolsDir, self.studyVCFFileName)
            matchOutputFileName = sampleName + "_variantmatch.gz"
            errFileName = sampleName + "_variantmatch.err"
            bcfVariantMatchCmd = "({0}/bin/bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\t%LINE' -T {1} {3} " \
                                 "| cut -f1,2,3,4,5,14,15 " \
                                 "| gzip) " \
                                 "1> {2} 2> {4}".format(
                self.bcfToolsDir, self.variantPositionFileName, matchOutputFileName, chosenFileToProcess, errFileName)

            os.system(bcfVariantMatchCmd)
            errlines = SSMergeCommonUtils.getErrFileContents(errFileName)
            if not errlines:
                numVariantsInSample = self.getNumVariantsInSample()
                # Default to the missing genotype if there were less than 50% of the total variants in the sample
                if (self.numTotVariants * 1.0 / numVariantsInSample) > 2:
                    self.processVariantMatchFile(matchOutputFileName, self.missingGenotype)
                else:
                    self.processVariantMatchFile(matchOutputFileName, self.defaultGenotype)
            else:
                returnErrMsg = "Error in processing file:{0}".format(
                    chosenFileToProcess) + os.linesep + os.linesep.join(errlines)
        except Exception, e:
            returnErrMsg = "Error in processing file:{0}".format(
                chosenFileToProcess) + os.linesep + e.message + os.linesep + traceback.format_exc()
        finally:
            self.cassandraCluster.shutdown()
            self.cassandraSession.shutdown()
            if returnErrMsg: return returnErrMsg
            return None


def matchVariantPosInStudyFiles(bcfToolsDir, studyName, studyVCFFileName, studyFilesInputDir, variantPositionFileName,
                                numTotVariants, defaultGenotype, missingGenotype, cassandraNodeIPs, keyspaceName,
                                sampleDefaultsTableName, variantTableName, sampleInsertLogTableName):
    """
    Match the universal set of unique variant positions from across all samples 
    (obtained from first pass, see SingleSampleMerge.py)

    :param bcfToolsDir: Directory containing the bcftools static binaries
    :type bcfToolsDir: str
    :param studyName: Study Name
    :type studyName: str
    :param studyVCFFileName: File Name for the specific sample being matched for variant positions
    :type studyVCFFileName: str
    :param studyFilesInputDir: Input directory where the single sample files reside
    :type studyFilesInputDir: str
    :param variantPositionFileName: Full path to file that has all the unique variant positions obtained form first pass
    :type variantPositionFileName: str
    :param numTotVariants: Total number of variants that were inserted into Cassandra during the first pass
    :type numTotVariants: long
    :param defaultGenotype: Default genotype to use for the sample
    :type defaultGenotype: str
    :param missingGenotype: Missing genotype to use for the sample
    :type missingGenotype: str
    :param cassandraNodeIPs: Set of IP addresses to connect to the Cassandra cluster
    :type cassandraNodeIPs: list[str]
    :param keyspaceName: Cassandra key space where the variant table resides
    :type keyspaceName: str
    :param sampleDefaultsTableName: Name of the Cassandra table that has the default genotypes per sample
    :type sampleDefaultsTableName: str
    :param variantTableName: Cassandra variant table
    :type variantTableName: str
    :param sampleInsertLogTableName: Cassandra sample insert log table to determine the number of variants in the sample
    :type sampleInsertLogTableName: str
    :return: Error messages, if any, during the function execution. None otherwise.
    :rtype: str
    """
    procVarMatch = ProcessVariantMatch(bcfToolsDir, studyName, studyVCFFileName, studyFilesInputDir,
                                       variantPositionFileName,
                                       numTotVariants, defaultGenotype, missingGenotype, cassandraNodeIPs, keyspaceName,
                                       sampleDefaultsTableName, variantTableName, sampleInsertLogTableName)
    return procVarMatch.matchVariantPosInStudyFiles()


def processStudyFiles(bcfToolsDir, studyName, studyVCFFileName, cassandraNodeIPs, keyspaceName, headerTableName,
                      variantTableName, sampleInsertLogTableName, uniquePosTableName):
    """
    Process a study file (runs on each core on each Apache Spark node)

    :param studyName: Name of the study (ex: PRJEB21300)
    :type studyName: str
    :param studyVCFFileName: Full path to the study file (ex: /mnt/shared_storage/PRJEB21300/StudySample1.vcf.gz)
    :type studyVCFFileName: str
    :param cassandraNodeIPs: Set of "seed" IPs of Cassandra nodes to initialize connection
    :type cassandraNodeIPs: list[str]
    :param keyspaceName: Cassandra key space (roughly analogous to a database) where the variant information should be written to
    :type keyspaceName: str
    :param headerTableName: Name of the Cassandra table that stores variant header information
    :type headerTableName: str
    :param variantTableName: Name of the Cassandra table that stores variant information
    :type variantTableName: str
    :param sampleInsertLogTableName: Name of the Cassandra table that stores the status of variant insertion from a sample
    :type sampleInsertLogTableName: str
    :param uniquePosTableName: Name of the Cassandra table that stores unique chromosome + start_pos information
    :type uniquePosTableName: str
    :param bcfToolsDir: Directory containing the bcftools static binaries (ex: /mnt/shared_storage/bcftools)
    :type bcfToolsDir: str
    :return: Tuple - Total variants in sample, Error messages if any 
    """
    ssMergeObj = SingleSampleMerge(bcfToolsDir, studyName, studyVCFFileName, cassandraNodeIPs, keyspaceName, headerTableName,
                                   variantTableName, sampleInsertLogTableName, uniquePosTableName)
    return ssMergeObj.processStudyFile()

def getChromChunkSet(cassandraSession, keyspaceName, uniquePosTableName):
    """
    Helper function - Get the set of chromosome, chunk combinations for the study
    
    :param cassandraSession: Cassandra session object
    :param keyspaceName: Cassandra keyspace that holds the table with the unique variant positions 
    :param uniquePosTableName: Name of the Cassandra table that holds the unique variant positions
    :return: Set of unique chromosome, chunk combinations
    """
    chromChunkSet = set()
    # Due to Cassandra's limitations, we cannot filter for in_src = 1. Instead we filter this manually as shown below.
    results = cassandraSession.execute("select distinct chrom, chunk, in_src from {0}.{1};"
                                       .format(keyspaceName, uniquePosTableName),
                                       timeout=SSMergeCommonUtils.LARGE_QUERY_TIMEOUT_IN_SECS)
    for result in results:
        if result.in_src == 1:
            chromChunkSet.add((result.chrom, result.chunk))
    return chromChunkSet


def getCassandraTotVarsCount(cassandraNodeIPs, keyspaceName, variantTableName, chrom, chunk):
    """
    Helper function - Get total variants count from Cassandra

    :param cassandraNodeIPs: IP addresses of a few Cassandra nodes 
    :param keyspaceName: Cassandra keyspace that holds the table with the variants 
    :param variantTableName: Name of the Cassandra variant table
    :param chrom: Chromosome to look for
    :param chunk: Chunk to look for
    :return: Number of records in Cassandra variant table for the given chromosome and chunk 
    """
    cluster = Cluster(cassandraNodeIPs)
    session = cluster.connect()
    totVarsCount = None
    results = session.execute ("select count(*) as num_variants from {0}.{1} where chrom = '{2}' and chunk = {3};"\
        .format(keyspaceName, variantTableName, chrom, chunk), timeout=SSMergeCommonUtils.LARGE_QUERY_TIMEOUT_IN_SECS)
    if results:
        results = iter(results)
        totVarsCount = results.next().num_variants
    session.shutdown()
    cluster.shutdown()
    if totVarsCount is not None: return totVarsCount
    return -1

def getCassandraDistinctVarPosCount(cassandraNodeIPs, keyspaceName, uniquePosTableName, chrom, chunk):
    """
    Helper function - Get total number of distinct variant positions from Cassandra
    
    :param cassandraNodeIPs: IP addresses of a few Cassandra nodes
    :param keyspaceName: Cassandra keyspace that holds the table with the unique variant positions
    :param uniquePosTableName: Name of the Cassandra table that holds the unique variant positions
    :param chrom: Chromosome to look for
    :param chunk: Chunk to look for
    :return: Number of distinct variant positions stored in Cassandra for the given chromosome and chunk
    """
    cluster = Cluster(cassandraNodeIPs)
    session = cluster.connect()
    totVarsCount = None
    results = session.execute("select count(*) as num_variants from {0}.{1} where chrom = '{2}' "
                              "and chunk = {3} and in_src = 1;" \
                              .format(keyspaceName, uniquePosTableName, chrom, chunk),
                              timeout=SSMergeCommonUtils.LARGE_QUERY_TIMEOUT_IN_SECS)
    if results:
        results = iter(results)
        totVarsCount = results.next().num_variants
    session.shutdown()
    cluster.shutdown()
    if totVarsCount is not None: return totVarsCount
    return -1

def validateRecordCount(cassandraNodeIPs, keyspaceName, local_session_var, results, sparkContext,
                        studyName, variantTableName, studyInfoTableName, uniquePosTableName):
    atleastOneError = False
    totVarsInSampleFiles = 0
    for (totVarsInSample, errMsg) in results:
        if errMsg:
            print(errMsg)
            atleastOneError = True
        else:
            # Count the number of variants in the sample VCF files
            totVarsInSampleFiles += totVarsInSample
    if not atleastOneError:
        print("Variant insertion into Cassandra complete!")
        print("Beginning record count validation.........")

        # Count the number of variants from the variants table
        chromChunkSet = getChromChunkSet(local_session_var, keyspaceName, uniquePosTableName)
        chromChunkSetRDD = sparkContext.parallelize(chromChunkSet)
        totVarsCountList = chromChunkSetRDD.map(lambda entry: getCassandraTotVarsCount(cassandraNodeIPs, keyspaceName,
                                                                                       variantTableName, entry[0],
                                                                                       entry[1])).collect()
        totVarsInCassandra = sum(totVarsCountList)
        totDistinctVarPosCountList = chromChunkSetRDD.map(
            lambda entry: getCassandraDistinctVarPosCount(cassandraNodeIPs, keyspaceName, uniquePosTableName, entry[0],
                                                          entry[1])).collect()
        totDistinctVarPosInCassandra = sum(totDistinctVarPosCountList)
        # Register counts of variants
        local_session_var.execute("insert into {0}.{1} (studyname, tot_num_variants, distinct_num_variant_pos) "
                                  "values ('{2}',{3}, {4})"
                                  .format(keyspaceName, studyInfoTableName, studyName, totVarsInCassandra,
                                          totDistinctVarPosInCassandra))

        print("******************************************************************************************************")
        print("Number of variants as counted from source files: {0}".format(str(totVarsInSampleFiles)))
        print("Number of variant inserts recorded as retrieved from Cassandra: {0}".format(str(totVarsInCassandra)))
        print("Number of distinct variant positions as retrieved from Cassandra: {0}".format(
            str(totDistinctVarPosInCassandra)))
        print("******************************************************************************************************")

        return totVarsInCassandra, totVarsInSampleFiles == totVarsInCassandra


def initSparkEnv(cassandraNodeIPs):
    conf = SparkConf().setMaster("spark://{0}:7077".format(SSMergeCommonUtils.get_ip_address())) \
        .setAppName("SingleSampleVCFMerge").set("spark.cassandra.connection.host", cassandraNodeIPs[0]) \
        .set("spark.cassandra.read.timeout_ms", SSMergeCommonUtils.SPARK_CASSANDRA_READ_TIMEOUT_IN_MS) \
        .set("spark.cassandra.connection.timeout_ms", SSMergeCommonUtils.SPARK_CASSANDRA_READ_TIMEOUT_IN_MS)
    sc = SparkContext(conf=conf)
    sc.setLogLevel("WARN")
    return conf, sc


def initCassandraEnv(cassandraNodeIPs, studyName):
    keyspaceName = "variant_ksp"
    variantTableName = "variants_{0}".format(studyName.lower())
    headerTableName = "headers_{0}".format(studyName.lower())
    sampleInsertLogTableName = "sample_insert_log_{0}".format(studyName.lower())
    studyInfoTableName = "study_info_{0}".format(studyName.lower())
    uniquePosTableName = "uniq_pos_{0}".format(studyName.lower())
    sampleDefaultsTableName = "sample_defaults_{0}".format(studyName.lower())
    local_cluster_var = Cluster(cassandraNodeIPs)
    local_session_var = local_cluster_var.connect()
    local_session_var.execute("create keyspace if not exists {0}".format(keyspaceName) +
                              " with replication = {'class': 'SimpleStrategy', 'replication_factor' : 1};"
                              , timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(samplename varchar, header varchar, "
                              "primary key(samplename));"
                              .format(keyspaceName, headerTableName)
                              , timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(chrom varchar, chunk int, start_pos bigint, ref varchar, alt varchar, qual varchar, "
                              "filter varchar, info varchar, sampleinfoformat varchar, "
                              "sampleinfo  varchar, var_id varchar, var_uniq_id varchar, sampleName varchar,  "
                              "primary key((chrom, chunk), start_pos, var_uniq_id, samplename));"
                              .format(keyspaceName, variantTableName)
                              , timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(samplename varchar, proc_status varchar, num_variants bigint, primary key(samplename));"
                              .format(keyspaceName, sampleInsertLogTableName)
                              , timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(studyname varchar, tot_num_variants bigint, distinct_num_variant_pos bigint, "
                              "primary key(studyname));"
                              .format(keyspaceName, studyInfoTableName)
                              , timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(chrom varchar, chunk int, in_src int, start_pos bigint, "
                              "primary key((chrom,chunk,in_src), start_pos));"
                              .format(keyspaceName, uniquePosTableName)
                              , timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    local_session_var.execute("create table if not exists {0}.{1} (samplename varchar, default_genotype varchar, "
                              "primary key(samplename));"
                              .format(keyspaceName, sampleDefaultsTableName)
                              , timeout=SSMergeCommonUtils.SMALL_QUERY_TIMEOUT_IN_SECS)
    return keyspaceName, local_cluster_var, local_session_var, studyInfoTableName, variantTableName, \
           headerTableName,sampleInsertLogTableName, uniquePosTableName, sampleDefaultsTableName

# Helper function to write the unique variant positions from the first pass into a single file
def genVarPosFile(sparkContext, studyFilesInputDir, keyspaceName, uniquePosTableName):
    variantPositionFileName = studyFilesInputDir + os.path.sep + "unique_variant_positions.gz"
    sql = SQLContext(sparkContext)
    variants = sql.read.format("org.apache.spark.sql.cassandra"). \
        load(keyspace=keyspaceName, table=uniquePosTableName)
    variants.registerTempTable("variantsTable")
    resultDF = sql.sql("select chrom,start_pos from variantsTable order by 1,2")
    iterator = resultDF.toLocalIterator()
    with gzip.open(variantPositionFileName, "wb") as variantPositionFileHandle:
        for result in iterator:
            # For some reason, if you don't store the output of the write, every write result is echoed to the screen!!!
            discardOutput = variantPositionFileHandle.write(result["chrom"] + "\t" + str(result["start_pos"])
                                                            + os.linesep)
    return variantPositionFileName

def mainFunc():
    if len(sys.argv) != 8:
        print("Usage: SingleSampleMerge.py <Study PRJ ID> <Default Genotype> <Missing Genotype> "
              "<Full Path to study files> <Cassandra node IP1> <Cassandra node IP2> <BCF Tools Directory>")
        sys.exit(1)

    startTime = datetime.datetime.now()
    # Parse arguments
    studyName, defaultGenotype, missingGenotype, studyFilesInputDir = sys.argv[1:5]
    cassandraNodeIPs, bcfToolsDir = [sys.argv[5], sys.argv[6]], sys.argv[7]

    # Get the list of study files
    studyVCFFileNames = SSMergeCommonUtils.getVCFFileNames(studyFilesInputDir)

    # region Create keyspace and the required tables for variants, headers, sample insertion log, study info and
    # unique chromosome position tables in Cassandra
    keyspaceName, local_cluster_var, local_session_var, studyInfoTableName, variantTableName, headerTableName, \
    sampleInsertLogTableName, uniquePosTableName, sampleDefaultsTableName = \
        initCassandraEnv(cassandraNodeIPs, studyName)
    # endregion

    # region Initialize Apache Spark
    sparkConf, sparkContext = initSparkEnv(cassandraNodeIPs)
    # endregion

    # region Process study files in parallel
    numPartitions = len(studyVCFFileNames)
    studyIndivRDD = sparkContext.parallelize(studyVCFFileNames, numPartitions)
    results = studyIndivRDD.map(lambda entry: processStudyFiles(bcfToolsDir, studyName,
                                                                entry,
                                                                cassandraNodeIPs, keyspaceName, headerTableName,
                                                                variantTableName, sampleInsertLogTableName,
                                                                uniquePosTableName)).collect()
    # endregion

    # region Validate Cassandra record counts against record counts in the source file
    totVarsInCassandra, recordCountMatched = validateRecordCount(cassandraNodeIPs, keyspaceName, local_session_var,
                                                                 results, sparkContext, studyName, variantTableName,
                                                                 studyInfoTableName, uniquePosTableName)
    if not recordCountMatched:
        print("ERROR: Record counts did not match between the source study files and the Cassandra variant table")
    else:
        # Generate a position file with the unique set of variant positions determined from the first pass
        variantPositionFileName = genVarPosFile(sparkContext, studyFilesInputDir, keyspaceName, uniquePosTableName)

        # Scan each sample file against the unique set of variant positions
        # to insert sample genotypes at these positions
        numPartitions = len(studyVCFFileNames)
        studyIndivRDD = sparkContext.parallelize(studyVCFFileNames, numPartitions)
        processResults = studyIndivRDD.map(lambda studyVCFFileName:
                                           matchVariantPosInStudyFiles(bcfToolsDir, studyName,
                                                                       studyVCFFileName, studyFilesInputDir,
                                                                       variantPositionFileName,
                                                                       totVarsInCassandra, defaultGenotype,
                                                                       missingGenotype, cassandraNodeIPs,
                                                                       keyspaceName,
                                                                       sampleDefaultsTableName,
                                                                       variantTableName,
                                                                       sampleInsertLogTableName)).collect()
        for result in processResults:
            if result:
                print(result)
    # endregion

    # region Print elapsed time
    endTime = datetime.datetime.now()
    print("Time taken to process the study {0}: {1}".format(studyName, endTime - startTime))
    # endregion

    # region Shutdown Cassandra and Apache Spark
    local_session_var.shutdown()
    local_cluster_var.shutdown()
    sparkContext.stop()
    # endregion





if __name__ == '__main__':
    mainFunc()
