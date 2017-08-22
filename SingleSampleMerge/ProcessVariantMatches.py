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

from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext
import os, gzip, glob, traceback, subprocess, sys, SSMergeCommonUtils
from cassandra.cluster import Cluster

class ProcessVariantMatch:
    def __init__(self,*args):

        self.bcfToolsDir, self.studyName, self.studyVCFFileName, self.studyFilesInputDir, self.variantPositionFileName,\
        self.numTotVariants, self.defaultGenotype, self.missingGenotype, self.cassandraNodeIPs, self.keyspaceName,\
        self.sampleDefaultsTableName, self.variantTableName, self.sampleInsertLogTableName = args

        self.sampleName = SSMergeCommonUtils.getSampleName(self.bcfToolsDir, self.studyVCFFileName)
        self.cassandraCluster = Cluster(self.cassandraNodeIPs)
        self.cassandraSession = self.cassandraCluster.connect(self.keyspaceName)
        self.defaultGenotypeInsertPrepStmt = self.cassandraSession.prepare(
            "insert into {0}.{1} (samplename, default_genotype) values (?,?)".format(self.keyspaceName,
                                                                                     self.sampleDefaultsTableName))
        self.variantInsertPrepStmt = self.cassandraSession.prepare(
            "insert into {0}.{1} (chrom,chunk,start_pos,ref,alt,samplename, sampleinfoformat, sampleinfo) "
            "values (?,?,?,?,?,?,?,?)".format(self.keyspaceName, self.variantTableName))

    def insertDefaultGenotypeToCassandra(self):
        """
        Insert default genotype, assumed for a given sample, into Cassandra.            
        """
        self.cassandraSession.execute(self.defaultGenotypeInsertPrepStmt.bind([self.sampleName, self.defaultGenotype]))

    @staticmethod
    def getVariantInfo(matchLine):
        """
        Helper function to get variant information from the match file generated by this program 
        
        :param matchLine: A line from the match file generated by this program
        :type matchLine: str
        """
        chromMatch, varPosMatch, ref, alt, genotype, formatMatch, sampleinfoMatch = matchLine.strip().split("\t")
        varPosMatch = long(float(varPosMatch))
        chunk = int(varPosMatch / SSMergeCommonUtils.CHR_POS_CHUNKSIZE)
        return chromMatch, varPosMatch, ref, alt, genotype, formatMatch, sampleinfoMatch, chunk

    def processVariantMatchFile(self, matchOutputFileName):
        """
        Process a variant match file that is generated for every sample. The variant match file contains only entries 
        at the variant positions (determined from the first pass. see SingleSampleMerge.py) for a given sample.
        
        :param matchOutputFileName: Output file that contains the result of the sample file matched against the variant
                                    position file
        :type matchOutputFileName: str
        """
        # Insert default genotype assumed for the sample into Cassandra
        self.insertDefaultGenotypeToCassandra()

        # Compare the universal variant list against the matches found at variant positions in the individual sample files
        with gzip.open(self.variantPositionFileName, "rb") as varListFile:
            with gzip.open(matchOutputFileName, "rb") as varMatchFile:
                for matchLine in varMatchFile:
                    chromMatch, varPosMatch, ref, alt, genotype, \
                    formatMatch, sampleinfoMatch, chunk = self.getVariantInfo(matchLine)
                    for varPosLine in varListFile:
                        varPosLine = varPosLine.strip()
                        chromToFind, varPosToFind  = varPosLine.split("\t")
                        varPosToFind = long(float(varPosToFind))

                        # region Handle an edge case THAT SHOULD NOT HAPPEN
                        # If a variant position from the universal variant list is farther than
                        # the position in the match file. This can happen if the match file has multiple entries for
                        # the same position. In such cases, forward to a position in the match file
                        while varPosToFind > varPosMatch:
                            matchLine = varMatchFile.readline()
                            if not matchLine: return
                            chromMatch, varPosMatch, ref, alt, genotype, \
                            formatMatch, sampleinfoMatch, chunk = self.getVariantInfo(matchLine)
                        # endregion

                        # If we assumed the default was a "missing" genotype (ex: ./.) for a sample, insert to Cassandra
                        # if we have evidence otherwise i.e.,a match was found at a variant position in the sample.
                        #
                        # If we assumed the default was a "non-missing" genotype (ex: 0/0) for a sample, insert to Cassandra
                        # if we have evidence otherwise i.e.,no variant was found at a position in the sample (or)
                        # the variant found at the position had a different genotype than the default genotype (ex: 0/1)
                        if (chromToFind, varPosToFind) == (chromMatch, varPosMatch):
                            if genotype == self.defaultGenotype or alt != '.': break
                            self.cassandraSession.execute(self.variantInsertPrepStmt.bind(
                                    [chromToFind, chunk, varPosToFind, ref, alt,
                                     self.sampleName, formatMatch, sampleinfoMatch]))
                            break
                        else:
                            if self.defaultGenotype != "./." and self.defaultGenotype != ".|.":
                                self.cassandraSession.execute(self.variantInsertPrepStmt.bind(
                                    [chromToFind, chunk, varPosToFind, ref, alt,
                                     self.sampleName, "GT", self.missingGenotype]))


    def getNumMatchedVariants(self):
        """
        Get number of matched variant positions in a given sample when matched 
        against the universal set of variant positions from the first pass (see SingleSampleMerge.py)
        
        :return: Number of variant positions matched for the given sample
        :rtype: int
        """
        rows = self.cassandraSession.execute ("select num_variants from {0}.{1} where studyname = '{2}' "
                                              "and samplename = '{3}'"
                                              .format(self.keyspaceName, self.sampleInsertLogTableName,
                                                      self.studyName, self.sampleName))
        if rows:
            rows = iter(rows)
            return rows.next().num_variants
        else:
            raise Exception("Could not find number of matched variants for the sample: {0}".format(self.sampleName))

    def matchVariantPosInStudyFiles(self):
        returnErrMsg = None
        chosenFileToProcess = self.studyVCFFileName
        try:
            baseDir = self.studyFilesInputDir
            os.chdir(baseDir)
            sampleNameCmd = subprocess.Popen("{0}/bin/bcftools query -l {1}"
                                             .format(self.bcfToolsDir,chosenFileToProcess),shell=True,
                                                        stdout=subprocess.PIPE)
            sampleName, err = sampleNameCmd.communicate()
            if err: return "Could not obtain sample name from the SNP file:{0}".format(chosenFileToProcess)
            sampleName = sampleName.strip()

            matchOutputFileName = sampleName + "_variantmatch.gz"
            errFileName = sampleName + "_variantmatch.err"
            bcfVariantMatchCmd = "({0}/bin/bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\t%LINE' -T {1} {3} " \
                                 "| cut -f1,2,3,4,5,14,15 " \
                                 "| gzip) " \
                                 "1> {2} 2> {4}".format(
                self.bcfToolsDir, self.variantPositionFileName, matchOutputFileName, chosenFileToProcess, errFileName)

            os.system(bcfVariantMatchCmd)
            errFileHandle = open(errFileName, "r")
            errlines = errFileHandle.readlines()
            errFileHandle.close()
            if not errlines:
                numMatchedVariants = self.getNumMatchedVariants()
                if (self.numTotVariants * 1.0 / numMatchedVariants) > 2:
                    self.processVariantMatchFile(matchOutputFileName)
                else:
                    self.processVariantMatchFile(variantPositionFileName, matchOutputFileName, sampleName, defaultGenotype,
                                            missingGenotype)
            else:
                returnErrMsg = "Error in processing file:{0}".format(
                    chosenFileToProcess) + os.linesep + os.linesep.join(errlines)
        except Exception, e:
            returnErrMsg = "Error in processing file:{0}".format(
                chosenFileToProcess) + os.linesep + traceback.format_exc()
        finally:
            cassandraCluster.shutdown()
            cassandraSession.shutdown()
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
    procVarMatch = ProcessVariantMatch(bcfToolsDir, studyName, studyVCFFileName, studyFilesInputDir, variantPositionFileName,
                                numTotVariants, defaultGenotype, missingGenotype, cassandraNodeIPs, keyspaceName,
                                sampleDefaultsTableName, variantTableName, sampleInsertLogTableName)
    procVarMatch.matchVariantPosInStudyFiles()

def mainFunc():
    if len(sys.argv) != 8:
        print("Usage: ProcessVariantMatches.py <Study PRJ ID> <Default Genotype> <Missing Genotype> <Full Path to study files> <Cassandra node IP1> <Cassandra node IP2> <BCF Tools Directory>")
        sys.exit(1)

    # region Parse arguments
    studyName = sys.argv[1]
    defaultGenotype = sys.argv[2]
    missingGenotype = sys.argv[3]
    studyFilesInputDir = sys.argv[4]
    cassandraNodeIPs = [sys.argv[5], sys.argv[6]]
    bcfToolsDir = sys.argv[7]
    # endregion

    # region Get the list of study files
    os.chdir(studyFilesInputDir)
    dirContents = glob.glob("*_filtervcf.gz")
    dirContents.sort()
    studyVCFFileNames = dirContents
    # endregion

    # region Initialize variant, header, sample insert log and study table names in Cassandra
    keyspaceName = "variant_ksp"
    variantTableName = "variants_{0}".format(studyName.lower())
    sampleInsertLogTableName = "sample_insert_log"
    sampleDefaultsTableName = "sample_defaults_{0}".format(studyName.lower())
    studyInfoTableName = "study_info_{0}".format(studyName.lower())
    uniquePosTableName = "uniq_pos_{0}".format(studyName.lower())
    local_cluster_var = Cluster(cassandraNodeIPs)
    local_session_var = local_cluster_var.connect()
    local_session_var.execute("create table if not exists {0}.{1} (samplename varchar, default_genotype varchar, primary key(samplename));"
                              .format(keyspaceName,sampleDefaultsTableName))
    # endregion

    # region Initialize Apache Spark
    conf = SparkConf().setMaster("spark://{0}:7077".format(SSMergeCommonUtils.get_ip_address())).setAppName("SingleSampleVCFMerge").set("spark.cassandra.connection.host", cassandraNodeIPs[0]).set("spark.scheduler.listenerbus.eventqueue.size", "100000").set("spark.cassandra.read.timeout_ms", 1200000).set("spark.cassandra.connection.timeout_ms", 1200000)
    sc = SparkContext(conf=conf)
    sc.setLogLevel("INFO")
    # endregion

    # region Obtain the total number of variants from across the samples as determined from the first pass
    variantPositionFileName = studyFilesInputDir + os.path.sep + "unique_variant_positions.gz"
    rows = local_session_var.execute("select tot_num_variants from {0}.{1};".format(keyspaceName, studyInfoTableName))
    numTotVariants = 0
    if rows:
        rows = iter(rows)
        numTotVariants = rows.next().tot_num_variants
    else:
        raise Exception("Could not obtain number of variants for the study: {0} from the table: {1}.{2}"
                        .format(studyName, keyspaceName, studyInfoTableName))
    # endregion

    # region Generate a position file with the unique set of variant positions determined from the first pass
    sql = SQLContext(sc)
    variants = sql.read.format("org.apache.spark.sql.cassandra").\
                   load(keyspace=keyspaceName, table=uniquePosTableName)
    variants.registerTempTable("variantsTable")
    resultDF = sql.sql("select chrom,start_pos from variantsTable order by 1,2")
    iterator = resultDF.toLocalIterator()

    variantPositionFileHandle = gzip.open(variantPositionFileName, "wb")
    for result in iterator:
        discardOutput = variantPositionFileHandle.write(result["chrom"] + "\t" + str(result["start_pos"]) + os.linesep)
    variantPositionFileHandle.close()
    # endregion

    # region Scan each sample file against the unique set of variant positions
    # to insert sample genotypes at these positions
    numPartitions = len(studyVCFFileNames)
    studyIndivRDD = sc.parallelize(studyVCFFileNames, numPartitions)
    processResults = studyIndivRDD.map(lambda studyVCFFileName: matchVariantPosInStudyFiles(bcfToolsDir, studyName,
                                                                                         studyVCFFileName, studyFilesInputDir,
                                                                                         variantPositionFileName,
                                                                                         numTotVariants, defaultGenotype,
                                                                                         missingGenotype, cassandraNodeIPs,
                                                                                         keyspaceName,
                                                                                         sampleDefaultsTableName,
                                                                                         variantTableName,
                                                                                         sampleInsertLogTableName)).collect()
    # endregion

    # region Print any error messages obtained from the process above
    for result in processResults:
        if result:
            print(result)
    # endregion

    # region Shutdown Cassandra and Apache Spark
    local_session_var.shutdown()
    local_cluster_var.shutdown()
    sc.stop()
    # endregion

if __name__ == "__main__":
    mainFunc()