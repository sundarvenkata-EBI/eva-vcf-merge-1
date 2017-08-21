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

import os, hashlib, sys, glob, gzip, SSMergeCommonUtils
import traceback, subprocess

from cassandra.cluster import Cluster
from cassandra.query import BatchStatement, BatchType
from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext

def stringDiffIndex(string1, string2):
    """
    Return index of difference between two strings
    For examples, see https://commons.apache.org/proper/commons-lang/apidocs/org/apache/commons/lang3/StringUtils.html#indexOfDifference-java.lang.CharSequence...-
    """
    if string1 is None or string2 is None: return -1
    string1Length = len(string1)
    string2Length = len(string2)
    shorterLength = min(string1Length, string2Length)
    longerLength = max(string1Length, string2Length)
    index = -1
    for i in xrange(longerLength):
        if (i == shorterLength) or (string1[i] != string2[i]):
            index = i
            break
    return index

def getNormalizedStartPos(start_pos, ref, alt):
    """
    Get a normalized start position for a given variant.
    This is a verbatim implementation of https://github.com/EBIvariation/eva-pipeline/blob/develop/src/main/java/uk/ac/ebi/eva/pipeline/io/mappers/VariantVcfFactory.java#L190
    
    :param start_pos: Variant start position
    :type start_pos: long
    :param ref: Variant reference allele
    :type ref: str
    :param alt: Variant alternate allele
    :type alt: str
    :return: The normalized start position
    :rtype: long
    """
    # Remove trailing bases
    refReversed = ref[::-1]
    altReversed = alt[::-1]
    indexOfDifference = stringDiffIndex(refReversed, altReversed)
    ref = refReversed[indexOfDifference:][::-1]
    alt = altReversed[indexOfDifference:][::-1]
    # Remove leading bases
    indexOfDifference = stringDiffIndex(ref,alt)
    start = start_pos + indexOfDifference
    return start

def getSampleName(bcfToolsDir, vcfFileName):
    """
    Get sample name from a VCF file
    
    :param bcfToolsDir: Directory containing the bcftools static binaries
    :type bcfToolsDir: str
    :param vcfFileName: Full path to the VCF file
    :type vcfFileName: str
    :return: Sample name
    :rtype: str
    """
    sampleNameCmd = subprocess.Popen("{0}/bin/bcftools query -l {1}".format(bcfToolsDir, vcfFileName), shell=True, stdout=subprocess.PIPE)
    sampleName, err = sampleNameCmd.communicate()
    if err: return err
    return sampleName.strip()

def writeVariantToCassandra(linesToWrite, sampleName):
    """
    Write a set of variant lines to Cassandra
    
    :param linesToWrite: Set of variant lines to write from a given sample
    :type linesToWrite: list[str]
    :param sampleName: Name of the sample that houses the variant lines
    :type sampleName: str
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
        chunk = int(position/SSMergeCommonUtils.CHR_POS_CHUNKSIZE)

        # Insert the chromosome position that was scanned
        batch.add(uniquePosInsertStmt.bind([chromosome, chunk, position]))
        # Also, insert the normalized chromosome positions for each alternate allele
        for altAllele in alt.split(","):
            revised_start_pos = getNormalizedStartPos(position, ref, altAllele)
            if revised_start_pos != position: batch.add(uniquePosInsertStmt.bind([chromosome, chunk, revised_start_pos]))

        # Unique Variant ID is chromosome_12-digit zero padded start_MD5(REF_ALT)_20-character samplename
        variantID = chromosome + "_" + str(position).zfill(12) + "_" + hashlib.md5(ref + "_" + alt).hexdigest() + "_" + sampleName.zfill(20)
        # Insert the variant information
        boundStmt = variantInsertStmt.bind([chromosome, chunk, position, ref, alt, qual, qualFilter, info, sampleInfoFormat, sampleInfo, rsID, variantID, sampleName])
        batch.add(boundStmt)

    # Batch execute the insert statements added in the loop above
    cassandraSession.execute(batch, timeout = 1200)

def writeHeaderToCassandra(headerLines, sampleName):
    """
    Write a set of header lines to Cassandra

    :param headerLines: Set of header lines to write from a given sample
    :type headerLines: str
    :param sampleName: Name of the sample that houses the header lines
    :type sampleName: str
    :return: None
    """
    cassandraSession.execute(headerInsertStmt.bind([sampleName, headerLines]), timeout=1200)

def cassandraInsert(vcfFileName):
    """
    Insert the contents of a VCF file: variants and headers into relevant Cassandra tables 
    
    :param vcfFileName: Full path to the VCF file
    :type vcfFileName: str
    :return: Total number of variants that were scanned from the VCF file
    :rtype: long
    """
    lineBatchSize = 50 # Write variants in batches
    totVarsInSample = 0
    headerLines, sampleName = "", ""
    linesToWrite = []

    with gzip.open(vcfFileName, 'rb') as vcfFileHandle:
        for line in vcfFileHandle:
            line = line.strip()
            if line.startswith("#"):
                headerLines += (line + os.linesep)
                if (line.startswith("#CHROM")):
                    sampleName = line.split("\t")[-1]
                    writeHeaderToCassandra(headerLines.strip(), sampleName)
                    break
        lineBatchIndex = 0
        for line in vcfFileHandle:
            totVarsInSample += 1
            line = line.strip()
            linesToWrite.append(line)
            lineBatchIndex += 1
            if lineBatchIndex == lineBatchSize:
                writeVariantToCassandra(linesToWrite, sampleName)
                lineBatchIndex = 0
                linesToWrite = []
        if linesToWrite:
            writeVariantToCassandra(linesToWrite, sampleName)

        return totVarsInSample

def getSampleProcessedStatus(keyspaceName, sampleInsertLogTableName, studyName, sampleName):
    """
    Check if the variants from a specific sample have already been inserted into the requisite Cassandra tables
    
    :param keyspaceName: Cassandra keyspace
    :type keyspaceName: str
    :param sampleInsertLogTableName: Cassandra Sample Insert log table name
    :type sampleInsertLogTableName: str
    :param studyName: Name of the study
    :type studyName: str
    :param sampleName: Name of the sample
    :type sampleName: str
    :return: True if variants were inserted for the given sample, False otherwise
    :rtype: bool
    """
    allrows = cassandraSession.execute("select proc_status from {0}.{1} where studyname = '{2}' and samplename = '{3}';"
                              .format(keyspaceName, sampleInsertLogTableName, studyName, sampleName))
    if allrows:
        allrows = iter(allrows)
        firstRow = allrows.next()
    else:
        return False
    if firstRow.proc_status: return firstRow.proc_status
    return ""


def getTotDistinctVarPosFromSampleFiles(vcfInputDirectory):
    """
    Get the total number of distinct variant positions from all the VCF files across all the samples
    
    :param vcfInputDirectory: Full path to the VCF directory
    :type vcfInputDirectory: str
    :return: Total number of distinct variants
    :rtype: long
    """
    procResult = subprocess.Popen(
        'zcat {0}/*_filtervcf.gz | grep -v "^#" | cut -f1,2 | sort | uniq | wc -l'.format(vcfInputDirectory),
        shell=True,
        stdout=subprocess.PIPE)
    totNumDistinctVarPos, errMsg = procResult.communicate()
    if errMsg: raise Exception("Could not get distinct number of variants due to error:{0}".format(errMsg))
    else: return long(float(totNumDistinctVarPos))

def update_processing_status(keyspaceName, sampleInsertLogTableName, sampleName,
                             studyName, totVarsInSample, proc_status):
    """
    Helper function - Update processing status as the sample files are being processed. 
    
    :param keyspaceName: Cassandra key space with the variant table
    :type keyspaceName: str
    :param sampleInsertLogTableName: Cassandra sample insert log table to record the number of variants in the sample
    :type sampleInsertLogTableName: str
    :param sampleName: Sample name
    :type sampleName: str
    :param studyName: Study name (ex: PRJEB21300)
    :type studyName: str
    :param totVarsInSample: Total number of variants in the sample
    :type totVarsInSample: long
    :param proc_status: Processing status ("variants_filtered" or "variants_inserted")
    :type proc_status: str
    :return: 
    """
    cassandraSession.execute(
        "insert into {0}.{1} (studyname, samplename, proc_status, num_variants) values ('{2}', '{3}', '{4}', {5})".format(
            keyspaceName, sampleInsertLogTableName, studyName, sampleName, proc_status, totVarsInSample))


def processStudyFiles(bcfToolsDir, studyName, studyFileName, cassandraNodeIPs, keyspaceName, headerTableName,
                      variantTableName, sampleInsertLogTableName, uniquePosTableName):
    """
    Process a study file (runs on each core on each Apache Spark node)
    
    :param studyName: Name of the study (ex: PRJEB21300)
    :type studyName: str
    :param studyFileName: Full path to the study file (ex: /mnt/shared_storage/PRJEB21300/StudySample1.vcf.gz)
    :type studyFileName: str
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
    :return: Tuple - (Total number of variants scanned, Error message if any) 
    """
    global cassandraCluster, cassandraSession, variantInsertStmt, headerInsertStmt, uniquePosInsertStmt
    samplePrefix, errFileContents, returnErrMsg, cassandraCluster, cassandraSession = None, None, None, None, None
    totVarsInSample = 0

    cassandraCluster = Cluster(cassandraNodeIPs)
    cassandraSession = cassandraCluster.connect(keyspaceName)

    sampleName = getSampleName(bcfToolsDir, studyFileName)

    sample_proc_status = getSampleProcessedStatus(keyspaceName, sampleInsertLogTableName, studyName, sampleName)
    if sample_proc_status == "variants_inserted": return (-1,"Already processed sample: {0} for study: {1}"
                                                          .format(sampleName, studyName))

    try:
        # Prepared statements to insert into variants, headers and unique positions tables in Cassandra
        variantInsertStmt = cassandraSession.prepare(
            "insert into {0}.{1} "
            "(chrom,chunk,start_pos,ref,alt,qual,filter,info,sampleinfoformat,sampleinfo,var_id,var_uniq_id,sampleName) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)".format(keyspaceName, variantTableName))
        headerInsertStmt = cassandraSession.prepare("insert into {0}.{1} (samplename, header) VALUES (?,?)"
                                                    .format(keyspaceName,headerTableName))
        uniquePosInsertStmt = cassandraSession.prepare("insert into {0}.{1} (chrom, chunk, start_pos) VALUES (?,?,?)"
                                                       .format(keyspaceName,uniquePosTableName))

        # Get the sample prefix and the base directory of the file that we are operating on
        samplePrefix = os.path.basename(studyFileName).split(".")[0]
        baseDir = os.path.dirname(studyFileName)
        os.chdir(baseDir)

        # region Run bcftools filter command to filter out monomorphic references
        if sample_proc_status != "variants_filtered":
            os.system("""({0}/bin/bcftools filter -e ALT=\\'.\\' {1} | gzip) 1>{2}_filtervcf.gz 2>{2}_filtering_err.txt"""
                      .format(bcfToolsDir,studyFileName, samplePrefix))
            errFileContents = SSMergeCommonUtils.getErrFileContents("{0}_filtering_err.txt".format(samplePrefix))
            if errFileContents:
                filterCommandResult = -1
            else:
                filterCommandResult = 0
                update_processing_status(keyspaceName, sampleInsertLogTableName, sampleName, studyName,
                                         0, "variants_filtered")
        else:
            filterCommandResult = 0
        # endregion

        # region Insert the variants into Cassandra from the filtered VCF files obtained from above
        if filterCommandResult != 0:
             returnErrMsg = "Failed to process {0} due to error: {1}".format(studyFileName, errFileContents)
        else:
            totVarsInSample = cassandraInsert(baseDir + os.path.sep + "{0}_filtervcf.gz".format(samplePrefix))
            os.system("echo {0} > {1}_filtered_variant_count.txt".format(str(totVarsInSample), samplePrefix))
        # endregion

        # region Update status to indicate that all variants have been inserted into Cassandra
        if not returnErrMsg and totVarsInSample != 0:
            update_processing_status(keyspaceName, sampleInsertLogTableName, sampleName, studyName,
                                     totVarsInSample, "variants_inserted")
        # endregion

    except Exception, ex:
        returnErrMsg = "Error in processing file:{0}".format(studyFileName) + os.linesep + ex.message \
                       + os.linesep + traceback.format_exc()
    finally:
        if cassandraCluster != None and cassandraSession != None:
            try:
                cassandraSession.shutdown()
                cassandraCluster.shutdown()
            except Exception, e:
                pass
        if returnErrMsg: return (-1, returnErrMsg)
        return (totVarsInSample, None)



if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("Usage: SingleSampleMerge.py <Study PRJ ID> <Full Path to study files> "
              "<Cassandra node IP1> <Cassandra node IP2> <BCF Tools Directory>")
        sys.exit(1)

    # region Parse arguments
    studyName = sys.argv[1]
    studyFilesInputDir = sys.argv[2]
    cassandraNodeIPs = [sys.argv[3], sys.argv[4]]
    bcfToolsDir = sys.argv[5]
    # endregion

    # region Get the list of study files
    os.chdir(studyFilesInputDir)
    dirContents = glob.glob("*.vcf.gz")
    dirContents.sort()
    studyFileNames = dirContents
    # endregion

    # region Create keyspace and the required tables for variants, headers, sample insertion log, study info and
    # unique chromosome position tables in Cassandra
    keyspaceName = "variant_ksp"
    variantTableName = "variants_{0}".format(studyName.lower())
    headerTableName = "headers_{0}".format(studyName.lower())
    sampleInsertLogTableName = "sample_insert_log"
    studyInfoTableName = "study_info_{0}".format(studyName.lower())
    uniquePosTableName = "uniq_pos_{0}".format(studyName.lower())
    local_cluster_var = Cluster(cassandraNodeIPs)
    local_session_var = local_cluster_var.connect()
    local_session_var.execute("create keyspace if not exists {0}".format(keyspaceName)  +
                              " with replication = {'class': 'SimpleStrategy', 'replication_factor' : 1};")
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(samplename varchar, header varchar, "
                              "primary key(samplename));"
                              .format(keyspaceName,headerTableName))
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(chrom varchar, chunk int, start_pos bigint, ref varchar, alt varchar, qual varchar, "
                              "filter varchar, info varchar, sampleinfoformat varchar, "
                              "sampleinfo  varchar, var_id varchar, var_uniq_id varchar, sampleName varchar,  "
                              "primary key((chrom, chunk), start_pos, ref, alt, samplename));"
                              .format(keyspaceName,variantTableName))
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(studyname varchar, samplename varchar, proc_status varchar, num_variants bigint, "
                              "primary key((studyname, samplename)));"
                              .format(keyspaceName,sampleInsertLogTableName))
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(studyname varchar, tot_num_variants bigint, distinct_num_variant_pos bigint, "
                              "primary key(studyname));"
                              .format(keyspaceName,studyInfoTableName))
    local_session_var.execute("create table if not exists {0}.{1} "
                              "(chrom varchar, chunk int, start_pos bigint, "
                              "primary key((chrom,chunk), start_pos));"
                              .format(keyspaceName,uniquePosTableName))
    # endregion

    # region Initialize Apache Spark
    conf = SparkConf().setMaster("spark://{0}:7077".format(SSMergeCommonUtils.get_ip_address())).setAppName(
        "SingleSampleVCFMerge").set("spark.cassandra.connection.host", cassandraNodeIPs[0]).set(
        "spark.cassandra.read.timeout_ms", 1200000).set("spark.cassandra.connection.timeout_ms", 1200000)
    sc = SparkContext(conf=conf)
    sc.setLogLevel("INFO")
    # endregion

    # region Process study files in parallel
    numPartitions = len(studyFileNames)
    studyIndivRDD = sc.parallelize(studyFileNames, numPartitions)
    results = studyIndivRDD.map(lambda entry: processStudyFiles(bcfToolsDir, studyName,
                                                                studyFilesInputDir + os.path.sep + entry,
                                                                cassandraNodeIPs, keyspaceName, headerTableName,
                                                                variantTableName, sampleInsertLogTableName,
                                                                uniquePosTableName)).collect()
    # endregion

    # region Validate Cassandra record counts against record counts in the source file
    atleastOneError = False
    totVarsInSampleFiles = 0

    for (totVarsInSample, errMsg) in results:
        if errMsg:
            print(errMsg)
            atleastOneError = True
        else:
            totVarsInSampleFiles += totVarsInSample

    totDistinctVarPosInSampleFiles = getTotDistinctVarPosFromSampleFiles(studyFilesInputDir)
    if not atleastOneError:
        # Count the number of variants from the variants table
        sql = SQLContext(sc)
        variants = sql.read.format("org.apache.spark.sql.cassandra"). \
            load(keyspace=keyspaceName,table=variantTableName)
        variants.registerTempTable("variantsTable")
        resultDF = sql.sql("select chrom,start_pos from variantsTable")
        totNumVariantsInserted = resultDF.count()
        resultDF = sql.sql("select chrom,start_pos from variantsTable group by 1,2")
        totNumDistinctVariantsInserted = resultDF.count()

        # Register counts of
        local_session_var.execute("insert into {0}.{1} (studyname, tot_num_variants, distinct_num_variant_pos) "
                                  "values ('{2}',{3}, {4})"
                                  .format(keyspaceName, studyInfoTableName, studyName, totNumVariantsInserted,
                                          totNumDistinctVariantsInserted))

        print("******************************************************************************************************")
        print("Number of variants as counted from source files: {0}".format(str(totVarsInSampleFiles)))
        print("Number of distinct variant positions as counted from source files: {0}".format(str(totDistinctVarPosInSampleFiles)))

        print("Number of variant inserts recorded as retrieved from Cassandra: {0}".format(str(totNumVariantsInserted)))
        print("Number of distinct variant positions as retrieved from Cassandra: {0}".format(str(totNumDistinctVariantsInserted)))
        print("******************************************************************************************************")
    # endregion

    # region Shutdown Cassandra and Apache Spark
    local_session_var.shutdown()
    local_cluster_var.shutdown()
    sc.stop()
    # endregion