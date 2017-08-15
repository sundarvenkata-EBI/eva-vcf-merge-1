# import ftplib
import os, hashlib, sys, glob, socket, gzip
import traceback, subprocess

from cassandra.cluster import Cluster
from cassandra.query import BatchStatement, BatchType
from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext

def getSampleName(bcfToolsDir, vcfFileName):
    sampleNameCmd = subprocess.Popen(
        "{0}/bin/bcftools query -l {1}".format(bcfToolsDir, vcfFileName), shell=True,
        stdout=subprocess.PIPE)
    sampleName, err = sampleNameCmd.communicate()
    if err: return err
    return sampleName.strip()

def get_ip_address():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    return s.getsockname()[0]

def writeVariantToCassandra(linesToWrite, sampleName):
        chunkSize = 1000000
        batch = BatchStatement(BatchType.UNLOGGED)
        for line in linesToWrite:
            lineComps = line.split("\t")
            chromosome = lineComps[0]
            position = long(float(lineComps[1]))
            rsID = lineComps[2].strip()
            ref = lineComps[3].strip()
            alt = lineComps[4].strip()
            qual = lineComps[5].strip()
            qualFilter = lineComps[6].strip()
            info = lineComps[7].strip()
            sampleInfoFormat = lineComps[8].strip()
            sampleInfo = lineComps[9].strip()
            chunk = int(position/chunkSize)
            variantID = chromosome + "_" + str(position).zfill(12) + "_" + hashlib.md5(ref + "_" + alt).hexdigest() + sampleName.zfill(12)
            boundStmt = variantInsertStmt.bind([chromosome, chunk, position, ref, alt, qual, qualFilter, info, sampleInfoFormat, sampleInfo, rsID, variantID, sampleName])
            batch.add(boundStmt)
        session.execute(batch, timeout = 1200)


def writeHeaderToCassandra(headerLines, sampleName):
    session.execute(headerInsertStmt.bind([sampleName, headerLines]), timeout=1200)


def cassandraInsert(vcfFileName):
    totVarsInSample = 0
    vcfFileHandle = gzip.open(vcfFileName, 'rb')
    headerLines = ""
    sampleName = ""
    lineBatchSize = 50
    linesToWrite = []
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
    vcfFileHandle.close()
    return totVarsInSample


def getErrFileContents(errFileName):
    errFileContents = None
    with open(errFileName, "r") as errFileHandle:
        errFileContents = errFileHandle.readlines()
    if errFileContents: return errFileContents
    return None

def getSampleProcessedStatus(keyspaceName, sampleInsertLogTableName, studyName, sampleName):
    allrows = session.execute("select insert_flag from {0}.{1} where studyname = '{2}' and samplename = '{3}';".format(keyspaceName, sampleInsertLogTableName, studyName, sampleName))
    if allrows:
        allrows = iter(allrows)
        firstRow = allrows.next()
    else:
        return False
    if firstRow.insert_flag == 1: return True
    return True

def processStudyFiles(studyName, studyFileName, cassandraNodeIPs, keyspaceName, headerTableName, variantTableName, sampleInsertLogTableName, bcfToolsDir):

    totVarsInSample = 0
    filterCommandResult = -1
    global cluster, session, variantInsertStmt, headerInsertStmt
    samplePrefix, errFileContents, returnErrMsg, cluster, session = None, None, None, None, None
    cluster = Cluster(cassandraNodeIPs)
    session = cluster.connect(keyspaceName)

    sampleName = getSampleName(bcfToolsDir, studyFileName)
    isSampleProcessed = getSampleProcessedStatus(keyspaceName, sampleInsertLogTableName, studyName, sampleName)
    if isSampleProcessed: return (-1, "Already processed sample: {0} for study: {1}".format(sampleName, studyName))

    try:
        variantInsertStmt = session.prepare(
            "insert into {0}.{1} (chrom,chunk,start_pos,ref,alt,qual,filter,info,sampleinfoformat,sampleinfo,var_id,var_uniq_id,sampleName) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)".format(keyspaceName, variantTableName))
        headerInsertStmt = session.prepare("INSERT INTO {0}.{1} (samplename, header) VALUES (?,?)".format(keyspaceName,headerTableName))

        samplePrefix = os.path.basename(studyFileName).split(".")[0]
        baseDir = os.path.dirname(studyFileName)
        filteredFileName = "{0}_filtervcf.gz".format(samplePrefix)
        os.chdir(baseDir)

        if not os.path.isfile(baseDir + os.path.sep + filteredFileName):
            os.system("""{0}/bin/bcftools filter -e ALT=\\'.\\' -o {2}_filtervcf.gz -O z {1} 2> {2}_filtering_err.txt""".format(bcfToolsDir,studyFileName, samplePrefix))
            errFileContents = getErrFileContents("{0}_filtering_err.txt".format(samplePrefix))
            if errFileContents: filterCommandResult = -1
            else: filterCommandResult = 0
        else:
            filterCommandResult = 0

        if filterCommandResult != 0:
             returnErrMsg = "Failed to process {0} due to error: {1}".format(studyFileName, errFileContents)
        else:
            totVarsInSample = cassandraInsert(baseDir + os.path.sep + "{0}_filtervcf.gz".format(samplePrefix))
            os.system("echo {0} > {1}_filtered_variant_count.txt".format(str(totVarsInSample), samplePrefix))

        if not returnErrMsg and totVarsInSample != 0:
            session.execute(
                "insert into {0}.{1} (studyname, samplename, insert_flag, num_variants) values ('{2}', '{3}', {4}, {5})".format(
                    keyspaceName, sampleInsertLogTableName, studyName, sampleName, 1, totVarsInSample))
    except Exception, ex:
        returnErrMsg = "Error in processing file:{0}".format(studyFileName) + os.linesep + traceback.format_exc()
    finally:
        if cluster != None and session != None:
            try:
                session.shutdown()
                cluster.shutdown()
            except Exception, e:
                pass
        if returnErrMsg: return (-1, returnErrMsg)
        return (totVarsInSample, None)



# studyFilesDir = "/pub/databases/eva/PRJEB13618/submitted_files"
# ftpSite = "ftp.ebi.ac.uk"
# ftpUserName = "anonymous"
# ftp = ftplib.FTP(ftpSite, ftpUserName)
# ftp.cwd(studyFilesDir)
if len(sys.argv) != 6:
    print("Usage: SingleSampleMerge.py <Study Name> <Full Path to study files> <Cassandra node IP1> <Cassandra node IP2> <BCF Tools Directory>")
    sys.exit(1)

# Parse arguments
studyName = sys.argv[1]
studyFilesInputDir = sys.argv[2]
cassandraNodeIPs = [sys.argv[3], sys.argv[4]]
bcfToolsDir = sys.argv[5]

# Get the list of study files
os.chdir(studyFilesInputDir)
dirContents = glob.glob("*.vcf.gz")
dirContents.sort()
studyFileNames = dirContents

# Create keyspace, variant, header and sample insertlog tables in Cassandra
keyspaceName = "variant_ksp"
variantTableName = "variants_{0}".format(studyName.lower())
headerTableName = "headers_{0}".format(studyName.lower())
sampleInsertLogTableName = "sample_insert_log"
studyInfoTableName = "study_info_{0}".format(studyName.lower())
local_cluster_var = Cluster(cassandraNodeIPs)
local_session_var = local_cluster_var.connect()
local_session_var.execute("create keyspace if not exists {0}".format(keyspaceName)  + " with replication = {'class': 'SimpleStrategy', 'replication_factor' : 1};")
local_session_var.execute("create table if not exists {0}.{1} (samplename varchar, header varchar, primary key(samplename));".format(keyspaceName,headerTableName))
local_session_var.execute("create table if not exists {0}.{1} (chrom varchar, chunk int, start_pos bigint, ref varchar, alt varchar, qual varchar, filter varchar, info varchar, sampleinfoformat varchar, sampleinfo  varchar, var_id varchar, var_uniq_id varchar, sampleName varchar,  primary key((chrom, chunk), start_pos, ref, alt, samplename));".format(keyspaceName,variantTableName))
local_session_var.execute("create table if not exists {0}.{1} (studyname varchar, samplename varchar, insert_flag int, num_variants bigint, primary key((studyname, samplename)));".format(keyspaceName,sampleInsertLogTableName))
local_session_var.execute("create table if not exists {0}.{1} (studyname varchar, tot_num_variants bigint, distinct_num_variants bigint, primary key(studyname));".format(keyspaceName,studyInfoTableName))

# Process study files in parallel
conf = SparkConf().setMaster("spark://{0}:7077".format(get_ip_address())).setAppName("SingleSampleVCFMerge").set("spark.cassandra.connection.host", cassandraNodeIPs[0]).set("spark.cassandra.read.timeout_ms", 1200000).set("spark.cassandra.connection.timeout_ms", 1200000)
sc = SparkContext(conf=conf)
sc.setLogLevel("INFO")
numPartitions = len(studyFileNames)
studyIndivRDD = sc.parallelize(studyFileNames, numPartitions)
results = studyIndivRDD.map(lambda entry: processStudyFiles(studyName, studyFilesInputDir + os.path.sep + entry, cassandraNodeIPs, keyspaceName, headerTableName, variantTableName, sampleInsertLogTableName, bcfToolsDir)).collect()


# region Validate Cassandra record counts against record counts in the source file
totNumVariantsInSourceFiles = 0
totNumVariantsInserted = 0

errMsg = None
for totVarsInSample, errMsg in results:
    if errMsg:
        print(errMsg)
    else:
        totNumVariantsInSourceFiles += totVarsInSample

if not errMsg:
    procResult = subprocess.Popen('zcat {0}/*.vcf.gz | grep -v "^#" | cut -f1,2 | sort | uniq | wc -l'.format(studyFilesInputDir), shell=True,
            stdout=subprocess.PIPE)
    totNumDistinctVariantsInSourceFiles, errMsg = procResult.communicate()

    if errMsg: print("Error retrieving unique set of variants from source files: {0}".format(errMsg))
    else:
        # Count the number of variants from the variants table
        sql = SQLContext(sc)
        variants = sql.read.format("org.apache.spark.sql.cassandra"). \
            load(keyspace=keyspaceName,table=variantTableName)
        variants.registerTempTable("variantsTable")
        resultDF = sql.sql("select chrom,start_pos from variantsTable")
        totNumVariantsInserted = resultDF.count()
        resultDF = sql.sql("select chrom,start_pos from variantsTable group by 1,2")
        totNumDistinctVariantsInserted = resultDF.count()
        local_session_var.execute("insert into {0}.{1} (studyname, tot_num_variants, distinct_num_variants) values ('{2}',{3}, {4})".format(keyspaceName, studyInfoTableName, studyName, totNumVariantsInserted, totNumDistinctVariantsInserted))

        print("******************************************************************************************************")
        print("Number of variants as counted from source files: {0}".format(str(totNumVariantsInSourceFiles)))
        print("Number of distinct variants as counted from source files: {0}".format(str(totNumDistinctVariantsInSourceFiles)))

        print("Number of variant inserts recorded as retrieved from Cassandra: {0}".format(str(totNumVariantsInserted)))
        print("Number of distinct variant inserts as retrieved from Cassandra: {0}".format(str(totNumDistinctVariantsInserted)))
        print("******************************************************************************************************")
    # endregion

local_session_var.shutdown()
local_cluster_var.shutdown()

sc.stop()