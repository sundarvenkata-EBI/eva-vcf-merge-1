from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext
import os, gzip, glob, traceback, subprocess, sys, socket
from cassandra.cluster import Cluster

def gzipfileLineCount(fname):
    count = 0
    with gzip.open(fname,"rb") as f:
        for line in f:
            count += 1
    return count

def insertDefaultGenotypeToCassandra(sampleName, defaultGenotype):
    session.execute(defaultGenotypeInsertPrepStmt.bind([sampleName, defaultGenotype]))

def get_ip_address():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    return s.getsockname()[0]

def processVariantMatchFile(uniqueVariantListFileName, matchOutputFileName, sampleName, defaultGenotype):
    insertDefaultGenotypeToCassandra(sampleName, defaultGenotype)
    with gzip.open(uniqueVariantListFileName, "rb") as varListFile:
        with gzip.open(matchOutputFileName, "rb") as varMatchFile:
            if defaultGenotype != "./.":
                for matchLine in varMatchFile:
                    chromMatch, varPosMatch, ref, alt, genotype, formatMatch, sampleinfoMatch = matchLine.strip().split("\t")
                    varPosMatch = long(varPosMatch)
                    chunk = int(varPosMatch/1000000)
                    for varPosLine in varListFile:
                        varPosLine = varPosLine.strip()
                        chromToFind, varPosToFind  = varPosLine.split("\t")
                        varPosToFind = long(varPosToFind)
                        while varPosToFind > varPosMatch:
                            matchLine = varMatchFile.readline()
                            if not matchLine: return
                            chromMatch, varPosMatch, ref, alt, genotype, formatMatch, sampleinfoMatch = matchLine.strip().split("\t")
                            varPosMatch = long(varPosMatch)
                            chunk = int(varPosMatch / 1000000)
                        if (chromToFind, varPosToFind) == (chromMatch,varPosMatch):
                            if genotype == defaultGenotype or alt != '.': break
                            session.execute(variantInsertPrepStmt.bind(
                                [chromToFind, chunk, varPosToFind, ref, alt, sampleName, formatMatch, sampleinfoMatch]))
                            break
                        else:
                            session.execute(variantInsertPrepStmt.bind(
                                [chromToFind, chunk, varPosToFind, ref, alt, sampleName, "GT", "./."]))
            else:
                for matchLine in varMatchFile:
                    chromMatch, varPosMatch, ref, alt, genotype, formatMatch, sampleinfoMatch = matchLine.strip().split("\t")
                    varPosMatch = long(varPosMatch)
                    chunk = int(varPosMatch / 1000000)
                    for varPosLine in varListFile:
                        varPosLine = varPosLine.strip()
                        chromToFind, varPosToFind  = varPosLine.split("\t")
                        varPosToFind = long(float(varPosToFind))
                        while varPosToFind > varPosMatch:
                            matchLine = varMatchFile.readline()
                            if not matchLine: return
                            chromMatch, varPosMatch, ref, alt, genotype, formatMatch, sampleinfoMatch = matchLine.strip().split("\t")
                            varPosMatch = long(varPosMatch)
                            chunk = int(varPosMatch / 1000000)
                        if (chromToFind, varPosToFind) == (chromMatch,varPosMatch):
                            if genotype == defaultGenotype or alt != '.': break
                            if genotype != defaultGenotype:
                                session.execute(variantInsertPrepStmt.bind(
                                    [chromToFind, chunk, varPosToFind, ref, alt, sampleName, formatMatch,
                                     sampleinfoMatch]))
                            break


def getNumMatchedVariants(keyspaceName, sampleInsertLogTableName, studyName, sampleName):
    rows = session.execute ("select num_variants from {0}.{1} where studyname = '{2}' and samplename = '{3}'".format(keyspaceName, sampleInsertLogTableName, studyName, sampleName))
    if rows:
        rows = iter(rows)
        return rows.next().num_variants
    else:
        raise Exception("Could not find number of matched variants for the sample: {0}".format(sampleName))


def matchVariantPosInStudyFiles(bcfToolsDir, studyName, defaultRefGenotype, studyFileName, numTotVariants, studyFilesInputDir, variantPositionFileName, cassandraNodeIPs, keyspaceName, sampleDefaultsTableName, variantTableName, sampleInsertLogTableName):
    global cluster, session, defaultGenotypeInsertPrepStmt, variantInsertPrepStmt
    cluster = Cluster(cassandraNodeIPs)
    session = cluster.connect(keyspaceName)
    defaultGenotypeInsertPrepStmt = session.prepare("insert into {0}.{1} (samplename, default_genotype) values (?,?)".format(keyspaceName, sampleDefaultsTableName))
    variantInsertPrepStmt = session.prepare("insert into {0}.{1} (chrom,chunk,start_pos,ref,alt,samplename, sampleinfoformat, sampleinfo) values (?,?,?,?,?,?,?,?)".format(keyspaceName, variantTableName))

    returnErrMsg = None
    chosenFileToProcess = studyFileName
    try:
        baseDir = studyFilesInputDir
        os.chdir(baseDir)
        sampleNameCmd = subprocess.Popen("{0}/bin/bcftools query -l {1}".format(bcfToolsDir, chosenFileToProcess),
                                         shell=True, stdout=subprocess.PIPE)
        sampleName, err = sampleNameCmd.communicate()
        if err: return "Could not obtain sample name from the SNP file:{0}".format(chosenFileToProcess)
        sampleName = sampleName.strip()

        matchOutputFileName = sampleName + "_variantmatch.gz"
        errFileName = sampleName + "_variantmatch.err"
        bcfVariantMatchCmd = "({0}/bin/bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\t%LINE' -T {1} {3} | cut -f1,2,3,4,5,14,15 | gzip) 1> {2} 2> {4}".format(bcfToolsDir, variantPositionFileName, matchOutputFileName, chosenFileToProcess, errFileName)

        result = os.system(bcfVariantMatchCmd)
        errFileHandle = open(errFileName, "r")
        errlines = errFileHandle.readlines()
        errFileHandle.close()
        if not errlines:
            numMatchedVariants = getNumMatchedVariants(keyspaceName, sampleInsertLogTableName, studyName, sampleName)
            if (numTotVariants*1.0/numMatchedVariants) > 2:
                processVariantMatchFile(variantPositionFileName, matchOutputFileName, sampleName, "./.")
            else:
                processVariantMatchFile(variantPositionFileName, matchOutputFileName, sampleName, defaultRefGenotype)
        else:
            returnErrMsg = "Error in processing file:{0}".format(chosenFileToProcess) + os.linesep + os.linesep.join(errlines)
    except Exception, e:
        returnErrMsg = "Error in processing file:{0}".format(chosenFileToProcess) + os.linesep + traceback.format_exc()
    finally:
        cluster.shutdown()
        session.shutdown()
        if returnErrMsg: return returnErrMsg
        return None


if len(sys.argv) != 7:
    print("Usage: ProcessVariantMatches.py <Study Name> <Default Reference Genotype> <Full Path to study files> <Cassandra node IP1> <Cassandra node IP2> <BCF Tools Directory>")
    sys.exit(1)

# Parse arguments
studyName = sys.argv[1]
defaultRefGenotype = sys.argv[2]
studyFilesInputDir = sys.argv[3]
cassandraNodeIPs = [sys.argv[4], sys.argv[5]]
bcfToolsDir = sys.argv[6]

# Get the list of study files
os.chdir(studyFilesInputDir)
dirContents = glob.glob("*_filtervcf.gz")
dirContents.sort()
studyFileNames = dirContents

# Create keyspace, variant, header and sample insertlog tables in Cassandra
keyspaceName = "variant_ksp"
variantTableName = "variants_{0}".format(studyName.lower())
headerTableName = "headers_{0}".format(studyName.lower())
sampleInsertLogTableName = "sample_insert_log"
sampleDefaultsTableName = "sample_defaults_{0}".format(studyName.lower())
studyInfoTableName = "study_info_{0}".format(studyName.lower())
local_cluster_var = Cluster(cassandraNodeIPs)
local_session_var = local_cluster_var.connect()
local_session_var.execute("create table if not exists {0}.{1} (samplename varchar, default_genotype varchar, primary key(samplename));".format(keyspaceName,sampleDefaultsTableName))

conf = SparkConf().setMaster("spark://{0}:7077".format(get_ip_address())).setAppName("SingleSampleVCFMerge").set("spark.cassandra.connection.host", cassandraNodeIPs[0]).set("spark.scheduler.listenerbus.eventqueue.size", "100000").set("spark.cassandra.read.timeout_ms", 1200000).set("spark.cassandra.connection.timeout_ms", 1200000)
sc = SparkContext(conf=conf)
sc.setLogLevel("INFO")

variantPositionFileName = studyFilesInputDir + os.path.sep + "unique_variant_positions.gz"
rows = local_session_var.execute("select tot_num_variants from {0}.{1};".format(keyspaceName, studyInfoTableName))
numTotVariants = 0
if rows:
    rows = iter(rows)
    numTotVariants = rows.next().distinct_num_variants
else:
    raise Exception("Could not obtain number of variants for the study: {0} from the table: {1}.{2}".format(studyName, keyspaceName, studyInfoTableName))


if not os.path.isfile(variantPositionFileName):
    sql = SQLContext(sc)
    variants = sql.read.format("org.apache.spark.sql.cassandra").\
                   load(keyspace=keyspaceName, table=variantTableName)
    variants.registerTempTable("variantsTable")
    resultDF = sql.sql("select chrom,start_pos from variantsTable group by 1,2 order by 1,2")
    iterator = resultDF.toLocalIterator()

    variantPositionFileHandle = gzip.open(variantPositionFileName, "wb")
    for result in iterator:
        discardOutput = variantPositionFileHandle.write(result["chrom"] + "\t" + str(result["start_pos"]) + os.linesep)
    variantPositionFileHandle.close()

numPartitions = len(studyFileNames)
studyIndivRDD = sc.parallelize(studyFileNames, numPartitions)
processResults = studyIndivRDD.map(lambda studyFileName: matchVariantPosInStudyFiles(bcfToolsDir, studyName, defaultRefGenotype, studyFileName, numTotVariants, studyFilesInputDir, variantPositionFileName, cassandraNodeIPs, keyspaceName, sampleDefaultsTableName, variantTableName, sampleInsertLogTableName)).collect()
for result in processResults:
    if result:
        print(result)

local_session_var.shutdown()
local_cluster_var.shutdown()
sc.stop()