import unittest, tempfile, SingleSampleMerge, ftplib, os, re, random, string, hashlib, traceback, SSMergeCommonUtils
from cassandra.cluster import Cluster


def downloadVCFFileFromFTP(ftpDir, studyFilesInputDir, vcfFileName):
    currDir = os.getcwd()
    try:
        os.chdir(studyFilesInputDir)
        exp = re.compile("//|/")
        topLevelFTPSite = re.split(exp, ftpDir)[1]
        ftpDir = "/" + "/".join(re.split(exp, ftpDir)[2:])
        ftp = ftplib.FTP(topLevelFTPSite)
        ftp.login()
        ftp.cwd(ftpDir)
        ftp.retrbinary("RETR {0}".format(vcfFileName), open(os.path.basename(vcfFileName), "wb").write)
    except Exception, ex:
        raise Exception(ex.message)
    finally:
        os.chdir(currDir)
        return True

class TestStringDiffIndex(unittest.TestCase):

    def setUp(self):
        self.ssMergeObj = self.createSingleSampleMergeObj()

    def test_stringDiffIndex(self):
        """
        For test cases, see https://commons.apache.org/proper/commons-lang/apidocs/org/apache/commons/lang3/StringUtils.html#indexOfDifference-java.lang.CharSequence...-
        """
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex(None, "something"), -1)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("something", None), -1)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("", ""), -1)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("", None), -1)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("", "abc"), 0)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("abc", ""), 0)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("abc", "abc"), -1)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("abc", "a"), 1)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("ab", "abxyz"), 2)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("abcde", "abxyz"), 2)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("abcde", "xyz"), 0)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("xyz","abcde"), 0)
        self.assertEqual(SSMergeCommonUtils.stringDiffIndex("i am a machine", "i am a robot"), 7)

    def test_getNormalizedStartPos(self):
        """
        For test cases, see https://github.com/EBIvariation/eva-pipeline/blob/develop/src/test/java/uk/ac/ebi/eva/pipeline/io/mappers/VariantVcfFactoryTest.java         
        """
        chromosome, start_pos, ref, alt = "1",1000,"",""
        with self.assertRaises(Exception):
            SSMergeCommonUtils.getNormalizedStartPos(chromosome, start_pos, ref, alt)
        chromosome, start_pos, ref, alt = "1", 1000, "CATAT", "CAT"
        self.assertEqual(SSMergeCommonUtils.getNormalizedStartPos(chromosome, start_pos, ref, alt), 1001)
        chromosome, start_pos, ref, alt = "1", 1000, "GTA", "GT"
        self.assertEqual(SSMergeCommonUtils.getNormalizedStartPos(chromosome, start_pos, ref, alt), 1002)
        chromosome, start_pos, ref, alt = "1", 1000, "AC", "ATC"
        self.assertEqual(SSMergeCommonUtils.getNormalizedStartPos(chromosome, start_pos, ref, alt), 1001)
        chromosome, start_pos, ref, alt = "1", 1000, "ATC", "AC"
        self.assertEqual(SSMergeCommonUtils.getNormalizedStartPos(chromosome, start_pos, ref, alt), 1001)
        chromosome, start_pos, ref, alt = "1", 1000000000, "AC", "ACT"
        self.assertEqual(SSMergeCommonUtils.getNormalizedStartPos(chromosome, start_pos, ref, alt), 1000000002)
        chromosome, start_pos, ref, alt = "1", 10040, "TGACGTAACGATT", "TGACGTAACGGTT"
        self.assertEqual(SSMergeCommonUtils.getNormalizedStartPos(chromosome, start_pos, ref, alt), 10050)
        chromosome, start_pos, ref, alt = "1", 10040, "TGACGTAACGATT", "TGACGTAATAC"
        self.assertEqual(SSMergeCommonUtils.getNormalizedStartPos(chromosome, start_pos, ref, alt), 10048)

    def test_getSampleName(self):
        studyFilesInputDir = tempfile.mkdtemp()
        vcfFileName = studyFilesInputDir + os.path.sep + "Proband-396.vcf.gz"
        downloadVCFFileFromFTP("ftp://ftp.ebi.ac.uk/pub/databases/eva/PRJEB21300/submitted_files", studyFilesInputDir, os.path.basename(vcfFileName))
        self.assertEqual(SSMergeCommonUtils.getSampleName(self.ssMergeObj.bcfToolsDir, vcfFileName), "Proband-396")

        studyFilesInputDir = tempfile.mkdtemp()
        vcfFileName = studyFilesInputDir + os.path.sep + "P1-CA1.combined.sorted.vcf.gz"
        downloadVCFFileFromFTP("ftp://ftp.ebi.ac.uk/pub/databases/eva/PRJEB10956/submitted_files", studyFilesInputDir,
                               os.path.basename(vcfFileName))
        self.assertEqual(SSMergeCommonUtils.getSampleName(self.ssMergeObj.bcfToolsDir, vcfFileName), "2")

    def test_writeVariantToCassandra(self):

        # region Write lines to cassandra
        linesToWrite = []
        expectedVariantResults = {}
        expectedPositionResults = set()

        line = "1\t2000\trs123\tTCACCC\tTGACGG\t.\t.\t."
        chromosome, position, rsID, ref, alt, qual, qualFilter, info = line.split("\t")
        linesToWrite.append(line)
        variantID = chromosome + "_" + str(position).zfill(12) + "_" + hashlib.md5(
            ref + "_" + alt).hexdigest() + "_" + self.ssMergeObj.sampleName.zfill(20)
        expectedVariantResults[variantID] = line + "\t\t\t0\t{0}\t{1}".format(variantID, self.ssMergeObj.sampleName)
        expectedPositionResults.add("1\t0\t2000")
        expectedPositionResults.add("1\t0\t2001")


        line = "1\t999\trs123\tGTCACCC\tG\t105\tFAIL\tAN=3;AC=1,2;AF=0.125,0.25;DP=63;NS=4;MQ=10685\tGT\t0|0"
        chromosome,position,rsID,ref,alt,qual,qualFilter,info,sampleInfoFormat,sampleInfo = line.split("\t")
        variantID = chromosome + "_" + str(position).zfill(12) + "_" + hashlib.md5(
            ref + "_" + alt).hexdigest() + "_" + self.ssMergeObj.sampleName.zfill(20)
        linesToWrite.append(line)
        expectedVariantResults[variantID] =  line + "\t0\t{0}\t{1}".format(variantID, self.ssMergeObj.sampleName)
        expectedPositionResults.add("1\t0\t999")
        expectedPositionResults.add("1\t0\t1000")

        line = "1\t1000040\trs123\tTGACGTAACGATT\tT,TGACGTAACGGTT,TGACGTAATAC\t110\tPASS\t.\tGT:AD:DP\t0/0:17:94"
        chromosome, position, rsID, ref, alt, qual, qualFilter, info, sampleInfoFormat, sampleInfo = line.split(
            "\t")
        variantID = chromosome + "_" + str(position).zfill(12) + "_" + hashlib.md5(
            ref + "_" + alt).hexdigest() + "_" + self.ssMergeObj.sampleName.zfill(20)
        linesToWrite.append(line)
        expectedVariantResults[variantID] = line + "\t1\t{0}\t{1}".format(variantID, self.ssMergeObj.sampleName)
        expectedPositionResults.add("1\t1\t1000040")
        expectedPositionResults.add("1\t1\t1000048")
        expectedPositionResults.add("1\t1\t1000050")

        self.ssMergeObj.writeVariantToCassandra(linesToWrite)
        # endregion

        resultRows = self.ssMergeObj.cassandraSession.execute("select * from {0}.{1}"
                                                              .format(self.ssMergeObj.keyspaceName,
                                                                      self.ssMergeObj.variantTableName))
        if resultRows:
            resultRows = iter(resultRows)
            for result in resultRows:
                resultString = "\t".join([str(elem) if elem is not None else '' for elem in
                                          [result.chrom, result.start_pos, result.var_id, result.ref, result.alt,
                                          result.qual, result.filter, result.info, result.sampleinfoformat,
                                          result.sampleinfo, result.chunk, result.var_uniq_id, result.samplename]])
                variantID = result.var_uniq_id
                self.assertEqual(resultString, expectedVariantResults[variantID])
        else:
            self.fail("No result rows")

        resultRows = self.ssMergeObj.cassandraSession.execute("select * from {0}.{1}"
                                                              .format(self.ssMergeObj.keyspaceName,
                                                                      self.ssMergeObj.uniquePosTableName))
        if resultRows:
            resultRows = iter(resultRows)
            for result in resultRows:
                resultString = "\t".join([str(elem) if elem is not None else '' for elem in
                                          [result.chrom, result.chunk, result.start_pos]])
                self.assertIn(resultString, expectedPositionResults)
        else:
            self.fail("No result rows")


    def createSingleSampleMergeObj(self):
        bcfToolsDir = os.getenv("HOME") + os.path.sep + "bcftools"
        keyspaceName = "variant_ksp"
        randomSuffix = "test"
        variantTableName = "variant_" + randomSuffix
        uniquePosTableName = "uniq_pos_" + randomSuffix
        headerTableName = "header_" + randomSuffix
        sampleInsertLogTableName = "sample_insert_log"
        cassandraNodeIPs = "172.22.70.148,172.22.70.139,172.22.70.150,172.22.68.19".split(",")
        studyName = "PRJEBTEST"
        studyVCFFileName = os.getenv("HOME") + os.path.sep + "mergevcfinput/PRJEBTEST" + os.path.sep + "Proband-1152.vcf.gz"
        studyInfoTableName = "study_info_{0}".format(studyName.lower())

        local_cluster_var = Cluster(cassandraNodeIPs)
        local_session_var = local_cluster_var.connect()
        local_session_var.execute("create keyspace if not exists {0}".format(keyspaceName) +
                                  " with replication = {'class': 'SimpleStrategy', 'replication_factor' : 1};")
        local_session_var.execute("create table if not exists {0}.{1} "
                                  "(samplename varchar, header varchar, "
                                  "primary key(samplename));"
                                  .format(keyspaceName, headerTableName))
        local_session_var.execute("create table if not exists {0}.{1} "
                                  "(chrom varchar, chunk int, start_pos bigint, ref varchar, alt varchar, qual varchar, "
                                  "filter varchar, info varchar, sampleinfoformat varchar, "
                                  "sampleinfo  varchar, var_id varchar, var_uniq_id varchar, sampleName varchar,  "
                                  "primary key((chrom, chunk), start_pos, var_uniq_id, samplename));"
                                  .format(keyspaceName, variantTableName))
        local_session_var.execute("create table if not exists {0}.{1} "
                                  "(studyname varchar, samplename varchar, proc_status varchar, num_variants bigint, "
                                  "primary key((studyname, samplename)));"
                                  .format(keyspaceName, sampleInsertLogTableName))
        local_session_var.execute("create table if not exists {0}.{1} "
                                  "(studyname varchar, tot_num_variants bigint, distinct_num_variant_pos bigint, "
                                  "primary key(studyname));"
                                  .format(keyspaceName, studyInfoTableName))
        local_session_var.execute("create table if not exists {0}.{1} "
                                  "(chrom varchar, chunk int, start_pos bigint, "
                                  "primary key((chrom,chunk), start_pos));"
                                  .format(keyspaceName, uniquePosTableName))

        for tableName in [headerTableName, variantTableName, studyInfoTableName, uniquePosTableName]:
            local_session_var.execute("truncate {0}.{1}".format(keyspaceName, tableName))

        local_session_var.shutdown()
        local_cluster_var.shutdown()

        return SingleSampleMerge.SingleSampleMerge(bcfToolsDir, studyName, studyVCFFileName, cassandraNodeIPs, keyspaceName,
                                            headerTableName,variantTableName, sampleInsertLogTableName,
                                            uniquePosTableName)



if __name__ == '__main__':
    unittest.main()