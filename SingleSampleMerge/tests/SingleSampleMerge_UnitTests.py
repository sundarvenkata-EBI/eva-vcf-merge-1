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
        bcfToolsDir = os.getenv("HOME") + os.path.sep + "bcftools"

        studyFilesInputDir = tempfile.mkdtemp()
        vcfFileName = studyFilesInputDir + os.path.sep + "Proband-396.vcf.gz"
        downloadVCFFileFromFTP("ftp://ftp.ebi.ac.uk/pub/databases/eva/PRJEB21300/submitted_files", studyFilesInputDir, os.path.basename(vcfFileName))
        self.assertEqual(SSMergeCommonUtils.getSampleName(bcfToolsDir, vcfFileName), "Proband-396")

        studyFilesInputDir = tempfile.mkdtemp()
        vcfFileName = studyFilesInputDir + os.path.sep + "P1-CA1.combined.sorted.vcf.gz"
        downloadVCFFileFromFTP("ftp://ftp.ebi.ac.uk/pub/databases/eva/PRJEB10956/submitted_files", studyFilesInputDir,
                               os.path.basename(vcfFileName))
        self.assertEqual(SSMergeCommonUtils.getSampleName(bcfToolsDir, vcfFileName), "2")

    def test_writeVariantToCassandra(self):
        global cassandraCluster, cassandraSession, variantInsertStmt, headerInsertStmt, uniquePosInsertStmt
        keyspaceName = "variant_ksp"
        randomSuffix = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
        variantTableName = "variant_test" + randomSuffix
        uniquePosTableName = "uniq_pos_" + randomSuffix
        sampleName = randomSuffix

        try:
            # region Setup Cassandra driver settings
            cassandraCluster = Cluster("172.22.70.148,172.22.70.139,172.22.70.150,172.22.68.19".split(","))
            cassandraSession = cassandraCluster.connect()
            cassandraSession.execute("create table if not exists {0}.{1} "
                                      "(chrom varchar, chunk int, start_pos bigint, ref varchar, alt varchar, qual varchar, "
                                      "filter varchar, info varchar, sampleinfoformat varchar, "
                                      "sampleinfo  varchar, var_id varchar, var_uniq_id varchar, sampleName varchar,  "
                                      "primary key((chrom, chunk), start_pos, var_uniq_id, samplename));"
                                      .format(keyspaceName, variantTableName))
            cassandraSession.execute("create table if not exists {0}.{1} "
                                      "(chrom varchar, chunk int, start_pos bigint, "
                                      "primary key((chrom,chunk), start_pos));"
                                      .format(keyspaceName, uniquePosTableName))

            variantInsertStmt = variantInsertStmt = cassandraSession.prepare(
                "insert into {0}.{1} "
                "(chrom,chunk,start_pos,ref,alt,qual,filter,info,sampleinfoformat,sampleinfo,var_id,var_uniq_id,sampleName) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)".format(keyspaceName, variantTableName))
            uniquePosInsertStmt = cassandraSession.prepare(
                "insert into {0}.{1} (chrom, chunk, start_pos) VALUES (?,?,?)"
                .format(keyspaceName, uniquePosTableName))
            # endregion

            # region Write lines to cassandra
            linesToWrite = []
            expectedResults = []

            line = "1\t1000\trs123\tTCACCC\tTGACGG\t.\t.\t."
            chromosome, position, rsID, ref, alt, qual, qualFilter, info = line.split("\t")
            linesToWrite.append(line)
            variantID = chromosome + "_" + str(position).zfill(12) + "_" + hashlib.md5(
                ref + "_" + alt).hexdigest() + "_" + sampleName.zfill(20)
            expectedResults.append(line + "\t\t\t0\t{0}\t{1}".format(variantID, sampleName))

            line = "1\t999\trs123\tGTCACCC\tG\t105\tFAIL\tAN=3;AC=1,2;AF=0.125,0.25;DP=63;NS=4;MQ=10685\tGT\t0|0"
            chromosome,position,rsID,ref,alt,qual,qualFilter,info,sampleInfoFormat,sampleInfo = line.split("\t")
            variantID = chromosome + "_" + str(position).zfill(12) + "_" + hashlib.md5(
                ref + "_" + alt).hexdigest() + "_" + sampleName.zfill(20)
            linesToWrite.append(line)
            expectedResults.append(line + "\t0\t{0}\t{1}".format(variantID, sampleName))

            line = "1\t1000040\trs123\tTGACGTAACGATT\tT,TGACGTAACGGTT,TGACGTAATAC\t110\tPASS\t.\tGT:AD:DP\t0/0:17:94"
            chromosome, position, rsID, ref, alt, qual, qualFilter, info, sampleInfoFormat, sampleInfo = line.split(
                "\t")
            variantID = chromosome + "_" + str(position).zfill(12) + "_" + hashlib.md5(
                ref + "_" + alt).hexdigest() + "_" + sampleName.zfill(20)
            linesToWrite.append(line)
            expectedResults.append(line + "\t0\t{0}\t{1}".format(variantID, sampleName))

            # endregion

        except Exception, ex:
            traceback.print_exc()
        finally:
            if keyspaceName and variantTableName:
                pass
                # cassandraSession.execute("drop table {0}.{1}".format(keyspaceName, variantTableName))


if __name__ == '__main__':
    unittest.main()