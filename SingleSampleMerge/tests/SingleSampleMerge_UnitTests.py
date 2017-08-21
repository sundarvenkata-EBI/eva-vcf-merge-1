import unittest, tempfile, SingleSampleMerge, ftplib, os, gzip

class TestStringDiffIndex(unittest.TestCase):
    def downloadVCFFileFromFTP(studyFilesInputDir, vcfFileName):
        currDir = os.getcwd()
        try:
            os.chdir(studyFilesInputDir)
            ftp = ftplib.FTP("ftp.ebi.ac.uk")
            ftp.login()
            ftp.cwd("/pub/databases/eva/PRJEB21300/submitted_files/")
            ftp.retrbinary("RETR {0}".format(vcfFileName), gzip.open("Proband-396.vcf.gz","wb").write)
        except Exception, ex:
            raise Exception(ex.message)
        finally:
            return True
            os.chdir(currDir)


    def test_stringDiffIndex(self):
        """
        For test cases, see https://commons.apache.org/proper/commons-lang/apidocs/org/apache/commons/lang3/StringUtils.html#indexOfDifference-java.lang.CharSequence...-
        """
        self.assertEqual(SingleSampleMerge.stringDiffIndex(None, "something"), -1)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("something", None), -1)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("", ""), -1)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("", None), -1)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("", "abc"), 0)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("abc", ""), 0)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("abc", "abc"), -1)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("abc", "a"), 1)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("ab", "abxyz"), 2)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("abcde", "abxyz"), 2)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("abcde", "xyz"), 0)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("xyz","abcde"), 0)
        self.assertEqual(SingleSampleMerge.stringDiffIndex("i am a machine", "i am a robot"), 7)

    def test_getNormalizedStartPos(self):
        """
        For test cases, see https://github.com/EBIvariation/eva-pipeline/blob/develop/src/test/java/uk/ac/ebi/eva/pipeline/io/mappers/VariantVcfFactoryTest.java         
        """
        chromosome, start_pos, ref, alt = "1",1000,"",""
        with self.assertRaises(Exception):
            SingleSampleMerge.getNormalizedStartPos(chromosome, start_pos, ref, alt)
        chromosome, start_pos, ref, alt = "1", 1000, "CATAT", "CAT"
        self.assertEqual(SingleSampleMerge.getNormalizedStartPos(chromosome, start_pos, ref, alt), 1001)
        chromosome, start_pos, ref, alt = "1", 1000, "GTA", "GT"
        self.assertEqual(SingleSampleMerge.getNormalizedStartPos(chromosome, start_pos, ref, alt), 1002)
        chromosome, start_pos, ref, alt = "1", 1000, "AC", "ATC"
        self.assertEqual(SingleSampleMerge.getNormalizedStartPos(chromosome, start_pos, ref, alt), 1001)
        chromosome, start_pos, ref, alt = "1", 1000, "ATC", "AC"
        self.assertEqual(SingleSampleMerge.getNormalizedStartPos(chromosome, start_pos, ref, alt), 1001)
        chromosome, start_pos, ref, alt = "1", 1000000000, "AC", "ACT"
        self.assertEqual(SingleSampleMerge.getNormalizedStartPos(chromosome, start_pos, ref, alt), 1000000002)
        chromosome, start_pos, ref, alt = "1", 10040, "TGACGTAACGATT", "TGACGTAACGGTT"
        self.assertEqual(SingleSampleMerge.getNormalizedStartPos(chromosome, start_pos, ref, alt), 10050)
        chromosome, start_pos, ref, alt = "1", 10040, "TGACGTAACGATT", "TGACGTAATAC"
        self.assertEqual(SingleSampleMerge.getNormalizedStartPos(chromosome, start_pos, ref, alt), 10048)

    def test_getSampleName(self):
        studyFilesInputDir = tempfile.mkdtemp()

        try:
            os.chdir(studyFilesInputDir)
            ftp = ftplib.FTP("ftp.ebi.ac.uk")
            ftp.login()
            ftp.cwd("/pub/databases/eva/PRJEB21300/submitted_files/")
            ftp.retrbinary("RETR Proband-396.vcf.gz", gzip.open("Proband-396.vcf.gz","wb").write)
        except Exception, ex:
            raise Exception(ex.message)
        finally:
            self.assertEqual(SingleSampleMerge.getSampleName())
            os.chdir(currDir)




if __name__ == '__main__':
    unittest.main()