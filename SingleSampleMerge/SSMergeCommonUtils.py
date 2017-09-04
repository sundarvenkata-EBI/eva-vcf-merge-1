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
import hashlib
import socket, glob, os, subprocess

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
    sampleNameCmd = subprocess.Popen("{0}/bin/bcftools query -l {1}".format(bcfToolsDir, vcfFileName), shell=True,
                                     stdout=subprocess.PIPE)
    sampleName, err = sampleNameCmd.communicate()
    if err: return err
    return sampleName.strip()


def getNormalizedStartPos(chromosome, start_pos, ref, alt):
    """
    Get a normalized start position for a given variant.
    This is almost a verbatim implementation of https://github.com/EBIvariation/eva-pipeline/blob/develop/src/main/java/uk/ac/ebi/eva/pipeline/io/mappers/VariantVcfFactory.java#L190 
    but removes any portions not related to calculating the normalized start position.

    :param chromosome: Variant chromosome
    :type chromosome: str
    :param start_pos: Variant start position
    :type start_pos: long
    :param ref: Variant reference allele
    :type ref: str
    :param alt: Variant alternate allele
    :type alt: str
    :return: The normalized start position
    :rtype: long
    """
    if (ref == alt):
        raise Exception("Alternate allele is identical to the reference in:{0}"
                        .format(",".join([chromosome, str(start_pos), ref, alt])))

    # Remove trailing bases
    refReversed = ref[::-1]
    altReversed = alt[::-1]
    indexOfDifference = stringDiffIndex(refReversed, altReversed)
    ref = refReversed[indexOfDifference:][::-1]
    alt = altReversed[indexOfDifference:][::-1]
    # Remove leading bases
    indexOfDifference = stringDiffIndex(ref, alt)
    start = start_pos + indexOfDifference
    return start

def getTotDistinctVarPosFromSampleFiles(vcfInputDirectory):
    """
    Get the total number of distinct variant positions from all the VCF files across all the samples

    :param vcfInputDirectory: Full path to the VCF directory
    :type vcfInputDirectory: str
    :return: Total number of distinct variants
    :rtype: long
    """
    procResult = subprocess.Popen(
        'zcat {0}/*_filtervcf.gz | grep -v "^#" | cut -f1,2 | sort -T {0} | uniq | wc -l'.format(vcfInputDirectory),
        shell=True,
        stdout=subprocess.PIPE)
    totNumDistinctVarPos, errMsg = procResult.communicate()
    if errMsg:
        raise Exception("Could not get distinct number of variants due to error:{0}".format(errMsg))
    else:
        return long(float(totNumDistinctVarPos))


def getVCFFileNames(studyFilesInputDir):
    """
    Given an input directory with all the submitted files for a study, return just the VCF files

    :param studyFilesInputDir: Input directory containing the submitted files
    :type studyFilesInputDir: str
    :return: The list of VCF files in the study directory
    :rtype: list[str]
    """
    dirContents = glob.glob(studyFilesInputDir + os.path.sep + "*.vcf.gz")
    dirContents.sort()
    return dirContents

def get_ip_address():
    """
    Get the current nodes' IP address
    
    :return: IP address of the current node
    :rtype: str
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    return s.getsockname()[0]

def getErrFileContents(errFileName):
    """
    Get contents, if present, from an error file
    
    :param errFileName: Full path to the error file
    :return: Contents of the error file if there were any errors, None otherwise
    :rtype: list
    """
    with open(errFileName, "r") as errFileHandle:
        errFileContents = errFileHandle.readlines()
        if errFileContents: return errFileContents
        return []

def getUniqueVariantID(sampleName, chromosome, position, ref, alt):
    variantID = chromosome + "_" + str(position).zfill(12) + "_" + hashlib.md5(ref + "_" + alt).hexdigest() \
                + "_" + sampleName.zfill(20)
    return variantID


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

# Variables

# Chunk size by which variant table is partitioned in Cassandra.
# Determines which variants are written to which nodes based on their chromosome and position.
# For ex: variant at chr 1 and position 1M will be written to a different node than the variant at chr 2 and position 2M.
CHR_POS_CHUNKSIZE = int(1e6)

SMALL_QUERY_TIMEOUT_IN_SECS = 120
LARGE_QUERY_TIMEOUT_IN_SECS = 12000
BATCH_WRITE_TIMEOUT_IN_SECS = 1200
SPARK_CASSANDRA_READ_TIMEOUT_IN_MS = 1200000