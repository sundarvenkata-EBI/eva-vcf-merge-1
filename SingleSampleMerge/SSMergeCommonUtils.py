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

import socket

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

# Variables

# Chunk size by which variant table is partitioned in Cassandra.
# Determines which variants are written to which nodes based on their chromosome and position.
# For ex: variant at chr 1 and position 1M will be written to a different node than the variant at chr 2 and position 2M.
CHR_POS_CHUNKSIZE = int(1e6)