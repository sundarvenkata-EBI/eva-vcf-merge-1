# Copyright 2021 EMBL - European Bioinformatics Institute
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

import os


def are_all_elements_unique(elements):
    """Check if there are any repeated element in the list of elements. If yes return False otherwise return True."""
    unique_elements = set()
    for element in elements:
        if element in unique_elements:
            return False
        unique_elements.add(element)
    return True


def write_files_to_list(files, alias, output_dir):
    """Write the list of files to a path built from alias and output_dir."""
    list_filename = os.path.join(output_dir, f"{alias}_files.list")
    os.makedirs(os.path.dirname(list_filename), exist_ok=True)
    with open(list_filename, "w") as handle:
        for filename in files:
            handle.write(filename + "\n")
    return list_filename
