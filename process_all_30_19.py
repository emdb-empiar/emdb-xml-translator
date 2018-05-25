#!/usr/bin/env python
"""
process_all_30_19.py

Wrapper script for emdb_xml_translate.py for converting v3.0 files in EMDB to v1.9

TODO:

Version history:

Copyright [2014-2016] EMBL - European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the
"License"); you may not use this file except in
compliance with the License. You may obtain a copy of
the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on
an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
"""
import glob
import os
import logging
import subprocess
from optparse import OptionParser
from emdb_settings import EMDBSettings

__author__ = 'Ardan Patwardhan, Sanja Abbott'
__email__ = 'ardan@ebi.ac.uk, sanja@ebi.ac.uk'
__date__ = '2017-07-27'

logging.basicConfig(level=EMDBSettings.log_level, format=EMDBSettings.log_format)


def process_all_30_19(file_path_template, out_dir):
    """
    TO DO
    """
    command_list_base = ['python', './emdb_xml_translate.py', '-v', '-i', '3.0', '-o', '1.9', '-f']
    emdb_files = glob.glob(file_path_template)
    num_errors = 0
    num_success = 0
    error_list = []
    i = 0
    for emdb_file in emdb_files:
        i += 1
        if i > 0:
            print i
            inf = os.path.basename(emdb_file)
            print inf
            outf = os.path.join(out_dir, inf)
            #logging.info("Input file: %s, output file: %s", emdb_file, outf)
            command_list = list(command_list_base)
            command_list.append(outf)
            command_list.append(emdb_file)
            cmd_text = ' '.join(command_list)
            #logging.info('Executing: %s', cmd_text)
            exit_code = subprocess.call(command_list)
            if exit_code != 0:
                num_errors += 1
                error_list.append(inf)
            else:
                num_success += 1
            # logging.warning('%d files successfully processed!', num_success)
            if num_errors > 0:
                logging.warning('%d errors!', num_errors)
                logging.warning('List of entries that were not translated')
                for entry in error_list:
                    logging.warning(entry)


def main():
    """
    Convert all EMDB XML 3.0 header files to XML 1.9 files
    """
    default_file_path_template = EMDBSettings.header_30_template
    default_out_dir = EMDBSettings.emdb_30_to_19_dir_out

    # Handle command line options
    usage = """
            python process_all_30_19.py [options]
            Convert all EMDB XML 3.0 header files to XML 1.9 files.

            Examples:
            python process_all_30_19.py

            Typical run:
            python process_all_30_19.py -t 'data/input/v3.0/EMD-*.xml' -o 'data/output/emdb30_to_19'
            data/input/v3.0/EMD-*.xml is the template used to glob all input 3.0 header files
            data/output/emdb30_to_19 is the output directory with the EMDB XML 1.9 files
            """
    version = "0.1"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-t", "--template", action="store", type="string", metavar="TEMPLATE", dest="file_path_template", default=default_file_path_template, help="Template used to glob all input 3.0 header files [default: %default]")
    parser.add_option("-o", "--out-dir", action="store", type="string", metavar="DIR", dest="out_dir", default=default_out_dir, help="Directory for EMDB XML 1.9 files [default: %default]")
    (options, args) = parser.parse_args()
    process_all_30_19(options.file_path_template, options.out_dir)


if __name__ == "__main__":
    main()
