#!/usr/bin/env python
"""
process_all_19_30.py

Wrapper script for emdb_xml_translate.py for converting v1.9 files
in EMDB to v3.0.

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
__date__ = '2017-06-14'

logging.basicConfig(level=EMDBSettings.log_level, format=EMDBSettings.log_format)


def process_all_19_30(file_path_template, out_dir):
    """
    TO DO
    """
    command_list_base = ['python', './emdb_xml_translate.py', '-f']
    emdb_files = glob.glob(file_path_template)
    num_errors = 0
    num_success = 0
    error_list = []
    i = 1
    for emdb_file in emdb_files:
        print i
        i = i + 1
        inf = os.path.basename(emdb_file)
        outf = os.path.join(out_dir, inf)
        logging.info("Input file: %s, output file: %s", emdb_file, outf)
        command_list = list(command_list_base)
        command_list.append(outf)
        command_list.append(emdb_file)
        cmd_text = ' '.join(command_list)
        logging.info('Executing: %s', cmd_text)
        exit_code = subprocess.call(command_list)
        if exit_code != 0:
            num_errors += 1
            error_list.append(inf)
        else:
            num_success += 1
    logging.warning('%d files successfully processed!', num_success)
    if num_errors > 0:
        logging.warning('%d errors!', num_errors)
        logging.warning('List of entries that were not translated')
        for entry in error_list:
            logging.warning(entry)


def main():
    """
    Convert all EMDB XML 1.9 header files to XML 3.0 files
    """
    default_file_path_template = EMDBSettings.header_19_from_30_template
    default_out_dir = EMDBSettings.emdb_19_to_30_dir_out

    # Handle command line options
    usage = """
            process_all_19_30.py [options]
            Convert all EMDB XML 1.9 header files to XML 3.x.x.x files.

            Examples:
            python process_all_19_30.py

            Typical run:
            python process_all_19_30.py -t '/data/emstaging/EMD-*/header/emd-*.xml' -o '/data/emdb30'
            /data/emstaging/EMD-*/header/emd-*.xml is the template used to glob all input 1.9 header files
            /data/emdb20 is the output directory with the EMDB XML 3.0 files

            python process_all_19_30.py -t '../emdb_xml_translator/emdb_xml_translator/data/output/v2.0_strict_to_1.9/EMD-*.xml' -o '../emdb_xml_translator/emdb_xml_translator/data/output/tmp/'

            """
    version = "0.1"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-t", "--template", action="store", type="string", metavar="TEMPLATE", dest="filePathTemplate", default=default_file_path_template, help="Template used to glob all input 1.9 header files [default: %default]")
    parser.add_option("-o", "--out-dir", action="store", type="string", metavar="DIR", dest="outDir", default=default_out_dir, help="Directory for EMDB XML 2.0 files [default: %default]")
    (options, args) = parser.parse_args()
    process_all_19_30(options.filePathTemplate, options.outDir)


if __name__ == "__main__":
    main()
