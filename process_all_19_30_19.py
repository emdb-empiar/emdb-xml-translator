#!/usr/bin/env python
"""
process_all_19_30_19.py

Wrapper script for emdb_xml_translate.py for converting existing v1.9 files into v3.0 files and back into v1.9

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


def process_all_19_30_19(file_path_template_v19, out_dir_19_to_30, out_dir_30_to_19):
    """
    Method for converting existing v1.9 files into v3.0 files and back into v1.9
    """
    # resulted v3.0 files should be validated with the relaxed schema (-v -r)
    command_list_19_to_30 = ['python', './emdb_xml_translate.py', '-v', '-r', '-p', '-f']
    command_list_30_to_19 = ['python', './emdb_xml_translate.py', '-v', '-r', '-p', '-i', '3.0', '-o', '1.9', '-f']
    emdb_19_files = glob.glob(file_path_template_v19)
    num_errors = 0
    num_success = 0
    error_list = []
    i = 0
    for emdb_19_file in emdb_19_files:
        i += 1
        if i > 0:
            #print i
            in_19_file = os.path.basename(emdb_19_file)
            if in_19_file == 'emd-0051.xml':
                print 'IN 1.9 FILE %s ' % emdb_19_file

                out_30_file = os.path.join(out_dir_19_to_30, in_19_file)
                print 'OUT 3.0 FILE %s ' % out_30_file

                command_list = list(command_list_19_to_30)
                command_list.append(out_30_file)
                command_list.append(emdb_19_file)
                cmd_text = ' '.join(command_list)
                print 'EXECUTING %s ' % cmd_text

                exit_code = subprocess.call(command_list)
                if exit_code != 0:
                    num_errors += 1
                    error_list.append(in_19_file)
                else:
                    num_success += 1

                # REVERSE
                in_30_file = out_30_file
                print 'IN 3.0 FILE %s' % in_30_file

                out_19_file = os.path.join(out_dir_30_to_19, in_19_file)
                print 'OUT 19 FILE %s' % out_19_file
                command_list = list(command_list_30_to_19)
                command_list.append(out_19_file)
                command_list.append(in_30_file)
                cmd_text = ' '.join(command_list)
                print 'EXECUTING %s ' % cmd_text
                exit_code = subprocess.call(command_list)

                if exit_code != 0:
                    num_errors += 1
                    error_list.append(in_30_file)
                else:
                    num_success += 1

def main():
    """
    Convert all EMDB XML 1.9 files to XML 3.0 relaxed header files and back to XML 1.9 files
    """
    file_path_template_v19 = EMDBSettings.header_19_19_template
    out_19_to_30_dir = EMDBSettings.emdb_19_to_30relax_dir_out
    out_30_to_19_dir = EMDBSettings.emdb_19_to_19_via_30relax_dir_out

    # Handle command line options
    usage = """
            python process_all_19_30_19.py [options]
            Convert all EMDB XML 3.0 header files to XML 1.9 files.

            Examples:
            python process_all_19_30_19.py

            """
    version = "0.1"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-t", "--template", action="store", type="string", metavar="TEMPLATE", dest="file_path_template_v19", default=file_path_template_v19, help="Template used to glob all input 1.9 header files [default: %default]")
    parser.add_option("-o", "--out-30-dir", action="store", type="string", metavar="DIR", dest="out_30_dir", default=out_19_to_30_dir, help="Directory for EMDB XML 3.0 relax files [default: %default]")
    parser.add_option("-f", "--final-out-19-dir", action="store", type="string", metavar="DIR", dest="out_19_dir", default=out_30_to_19_dir, help="Directory for EMDB XML 1.9 files [default: %default]")
    (options, args) = parser.parse_args()
    process_all_19_30_19(options.file_path_template_v19, options.out_30_dir, options.out_19_dir)


if __name__ == "__main__":
    main()
