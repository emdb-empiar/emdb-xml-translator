#!/usr/bin/env python
"""
diff_all.py

Diff v1.9 files with v1.9 files translated from v3.0 files

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


def diff_all(v_19_dir, v_30_to_19_dir, out_dir):
    """
    Creates a txt file that contains any differences between
    a v1.9 (cannonical) file created from a v1.9 using emdb_xml_translate.py
    and
    a v1.9 file that is a result of a v1.9 -> v3.0 -> v1.9 roundtrip translations
    """
    command_list_base = ['diff', '-b']
    file_path_template = os.path.join(v_30_to_19_dir, 'emd-*.xml')
    emdb_files = glob.glob(file_path_template)
    num_errors = 0
    num_success = 0
    error_list = []
    for emdb_file in emdb_files:
        inf = os.path.basename(emdb_file)
        outf_base = inf[0:-4] + '.txt'
        outf = os.path.join(out_dir, outf_base)
        v_19_file = os.path.join(v_19_dir, inf)
        v_30_to_19_file = os.path.join(v_30_to_19_dir, inf)
        logging.info("v1.9 file: %s, v3.0 to 1.9 file: %s, output file: %s", v_19_file, v_30_to_19_file, outf)
        command_list = list(command_list_base)
        command_list.append(v_19_file)
        command_list.append(v_30_to_19_file)
        cmd_text = ' '.join(command_list)
        logging.info('Executing: %s', cmd_text)
        with open(outf, 'w') as out_f:
            exit_code = subprocess.call(command_list, stdout=out_f)
            if exit_code > 1:
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
    Diff v1.9 files with v1.9 files translated from v3.0 files
    """
    default_v_19_dir = EMDBSettings.emdb_19_to_19_dir_out
    default_v_30_to_19_dir = EMDBSettings.emdb_19_to_19_via_30relax_dir_out
    default_out_dir = EMDBSettings.diff_19s_dir

    # Handle command line options
    usage = """
            python diff_all.py [options]
            Diff v1.9 files with v1.9 files translated from 3.0 files

            Examples:
            python diff_all.py

            Typical run:
            python diff_all.py -i '/data/emdb19_to_19' -j '/data/emdb30relax_to_19' -o '/data/emdb_diff'
            /data/emdb19_to_19 is the input directory with the canonical EMDB 1.9 XML files (one part of the diff)
            /data/emdb30relax_to_19 is the input directory with EMDB 1.9 XML files converted from EMDB XML 3.0 (the other part of the diff)
            /data/emdb_diff is the output directory with the diff files with the same name as the entry but with the suffix .txt
            """
    version = "0.1"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-i", "--v19dir", action="store", type="string", metavar="DIR", dest="v19Dir", default=default_v_19_dir, help="Input directory with the canonical EMDB 1.9 XML files (one part of the diff) [default: %default]")
    parser.add_option("-j", "--v30to19dir", action="store", type="string", metavar="DIR", dest="v30To19Dir", default=default_v_30_to_19_dir,
                      help="Input directory with EMDB 1.9 XML files converted from EMDB XML 3.0 (the other part of the diff) [default: %default]")
    parser.add_option("-o", "--out-dir", action="store", type="string", metavar="DIR", dest="outDir", default=default_out_dir,
                      help="Output directory with the diff files with the same name as the entry but with the suffix .txt [default: %default]")
    (options, args) = parser.parse_args()
    print args
    diff_all(options.v19Dir, options.v30To19Dir, options.outDir)


if __name__ == "__main__":
    main()
