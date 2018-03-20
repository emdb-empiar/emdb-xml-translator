#!/usr/bin/env python
"""
process_all_19_30_19.py

Wrapper script for emdb_xml_translate.py for converting v 1.9 files
in EMDB to v1.9. The only reason for doing this is that it puts elements
in a canonical order that makes comparison with 1.9 -> 3.0 -> 1.9
translation easier.

TODO:

Version history:
0.2, 2015-11-12, Ardan Patwardhan: Minor changes associated with moving file to new project structure
0.3, 2015-11-18, Ardan Patwardhan: Adding mechanism to exclude empty tags which should make comparison easier
0.4, 2017-11-14, Sanja Abbott: Added mechanism to detect open and close tags


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
import re
import logging
import subprocess
from dateutil import parser as dtp
from optparse import OptionParser
from emdb_settings import EMDBSettings

__author__ = 'Ardan Patwardhan, Sanja Abbott'
__email__ = 'ardan@ebi.ac.uk, sanja@ebi.ac.uk'
__date__ = '2017-06-14'

logging.basicConfig(level=EMDBSettings.log_level, format=EMDBSettings.log_format)


def process_all_19_19(file_path_template, out_dir):
    """
    Take a v1.9 file and read and write it using emdb_xml_translate.py to put it in a canonical form.
    Some post processing is also done to remove empty tags etc

    Parameters:
    @param file_path_template: Regular expression that is passed to a glob function to extract a list of input files
    @param out_dir: The canonical files will be written to this directory
    """
    command_list_base = ['python', './emdb_xml_translate.py', '-i', '1.9', '-o', '1.9', '-f']
    pat = re.compile(r'^[ ]*<(sliceSet|fscSet|maskSet|figureSet|fitting|externalReferences|pdbEntryIdList|imageAcquisition|specimenPreparation|helicalParameters)/>[ ]*$')
    space_tags = 'title|articleTitle|software|resolutionMethod|algorithm|name|timeResolvedState|molWtMethod|details|date|file|contourLevel|targetCriteria'
    pat_spaces_start = re.compile(r'^(.*<(%s)>)(((?!</(%s)>).)*)((</(%s)>)?(((?!<(/%s>)).)*))' % (space_tags, space_tags, space_tags, space_tags))
    pat_spaces_end = re.compile(r'^(((?!</(%s>)).)*)(</(%s)>.*)' % (space_tags, space_tags))
    pat_spaces_sub = re.compile(r'\s+')
    pat_spaces_rep = ' '

    open_pat = re.compile(r'^[ ]*<(supplement|fitting)>[ ]*$')
    close_pat = re.compile(r'^[ ]*</(supplement|fitting)>[ ]*$')

    open_close_pat = re.compile(r'^[ ]*<(volume|firstPage|lastPage|year)></(volume|firstPage|lastPage|year)>[ ]*$')
    open_close_details_pat = re.compile(r'^[ ]*<details></details>[ ]*$')

    def clean_date(date_content):
        """
        """
        return_date = None
        known_date_issues = {'2001-19-12': '2001-12-19', '1997-': '1997-01-01', 'None': '2000-01-01'}
        if date_content in known_date_issues:
            #print 'known'
            date_corr = known_date_issues.get(date_content)
            #print date_corr
            return_date = dtp.parse(date_corr).date()
        else:
            try:
                return_date = dtp.parse(date_content).date()
            except:
                pass
            #print 'input date %s' % date_content
        #print 'return_date %s ' % return_date
        return str(return_date)

    def open_close_tags_on_next_line(open_pat, close_pat, line, index, lines):
        """
        """
        return_value = False
        open_match = re.match(open_pat, line)
        close_match = None
        if index + 1 < len(lines):
            next_line = lines[index + 1]
            close_match = re.match(close_pat, next_line)
        if open_match is not None and close_match is not None:
            return_value = True

        return return_value

    def clean_spaces(inf_hd, outf_hd):
        """
        Reduce spaces, new-lines and tabs to single spaces
        Note: will NOT handle nested tags!!

        Parameters:
        @param inf_hd: Input file handle
        @param outf_hd: Output file handle
        """
        inf_hd.seek(0)
        outf_hd.seek(0)
        start_tag_found = False
        lines = inf_hd.readlines()
        index = -1
        while index < len(lines) - 1:
            index = index + 1
            line = lines[index]
            if start_tag_found is False:
                line_match = re.match(pat_spaces_start, line)
                if line_match is not None:
                    start_groups = line_match.groups()
                    if start_groups is not None:
                        prefix = start_groups[0]
                        start_tag = start_groups[1]
                        tag_content = start_groups[2]
                        if start_tag == 'date':
                            tag_content = clean_date(tag_content)
                        if start_tag == 'contourLevel':
                            tag_content = str(float(tag_content))
                        suffix = start_groups[5]
                        end_tag = start_groups[7]
                        if end_tag is not None and end_tag == start_tag:
                            empty_details_match = re.match(open_close_details_pat, line)
                            if empty_details_match is None:
                                cleaned_content = re.sub(pat_spaces_sub, pat_spaces_rep, tag_content).strip()
                                outf_hd.write('%s%s%s\n' % (prefix, cleaned_content, suffix))
                        else:
                            start_tag_found = True
                else:
                    line_match = re.match(pat, line)
                    if line_match is None:
                        if open_close_tags_on_next_line(open_pat, close_pat, line, index, lines):
                            index = index + 1
                        else:
                            line_match = re.match(open_close_pat, line)
                            if line_match is None:
                                outf_hd.write(line)

            else:
                line_match = re.match(pat_spaces_end, line)
                if line_match is not None:
                    end_groups = line_match.groups()
                    if end_groups is not None:
                        tag_content += ' ' + end_groups[0]
                        end_tag = end_groups[4]
                        suffix = end_groups[3]
                        if end_tag is not None and end_tag == start_tag:
                            cleaned_content = re.sub(pat_spaces_sub, pat_spaces_rep, tag_content).strip()
                            outf_hd.write('%s%s%s\n' % (prefix, cleaned_content, suffix))
                            start_tag_found = False
                    else:
                        tag_content += ' ' + line
                else:
                    tag_content += ' ' + line

    def clean_file(file_in, file_out):
        """
        """
        file_in_hd = open(file_in, 'r')
        outf_hd = open(file_out, 'w')
        """
        for line in tmpf_hd:
            m = re.match(pat, line)
            if m is None:
                outf_hd.write(line)
        """
        clean_spaces(file_in_hd, outf_hd)
        outf_hd.close()
        file_in_hd.close()

    emdb_files = glob.glob(file_path_template)
    num_errors = 0
    num_success = 0
    error_list = []
    i = 0
    for emdb_file in emdb_files:
        i = i + 1
        if i > 0:
            print i
            inf = os.path.basename(emdb_file)
            outf = os.path.join(out_dir, inf)
            tmpf = os.path.join(out_dir, inf + '.tmp')
            tmpf2 = os.path.join(out_dir, inf + '.tmp2')
            logging.info("Input file: %s, output file: %s", emdb_file, outf)
            command_list = list(command_list_base)
            command_list.append(tmpf)
            command_list.append(emdb_file)
            cmd_text = ' '.join(command_list)
            #logging.info('Executing: %s', cmd_text)
            #print 'commands %s' % command_list
            exit_code = subprocess.call(command_list)
            #print 'done'
            if exit_code != 0:
                num_errors += 1
                error_list.append(inf)
            else:
                clean_file(tmpf, tmpf2)
                clean_file(tmpf2, outf)
                os.remove(tmpf)
                os.remove(tmpf2)
                num_success += 1
    #logging.warning('%d files successfully processed!', num_success)
    if num_errors > 0:
        logging.warning('%d errors!', num_errors)
        logging.warning('List of entries that were not translated')
        for entry in error_list:
            logging.warning(entry)


def main():
    """
    Convert all EMDB XML 1.9 header files to canonical XML 1.9 files.
    This makes comparison with output from round-trip conversion (1.9 -> 3.0 -> 1.9) easier.
    """
    default_file_path_template = EMDBSettings.header_19_template
    default_out_dir = EMDBSettings.emdb_19_to_19_dir_out

    # Handle command line options
    usage = """
            python process_all_19_19.py [options]
            Convert all EMDB XML 1.9 header files to canonical XML 1.9 files.
            This makes comparison with output from round-trip conversion (1.9 -> 3.0 -> 1.9) easier.

            Examples:
            python process_all_19_19.py

            Typical run:
            python process_all_19_19.py -t '/data/emstaging/EMD-*/header/emd-*.xml' -o '/data/emdb19_to_19'
            /data/emstaging/EMD-*/header/emd-*.xml is the template used to glob all input 1.9 header files
            /data/emdb19_to_19 is the output directory with the canonical EMDB XML 1.9 files
            """
    version = "0.4"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-t", "--template", action="store", type="string", metavar="TEMPLATE", dest="filePathTemplate", default=default_file_path_template, help="Template used to glob all input 1.9 header files [default: %default]")
    parser.add_option("-o", "--out-dir", action="store", type="string", metavar="DIR", dest="outDir", default=default_out_dir, help="Directory for canonical EMDB 1.9 files [default: %default]")
    (options, args) = parser.parse_args()
    # print args
    process_all_19_19(options.filePathTemplate, options.outDir)


if __name__ == "__main__":
    main()
