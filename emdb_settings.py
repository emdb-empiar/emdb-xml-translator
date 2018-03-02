"""
emdb_settings.py

Global settings for project wrapped as a class to confine namespace

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

import logging
import os

__author__ = 'Ardan Patwardhan, Sanja Abbott'
__email__ = 'ardan@ebi.ac.uk, sanja@ebi.ac.uk'
__date__ = '2017-06-14'


class EMDBSettings(object):
    """
    This class holds some settings needed for translators
    """
    log_level = logging.INFO
    log_format = '%(asctime)s [%(levelname)s] %(name)s: %(message)s'

    # VERSION 3.0
    # input directory: contains EMDB XML 3.0 files
    emdb_30_dir_in = 'data/input/v3.0/'
    header_30_template = os.path.join(emdb_30_dir_in, '*.xml')
    # output directory: contains back translated EMDB XML 1.9 files
    emdb_30_to_19_dir_out = 'data/output/emdb30_to_19'
    # output directory: contains EMDB XML 3.0 translated back from translated v1.9
    emdb_19_to_30_dir_out = 'data/output/emdb19_to_30'
    # input directory: contains EMDB XML 1.9 files translated from v3.0
    emdb_19_from_30_dir_in = emdb_30_to_19_dir_out
    header_19_from_30_template = os.path.join(emdb_19_from_30_dir_in, '*.xml')
    # input directory: contains EMBD 1.9 from emprepare (EmDep)
    emdb_19_dir_in = 'data/input/v1.9/'
    header_19_template = os.path.join(emdb_19_dir_in, '*.xml')
    # output directory 1.9 to 1.9
    emdb_19_to_19_dir_out = 'data/output/emdb19_to_19'
    header_19_19_template = os.path.join(emdb_19_to_19_dir_out, '*.xml')
    emdb_19_to_30relax_dir_out = 'data/output/emdb19_to_30relax'
    emdb_19_to_19_via_30relax_dir_out = 'data/output/emdb30relax_to_19'
    # directory to outbut diffs between back-converted 1.9 files and the canonical ones
    diff_19s_dir = 'data/output/emdb_19s_diff'
    # schemas:
    schema30 = "emdb30.xsd"
    schema19 = "emdb19.xsd"
    schema30relaxed = "emdb30_relaxed.xsd"
