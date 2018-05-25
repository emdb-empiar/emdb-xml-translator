#!/usr/bin/env python
"""
emdb_xml_translate.py

Convert EMDB XML files from one schema version to another.

Copyright [2014-2016] EMBL - European Bioinformatics Institute
Licensed under the Apache License, Version 3.0 (the
"License"); you may not use this file except in
compliance with the License. You may obtain a copy of
the License at
http://www.apache.org/licenses/LICENSE-3.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on
an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
"""
import sys
import logging
import traceback
import string
import re

from optparse import OptionParser
from lxml import etree
from dateutil import parser as dtp
from emdb_settings import EMDBSettings

import emdb_30
import emdb_30relaxed
import emdb_19

__author__ = 'Ardan Patwardhan, Sanja Abbott'
__email__ = 'ardan@ebi.ac.uk, sanja@ebi.ac.uk'
__date__ = '2017-06-14'


class EMDBXMLTranslator(object):
    """
    Class for translating EMDB files 3.0 <-> 1.9
    """

    class Constants(object):
        """
        There are many constants in use for the translation. They have been collected here for  ease of use.
        """
        # Global constants
        LATEST_XML_VERSION = '3.0.0.0'
        # Used to identify sample supramolecule
        EM_SAMPLE_ID = 1000
        EM_DATE_FORMAT = '%d-%b-%Y'
        EM_UNIDENTIFIED_TAXID = 32644
        EMDB_PAT = re.compile(r'(?i)(EMD-){0,1}(\d{4,})')
        EMDB_PREFIX = 'EMD-'
        EMDB_DUMMY_CODE = 'EMD-0000'
        PDB_CHAIN_PAT = re.compile(r'(\d[\dA-Za-z]{3})([-_:; ]?)([A-Za-z0-9]+)')
        EUL_ANG_START_TAG = '{eulerAngleDetails}'
        EUL_ANG_END_TAG = '{/eulerAngleDetails}'
        EUL_ANG_PAT = re.compile(r'(.*)%s(.*)%s(.*)' % (EUL_ANG_START_TAG, EUL_ANG_END_TAG))
        HEL_TAG = '{helical/}'
        SP_TAG = '{singleParticle/}'
        HEL_SP_PAT = re.compile(r'(.*){(helical|singleParticle)/}(.*)')

        # EM methods
        EMM_EC = 'electronCrystallography'
        EMM_HEL = 'helical'
        EMM_SP = 'singleParticle'
        EMM_STOM = 'subtomogramAveraging'
        EMM_TOM = 'tomography'

        # Units
        U_ANG = u'\u212B'
        U_DEG = 'deg'
        U_DEGF = 'degrees'
        U_FIB_DOSE_RATE = 'ions/nm^2/s'
        U_KDA_NM = 'kDa/nm'
        U_KEL = 'K'
        U_KELF = 'Kelvin'
        U_KVOLT = 'kV'
        U_MCRN = 'microns'
        U_MDA = 'MDa'
        U_MG_ML = 'mg/mL'
        U_MM = 'mm'
        U_NM = 'nm'
        U_PAMP = 'pA'
        U_SEC = 's'
        U_MICROM = u'\u00B5' + 'm'
        U_EOVERANGSQR = 'e/' + U_ANG + '^2'
        U_PERCENTAGE = 'percentage'
        U_EL_A2 = 'e/A**2'

        # Status
        STS_REL = 'REL'
        STS_HPUB = 'HPUB'
        STS_HOLD = 'HOLD'
        STS_HOLD1 = 'HOLD1'
        STS_OBS = 'OBS'

        # Extension types
        EXT_BASE_MICROSCOPY_TYPE = 'base_microscopy_type'
        EXT_TOMOGRAPHY_MICROSCOPY_TYPE = 'tomography_microscopy_type'
        EXT_CRYST_MICROSCOPY_TYPE = 'crystallography_microscopy_type'
        EXT_BASE_PREPARATION_TYPE = 'base_preparation_type'
        EXT_TOMOGRAPHY_PREPARATION_TYPE = 'tomography_preparation_type'
        EXT_CRYST_PREPARATION_TYPE = 'crystallography_preparation_type'

        # Dictionaries for translations
        DATA_TYPE_DICT_19_TO_30 = {"Envelope stored as signed bytes": "IMAGE STORED AS SIGNED BYTE",
                                   "Image stored as Integer*2": "IMAGE STORED AS SIGNED INTEGER (2 BYTES)",
                                   "Image stored as Reals": "IMAGE STORED AS FLOATING POINT NUMBER (4 BYTES)"}

        # Cleaning up dictionaries for translation from 20 to 19
        PROC_SITE_30_TO_19 = {'pdbe': 'PDBe', 'rcsb': 'RCSB', 'pdbj': 'PDBj'}

        FITTING_30_to_19 = {'AB INITIO MODEL': 'flexible',
                            'BACKBONE TRACE': 'flexible',
                            'FLEXIBLE FIT': 'flexible',
                            'OTHER': 'flexible',
                            'RIGID BODY FIT': 'rigid body'}

        FITTING_19_to_30 = {'flexible': 'FLEXIBLE FIT',
                            'rigid body': 'RIGID BODY FIT'}

        SPECIMEN_HOLDER_30_to_19 = {'FEI TITAN KRIOS AUTOGRID HOLDER': 'FEI TITAN KRIOS AUTOGRID HOLDER',
                                    'GATAN 626 SINGLE TILT LIQUID NITROGEN CRYO TRANSFER HOLDER': 'GATAN LIQUID NITROGEN',
                                    'GATAN 910 MULTI-SPECIMEN SINGLE TILT CRYO TRANSFER HOLDER': 'GATAN LIQUID NITROGEN',
                                    'GATAN 914 HIGH TILT LIQUID NITROGEN CRYO TRANSFER TOMOGRAPHY HOLDER': 'GATAN LIQUID NITROGEN',
                                    'GATAN 915 DOUBLE TILT LIQUID NITROGEN CRYO TRANSFER HOLDER': 'GATAN LIQUID NITROGEN',
                                    'GATAN CHDT 3504 DOUBLE TILT HIGH RESOLUTION NITROGEN COOLING HOLDER': 'GATAN LIQUID NITROGEN',
                                    'GATAN CT3500 SINGLE TILT LIQUID NITROGEN CRYO TRANSFER HOLDER': 'GATAN LIQUID NITROGEN',
                                    'GATAN CT3500TR SINGLE TILT ROTATION LIQUID NITROGEN CRYO TRANSFER HOLDER': 'GATAN LIQUID NITROGEN',
                                    'GATAN HC 3500 SINGLE TILT HEATING/NITROGEN COOLING HOLDER': 'GATAN LIQUID NITROGEN',
                                    'GATAN HCHDT 3010 DOUBLE TILT HIGH RESOLUTION HELIUM COOLING HOLDER': 'GATAN HELIUM',
                                    'GATAN HCHST 3008 SINGLE TILT HIGH RESOLUTION HELIUM COOLING HOLDER': 'GATAN HELIUM',
                                    'GATAN HELIUM': 'GATAN HELIUM',
                                    'GATAN LIQUID NITROGEN': 'GATAN LIQUID NITROGEN',
                                    'GATAN UHRST 3500 SINGLE TILT ULTRA HIGH RESOLUTION NITROGEN COOLING HOLDER ': 'GATAN LIQUID NITROGEN',
                                    'GATAN ULTDT ULTRA LOW TEMPERATURE DOUBLE TILT HELIUM COOLING HOLDER': 'GATAN HELIUM',
                                    'GATAN ULTST ULTRA LOW TEMPERATURE SINGLE TILT HELIUM COOLING HOLDER': 'GATAN HELIUM',
                                    'HOME BUILD': 'HOME BUILD',
                                    'JEM3400FSC CRYOHOLDER': 'JEOL',
                                    'JEOL ': 'JEOL',
                                    'JEOL 3200FSC CRYOHOLDER': 'JEOL 3200FSC CRYOHOLDER',
                                    'OTHER': 'OTHER',
                                    'PHILIPS ROTATION HOLDER': 'PHILIPS ROTATION HOLDER',
                                    'SIDE ENTRY, EUCENTRIC': 'SIDE ENTRY, EUCENTRIC',
                                    'JEOL3200FSC CRYOHOLDER': 'JEOL 3200FSC CRYOHOLDER'}

        SPECIMEN_STATE_30_to_19 = {'particle': 'particle',
                                   'filament': 'filament',
                                   'twodarray': 'twoDArray',
                                   'threedarray': 'threeDArray',
                                   'helicalarray': 'helicalArray',
                                   'tissue': 'tissue',
                                   'cell': 'cell'}

    def __init__(self):
        # 0 = min, 3 = max
        self.warning_level = 1
        self.validate_xml_out = False
        self.relaxed = False
        self.roundtrip = False
        logging.basicConfig(level=EMDBSettings.log_level, format=EMDBSettings.log_format)

    def set_roundtrip(self, is_roundtrip):
        """
        Set flag to determine v3.0 schema
        """
        self.roundtrip = is_roundtrip

    def set_v30_schema(self, is_relaxed):
        """
        Set flag to determine v3.0 schema
        """
        self.relaxed = is_relaxed

    def set_validate(self, to_validate):
        """
        Set validation flag

        Parameters:
        @param to_validate: True or False (default)
        """
        self.validate_xml_out = to_validate

    def set_warning_level(self, level):
        """
        Set the level of logging warnings. 0 = no warnings, 3 = max warnings, 1 = default

        Parameters
        @param level: warning level 0 -> 3
        """
        if level <= 0:
            self.warning_level = 0
        elif level >= 3:
            self.warning_level = 3
        else:
            self.warning_level = level

    def warn(self, level, msg):
        """
        Log a warning message but take into account the warning_level

        Parameters:
        @param level: only messages with level >= warning_level are displayed
        @param msg: warning message
        """
        if level <= self.warning_level:
            logging.warning(msg)

    def check_set(self, get_value, set_value, transform=None):
        """
        Call set_value only if get_value does not return None

        Parameters:
        @param get_value: getter function that must return value
        @param set_value: setter function
        @param transform: Apply transform(x) before calling setter function
        """
        value = get_value()
        if value is not None:
            if transform is not None:
                try:
                    original_value = value
                    value = transform(original_value)
                except Exception as exp:
                    self.warn(3, "function check_set: Transform function did not work: %s(%s).  Error: %s" % (transform, original_value, exp))
                    self.warn(3, traceback.format_exc())
                    return
            try:
                set_value(value)
            except Exception as exp:
                self.warn(3, "function check_set: Setter function did not work: %s(%s). Error: %s" % (set_value, value, exp))
                self.warn(3, traceback.format_exc())

    def set_value_and_units(self, getter, setter, constructor, units=None, transform=None):
        """
        There are several elements that take a value and have an units attribute.
        This function makes it easier to copy over these elements

        Parameters:
        @param getter: Getter function to get (value,units)
        @param setter: Setter function to set (value,units)
        @param constructor: Constructor for object that takes (units, valueOf_) as params
        @param units: If this is not None then set units based on this otherwise transfer value
        @param transform: Apply transform(x) before calling setter function
        """
        x_value = getter()
        if x_value is not None:
            if units is None:
                units_value = x_value.get_units()
            else:
                units_value = units

            y_value = constructor(valueOf_=x_value.get_valueOf_(), units=units_value)
            if transform is not None:
                try:
                    z_value = transform(y_value)
                except Exception as exp:
                    self.warn(3, "function set_value_and_units: Transform function did not work: %s(%s). Error: %s" % (transform, y_value, exp))
                    self.warn(3, traceback.format_exc())
            else:
                z_value = y_value
            try:
                setter(z_value)
            except Exception as exp:
                self.warn(3, "function set_value_and_units: Setter function did not work: %s(%s). Error: %s" % (setter, z_value, exp))
                self.warn(3, traceback.format_exc())

    def format_emdb_code(self, code_in, number_only=False):
        """
        Format code so that it is either in the form EMD-XXXX or just XXXX

        Parameters:
        @param code_in: string representing EMDB code in some legal format, e.g. emd-xxxx, EMD-xxxx or xxxx
        @param number_only: Boolean if true - return only xxxx otherwise return EMD-xxxx
        @return: EMDB accession code in specified format
        """
        access_code = self.Constants.EMDB_DUMMY_CODE
        mtch = re.match(self.Constants.EMDB_PAT, code_in)
        if mtch is None:
            self.warn(1, 'EMDB accession code: %s does not match any standards. Using dummy code: %s' % (code_in, access_code))
            mtch = re.match(self.Constants.EMDB_PAT, access_code)
        match_groups = mtch.groups()
        if number_only:
            return match_groups[1]
        else:
            return '%s%s' % (self.Constants.EMDB_PREFIX, match_groups[1])

    def translate_1_9_to_3_0(self, input_file, output_file):
        """
        Convert input file from 1.9 to 3.0 schema

        Parameters:
        @param input_file: Name of input file
        @param output_file: Name of output file
        """
        const = self.Constants
        emdb30 = emdb_30
        if self.relaxed:
            emdb30 = emdb_30relaxed

        def make_software_list(soft_in):
            """
            Takes a string representing software and create a software list (3.0 construct).
            Convenience function for translating from 1.9 to 3.0

            Parameters:
            @param soft_in: software represented as string
            @return: software list as software_list_type (3.0)
            """
            if soft_in is not None:
                soft_list = emdb30.software_list_type()
                soft = emdb30.software_type()
                soft.set_name(soft_in)
                soft_list.add_software(soft)
                return soft_list
            else:
                return None

        def add_external_references(ref_in, ref_out):
            """
            Copy over reference list for journals or non-journals

            Parameters:
            @param ref_in: Input citation with reference list
            @param ref_out: Output citation to which reference list is added.
            """
            ext_ref_in = ref_in.get_externalReference()
            for ref in ext_ref_in:
                ext_ref_out = emdb30.external_referencesType()
                # XSD: <xs:attribute name="type" use="required">
                ext_ref_out.set_type(ref.get_type().upper())
                ext_ref_out.set_valueOf_(ref.get_valueOf_())
                ref_out.add_external_references(ext_ref_out)

        def copy_authors(get_authors, add_author, simple=False):
            """
            Copy authors from 1.9 -> 3.0 while reformatting them

            Parameters
            @param get_authors: getter function for getting authors from jrnl/nonjrnl object of 1.9
            @param add_author: adding (setter) function for adding an author to the list of jrnl/nonjrnl object authors
            @param simple: boolean - True means that the authors in 3.0 are simple strings, otherwise they are journal authors
            """
            authors_in = get_authors()
            if authors_in is not None:
                auth_in = authors_in.replace(', ', ',').split(',')
                index = 1
                for auth_str_in in auth_in:
                    if auth_str_in != '':
                        if simple is False:
                            author = emdb30.author_order_type()
                            author.set_valueOf_(auth_str_in)
                            author.set_order(index)
                            if author.hasContent_():
                                add_author(author)
                                index += 1
                        else:
                            add_author(auth_str_in)

        def copy_citation(cit_in, cit_out):
            """
            Copy over citation from 1.9 to 3.0

            Parameters:
            @param cit_in: Input citation in 1.9 schema
            @param cit_out: Output citation in 3.0 schema
            """
            jrnl_in = cit_in.get_journalArticle()
            if jrnl_in:
                # XSD: <xs:element name="journal_citation" substitutionGroup="citation_type"> has 13 elements and 1 attribute
                jrnl = emdb30.journal_citation()
                jrnl.original_tagname_ = 'journal_citation'
                # attribute 1 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:attribute name="published" type="xs:boolean" use="required"/>
                jrnl.set_published(cit_in.get_published())
                # element 1 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="author" type="author_order_type" maxOccurs="unbounded"/>
                copy_authors(jrnl_in.get_authors, jrnl.add_author)
                # element 2 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="title" type="xs:token"/>
                jrnl.set_title(jrnl_in.get_articleTitle())
                # element 3 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="journal" type="xs:token" minOccurs="0"/>
                jrnl.set_journal(jrnl_in.get_journal())
                # element 4 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="journal_abbreviation" type="xs:token">
                #jrnl.set_journal_abbreviation('')
                # element 5 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="country" type="xs:token" minOccurs="0"/>
                # element 6 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="issue" type="xs:positiveInteger" minOccurs="0"/>
                # element 7 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="volume" type="xs:string" nillable="true" minOccurs="0"/>
                # This is a fix because of bad data - emd-1648.xml has an empty volume tag!
                vol = jrnl_in.get_volume()
                if vol is not None and len(vol) > 0:
                    jrnl.set_volume(vol)
                # element 8 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="first_page" type="page_type" nillable="false" minOccurs="0"/>
                first_page = jrnl_in.get_firstPage()
                if first_page is not None and len(first_page) > 0:
                    jrnl.set_first_page(first_page)
                # element 9 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="last_page" type="page_type" minOccurs="0"/>
                last_page = jrnl_in.get_lastPage()
                if last_page is not None and len(last_page) > 0:
                    jrnl.set_last_page(last_page)
                # element 10 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="year" minOccurs="0">
                year = jrnl_in.get_year()
                if year is not None and len(year) > 0:
                    jrnl.set_year(year)
                # element 11 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="language" type="xs:language" minOccurs="0"/>
                # element 12 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                add_external_references(jrnl_in, jrnl)
                # element 13 - <xs:element name="journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>

                cit_out.set_citation_type(jrnl)
            else:
                non_jrnl_in = cit_in.get_nonJournalArticle()
                # XSD: <xs:element name="non_journal_citation" substitutionGroup="citation_type"> has 14 elements and 1 attribute
                non_jrnl = emdb30.non_journal_citation()
                non_jrnl.original_tagname_ = 'non_journal_citation'
                # attribute 1 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:attribute name="published" type="xs:boolean" use="required"/>
                non_jrnl.set_published(cit_in.get_published())
                # element 1 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="author" type="author_order_type" maxOccurs="unbounded"/>
                copy_authors(non_jrnl_in.get_authors, non_jrnl.add_author)
                # element 2 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="editor" type="author_order_type" minOccurs="0" maxOccurs="unbounded"/>
                copy_authors(non_jrnl_in.get_editor, non_jrnl.add_editor)
                # element 3 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="title" type="xs:token"/>
                self.check_set(non_jrnl_in.get_book, non_jrnl.set_title)
                # element 4 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="thesis_title" type="xs:token" minOccurs="0"/>
                self.check_set(non_jrnl_in.get_thesisTitle, non_jrnl.set_thesis_title)
                # element 5 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="chapter_title" type="xs:token" minOccurs="0"/>
                self.check_set(non_jrnl_in.get_chapterTitle, non_jrnl.set_chapter_title)
                # element 6 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="volume" type="xs:string" minOccurs="0"/>
                self.check_set(non_jrnl_in.get_volume, non_jrnl.set_volume)
                # element 7 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="publisher" type="xs:token" minOccurs="0"/>
                self.check_set(non_jrnl_in.get_publisher, non_jrnl.set_publisher)
                # element 8 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="publisher_location" type="xs:token" minOccurs="0"/>
                self.check_set(non_jrnl_in.get_publisherLocation, non_jrnl.set_publisher_location)
                # element 9 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="first_page" type="page_type" minOccurs="0"/>
                self.check_set(non_jrnl_in.get_firstPage, non_jrnl.set_first_page)
                # element 10 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="last_page" type="page_type" minOccurs="0"/>
                self.check_set(non_jrnl_in.get_lastPage, non_jrnl.set_last_page)
                # element 11 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="year">
                self.check_set(non_jrnl_in.get_year, non_jrnl.set_year)
                # element 12 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="language" type="xs:language" minOccurs="0"/>
                # element 13 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">
                add_external_references(non_jrnl_in, non_jrnl)
                # element 14 - <xs:element name="non_journal_citation" substitutionGroup="citation_type">
                # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>

                cit_out.set_citation_type(non_jrnl)

        def copy_recombinant_source(eng_in, set_source_func):
            """
            Copy engineered source from 1.9 to recombinant source in 3.0

            Parameters:
            @param eng_in: Engineered source object (input)
            @param set_source_func: Method belonging to molecule/supramolecule object that sets its expression system, e_value.g. set_source_func = mol.set_recombinant_expression for protein
            """
            # XSD: <xs:complexType name="recombinant_source_type"> has 5 elements and 1 attribute
            if eng_in is not None:
                recs = emdb30.recombinant_source_type()
                # attribute 1 - <xs:complexType name="recombinant_source_type">
                # XSD: <xs:attribute name="database" use="required">
                recs.set_database('NCBI')
                # element 1 - <xs:complexType name="recombinant_source_type">
                # XSD: <xs:element name="recombinant_organism" type="organism_type"/>
                exp_sys_in = eng_in.get_expSystem()
                if exp_sys_in is not None:
                    # XSD: <xs:complexType name="organism_type">; extension base="xs:token">
                    org = emdb30.organism_type()
                    self.check_set(exp_sys_in.get_valueOf_, org.set_valueOf_)
                    # XSD: <xs:attribute name="ncbi" type="xs:positiveInteger"/>
                    self.check_set(exp_sys_in.get_ncbiTaxId, org.set_ncbi)
                    recs.set_recombinant_organism(org)
                # element 2 - <xs:complexType name="recombinant_source_type">
                # XSD: <xs:element name="recombinant_strain" type="xs:token" minOccurs="0"/>
                self.check_set(eng_in.get_expSystemStrain, recs.set_recombinant_strain)
                # element 3 - <xs:complexType name="recombinant_source_type">
                # XSD: <xs:element name="recombinant_cell" type="xs:token" minOccurs="0"/>
                self.check_set(eng_in.get_expSystemCell, recs.set_recombinant_cell)
                # element 4 - <xs:complexType name="recombinant_source_type">
                # XSD: <xs:element name="recombinant_plasmid" type="xs:token" minOccurs="0"/>
                self.check_set(eng_in.get_vector, recs.set_recombinant_plasmid)
                # element 5 - <xs:complexType name="recombinant_source_type">
                # XSD: <xs:element name="recombinant_synonym_organism" type="xs:token" minOccurs="0"/>
                set_source_func(recs)

        def set_sci_name(c_in, c_out):
            """
            Create v3.0 supramolecule/macromolecule name object

            Parameters:
            @param c_in: v1.9 component with getter functions for name
            @param c_out: v3.0 component with setter functions for name
            """
            # XSD: <xs:complexType name="sci_name_type">, extension of sx:token, has 1 attribute
            name = emdb30.sci_name_type()
            sci_name = c_in.get_sciName()
            if sci_name is not None:
                name.set_valueOf_(sci_name)
            # attribute 1 - <xs:complexType name="sci_name_type">
            # XSD: <xs:attribute name="synonym" type="xs:token"/>
            syn_name = c_in.get_synName()
            if syn_name is not None:
                name.set_synonym(syn_name)
            c_out.set_name(name)

        def set_mol_weight(setter_func, wt_theo_in=None, wt_exp_in=None, wt_meth_in=None):
            """
            Set molecular weight if provided

            Parameters:
            @param setter_func: Function to set molecular weight: fig_in(mol_wt)
            @param wt_theo_in: Theoretical molecular weight
            @param wt_exp_in: Experimental molecular weight
            @param wt_meth_in: Method used for calculating experimental weight
            """
            if wt_theo_in is not None or wt_exp_in is not None or wt_meth_in is not None:
                mol_wt = emdb30.molecular_weight_type()
                setter_func(mol_wt)
                if wt_exp_in is not None:
                    mol_wt.set_experimental(emdb30.experimentalType(valueOf_=wt_exp_in.get_valueOf_(), units=wt_exp_in.get_units()))
                if wt_theo_in is not None:
                    mol_wt.set_theoretical(emdb30.theoreticalType(valueOf_=wt_theo_in.get_valueOf_(), units=wt_theo_in.get_units()))
                if wt_meth_in is not None:
                    mol_wt.set_method(wt_meth_in)

        def add_reference(x_ref_in, adder_func):
            """
            """
            ref_uni_prots = x_ref_in.get_refUniProt()
            for ref_uni_prot in ref_uni_prots:
                if ref_uni_prot is not None:
                    adder_func(emdb30.external_referencesType(valueOf_=ref_uni_prot.upper(), type_='UNIPROTKB'))
            ref_gos = x_ref_in.get_refGo()
            for ref_go in ref_gos:
                if ref_go is not None:
                    adder_func(emdb30.external_referencesType(valueOf_=ref_go, type_='GO'))
            ref_inters = x_ref_in.get_refInterpro()
            for ref_inter in ref_inters:
                if ref_inter is not None:
                    adder_func(emdb30.external_referencesType(valueOf_=ref_inter, type_='INTERPRO'))

        def add_mol_references(adder_func, x_refs_in):
            """
            Add external references to molecule sequence

            Parameters:
            @param adder_func: Adder function to add external references, e_value.g.seq.add_external_references(map_in_file)
            @param x_refs_in: v1.9 externalReferences object
            """
            if x_refs_in is not None:
                if isinstance(x_refs_in, list):
                    for x_ref_in in x_refs_in:
                        add_reference(x_ref_in, adder_func)
                else:
                    add_reference(x_refs_in, adder_func)

        def copy_map_19_to_30(map_in, map_out, is_map=False):
            """
            Copy map from 1.9 to 3.0

            Parameters:
            @param map_in: input 1.9 map
            @param map_out: output 3.0 map
            @param spec_prep_in: specimen preparation from 1.9
            """
            # Set file and related attributes
            map_in_file = map_in.get_file()
            # map_out is of <xs:complexType name="map_type"> has 14 elements and 2 attributes
            # attribute 1 - <xs:complexType name="map_type">
            # XSD: <xs:attribute name="format" fixed="CCP4" use="required"/>
            map_out.set_format(map_in_file.get_format())
            # attribute 2 - <xs:complexType name="map_type">
            # XSD: <xs:attribute name="size_kbytes" type="xs:positiveInteger" use="required"/>
            map_out.set_size_kbytes(map_in_file.get_sizeKb())
            # element 1 - <xs:complexType name="map_type">
            # XSD: <xs:element name="file">
            map_filename = map_in_file.get_valueOf_()
            if map_filename != '':
                map_out.set_file(map_filename)
            else:
                # this is an annotation problem
                self.warn(1, "No file name given for a map!")
                map_out.set_file('emd_0000.map.gz')
            # element 2 - <xs:complexType name="map_type">
            # XSD: <xs:element name="symmetry" type="applied_symmetry_type" minOccurs="0"> has 1 choice of 3 elements
            sym = emdb30.applied_symmetry_type()
            # choice 1 - <xs:element name="symmetry" type="applied_symmetry_type" minOccurs="0">
            # XSD: <xs:element name="space_group" type="xs:token"/>
            grp_num = map_in.get_spaceGroupNumber()
            if grp_num is not None:
                sym.set_space_group(grp_num)
            # choice 2 - <xs:element name="symmetry" type="applied_symmetry_type" minOccurs="0">
            # XSD: <xs:element name="point_group">
            # choice 3 - <xs:element name="symmetry" type="applied_symmetry_type" minOccurs="0">
            # XSD: <xs:element name="helical_parameters" type="helical_parameters_type">
#             if spec_prep_in is not None:
#                 # XSD: <xs:complexType name="helical_parameters_type"> has 3 elements
#                 hel = emdb30.helical_parameters_type()
#                 hel_in = spec_prep_in.get_helicalParameters()
#                 if hel_in is not None:
#                     # sym = emdb30.applied_symmetry_type()
#                     # element 1 - <xs:complexType name="helical_parameters_type">
#                     # XSD: <xs:element name="delta_z">
#                     self.set_value_and_units(hel_in.get_deltaZ, hel.set_delta_z, emdb30.delta_zType, units=const.U_ANG)
#                     # element 2 - <xs:complexType name="helical_parameters_type">
#                     # XSD: <xs:element name="delta_phi">
#                     self.set_value_and_units(hel_in.get_deltaPhi, hel.set_delta_phi, emdb30.delta_phiType, units=const.U_DEG)
#                     # element 3 - <xs:complexType name="helical_parameters_type">
#                     # XSD: <xs:element name="axial_symmetry">
#                     self.check_set(hel_in.get_axialSymmetry, hel.set_axial_symmetry)
#                     # not in the schema anymore
#                     # self.check_set(hel_in.get_hand, hel.set_hand)
#                     if hel.hasContent_():
#                         sym.set_helical_parameters(hel)

            map_out.set_symmetry(sym)
            # element 3 - <xs:complexType name="map_type">
            # XSD: <xs:element name="data_type" type="map_data_type"/>
            data_type_in = map_in.get_dataType()
            data_type_out = const.DATA_TYPE_DICT_19_TO_30.get(data_type_in, None)
            map_out.set_data_type(data_type_out)
            # element 4 - <xs:complexType name="map_type">
            # XSD: <xs:element name="dimensions" type="integer_vector_map_type"/>
            dim_in = map_in.get_dimensions()
            # XSD: <xs:complexType name="integer_vector_map_type"> has 3 elements
            # element 1 - <xs:element name="col" type="xs:positiveInteger"/>
            # element 2 - <xs:element name="row" type="xs:positiveInteger"/>
            # element 3 - <xs:element name="sec" type="xs:positiveInteger"/>
            if dim_in is not None:
                dim_out = emdb30.integer_vector_map_type(row=dim_in.get_numRows(), col=dim_in.get_numColumns(), sec=dim_in.get_numSections())
                map_out.set_dimensions(dim_out)
            # element 5 - <xs:complexType name="map_type">
            # XSD: <xs:element name="origin"> has 3 elements
            orig_in = map_in.get_origin()
            # element 1 - <xs:element name="col" type="xs:integer"/>
            # element 2 - <xs:element name="row" type="xs:integer"/>
            # element 3 - <xs:element name="sec" type="xs:integer"/>
            orig_out = emdb30.originType(col=orig_in.get_originCol(), row=orig_in.get_originRow(), sec=orig_in.get_originSec())
            if orig_in is not None:
                map_out.set_origin(orig_out)
            # element 6 - <xs:complexType name="map_type">
            # XSD: <xs:element name="spacing"> has 3 elements
            spc_in = map_in.get_spacing()
            # element 1 - <xs:element name="x" type="xs:positiveInteger"/>
            # element 2 - <xs:element name="y" type="xs:nonNegativeInteger"/>
            # element 3 - <xs:element name="z" type="xs:nonNegativeInteger"/>
            if spc_in is not None:
                spc_out = emdb30.spacingType(spc_in.get_spacingRow(), spc_in.get_spacingCol(), spc_in.get_spacingSec())
                map_out.set_spacing(spc_out)
            # element 7 - <xs:complexType name="map_type">
            # XSD: <xs:element name="cell"> has 6 elements
            cell_in = map_in.get_cell()
            # element 1 - <xs:element name="a" type="cell_type"/>
            # element 2 - <xs:element name="b" type="cell_type"/>
            # element 3 - <xs:element name="c" type="cell_type"/>
            # element 4 - <xs:element name="alpha" type="cell_angle_type"/>
            # element 5 - <xs:element name="beta" type="cell_angle_type"/>
            # element 6 - <xs:element name="gamma" type="cell_angle_type"/>
            if cell_in is not None:
                cell_out = emdb30.cellType(a=emdb30.cell_type(valueOf_=cell_in.get_cellA().get_valueOf_(), units=const.U_ANG),
                                           b=emdb30.cell_type(valueOf_=cell_in.get_cellB().get_valueOf_(), units=const.U_ANG),
                                           c=emdb30.cell_type(valueOf_=cell_in.get_cellC().get_valueOf_(), units=const.U_ANG),
                                           alpha=emdb30.cell_angle_type(valueOf_=cell_in.get_cellAlpha().get_valueOf_(), units=const.U_DEG),
                                           beta=emdb30.cell_angle_type(valueOf_=cell_in.get_cellBeta().get_valueOf_(), units=const.U_DEG),
                                           gamma=emdb30.cell_angle_type(valueOf_=cell_in.get_cellGamma().get_valueOf_(), units=const.U_DEG))
                map_out.set_cell(cell_out)
            # element 8 - <xs:complexType name="map_type">
            # XSD: <xs:element name="axis_order"> has 3 elements
            axis_in = map_in.get_axisOrder()
            # element 1 - <xs:element name="fast">
            # element 2 - <xs:element name="medium">
            # element 3 - <xs:element name="slow">
            if axis_in is not None:
                axis_out = emdb30.axis_orderType(fast=axis_in.get_axisOrderFast(), medium=axis_in.get_axisOrderMedium(), slow=axis_in.get_axisOrderSlow())
                map_out.set_axis_order(axis_out)
            # element 9 - <xs:complexType name="map_type">
            # XSD: <xs:element name="statistics" type="map_statistics_type" minOccurs="0"/>
            map_out.set_statistics(map_in.get_statistics())
            # element 10 - <xs:complexType name="map_type">
            # XSD: <xs:element name="pixel_spacing"> has 3 elements
            pix_in = map_in.get_pixelSpacing()
            # element 1 - <xs:element name="x" type="pixel_spacing_type"/>
            # element 2 - <xs:element name="y" type="pixel_spacing_type"/>
            # element 3 - <xs:element name="z" type="pixel_spacing_type"/>
            if pix_in is not None:
                pix_out = emdb30.pixel_spacingType(emdb30.pixel_spacing_type(valueOf_=pix_in.get_pixelX().get_valueOf_(), units=const.U_ANG),
                                                   y=emdb30.pixel_spacing_type(valueOf_=pix_in.get_pixelY().get_valueOf_(), units=const.U_ANG),
                                                   z=emdb30.pixel_spacing_type(valueOf_=pix_in.get_pixelZ().get_valueOf_(), units=const.U_ANG))
                map_out.set_pixel_spacing(pix_out)
            # element 11 - <xs:complexType name="map_type">
            # XSD: <xs:element name="contour_list" minOccurs="0"> has 1 element
            # masks do not have contour level
            add_contour_level_details = ''
            if is_map:
                cntr_in = map_in.get_contourLevel()
                if cntr_in is not None:
                    cntr_list = emdb30.contour_listType()
                    # element 1 - <xs:element name="contour_list" minOccurs="0">
                    # XSD: <xs:element name="contour" maxOccurs="unbounded"> has 2 elements and 1 attribute
                    cntr = emdb30.contourType()
                    # attribute 1 - <xs:element name="contour" maxOccurs="unbounded">
                    # XSD: <xs:attribute name="primary" type="xs:boolean" use="required"/>
                    cntr.set_primary(True)
                    # element 1 - <xs:element name="contour" maxOccurs="unbounded">
                    # XSD: <xs:element name="level" type="xs:float">
                    level = cntr_in.get_valueOf_()
                    if level is not None:
                        if level.find('.') == -1:
                            # this is a whole number; note in details
                            add_contour_level_details = '{level is a whole number}'
                        cntr.set_level(float(level))
                    # element 2 - <xs:element name="contour" maxOccurs="unbounded">
                    # XSD: <xs:element name="source" minOccurs="0">
                    source_in = cntr_in.get_source()
                    if source_in is not None:
                        source_upper = source_in.upper()
                        cntr.set_source(source_upper)
                    cntr_list.add_contour(cntr)
                    map_out.set_contour_list(cntr_list)

            # element 12 - <xs:complexType name="map_type">
            # XSD: <xs:element name="label" type="xs:token" minOccurs="0"/>
            #self.check_set(map_in.get_label, map_out.set_label)
            # element 13 - <xs:complexType name="map_type">
            # XSD: <xs:element name="annotation_details" type="xs:string" minOccurs="0"/>
            self.check_set(map_in.get_annotationDetails, map_out.set_annotation_details)
            # element 14 - <xs:complexType name="map_type">
            # XSD: <xs:element name="details" type="xs:string" minOccurs="0">
            map_details = map_in.get_details()
            all_details = add_contour_level_details
            if map_details is not None:
                all_details = map_details + add_contour_level_details
            if all_details != '':
                map_out.set_details(all_details)

        def set_base_preparation(prep, spec_prep_in, vitr_in):
            """
            Method that sets the elements common to all preparation types

            Parameters:
            @params: prep - the object of one of 5 preparation methods
            @param: spec_prep_in: specimen preparation from 1.9
            @param: vitr_in: vitrification object in experiment
            """
            # XSD: <xs:complexType name="base_preparation_type"> has 8 elements
            # element 1 - <xs:complexType name="base_preparation_type">
            # XSD: <xs:element name="concentration" minOccurs="0"> has 1 attribute
            if spec_prep_in is not None:
                conc_in = spec_prep_in.get_specimenConc()
                if conc_in is not None:
                    conc = emdb30.concentrationType()
                    conc.set_valueOf_(conc_in.get_valueOf_())
                    # attribute 1
                    # XSD: <xs:attribute name="units" use="required">
                    conc.set_units('mg/mL')
                    prep.set_concentration(conc)
                # element 2 - <xs:complexType name="base_preparation_type">
                # XSD: <xs:element name="buffer" type="buffer_type" minOccurs="0">
                buf_in = spec_prep_in.get_buffer()
                if buf_in is not None:
                    # XSD: <xs:complexType name="buffer_type"> has 3 elements
                    buf = emdb30.buffer_type()
                    # element 1 - <xs:complexType name="buffer_type">
                    # XSD: <xs:element name="ph">
                    self.check_set(buf_in.get_ph, buf.set_ph)
                    # element 2 - <xs:complexType name="buffer_type">
                    # XSD: <xs:element name="component" maxOccurs="unbounded" type="buffer_component_type" minOccurs="0">
                    # element 3 - <xs:complexType name="buffer_type">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                    self.check_set(buf_in.get_details, buf.set_details)
                    prep.set_buffer(buf)
                # element 3 - <xs:complexType name="base_preparation_type">
                # XSD: <xs:element name="staining" minOccurs="0"> has 3 elements
                stain_in = spec_prep_in.get_staining()
                if stain_in is not None:
                    stain = emdb30.stainingType()
                    # element 1 - <xs:element name="staining" minOccurs="0">
                    # XSD: <xs:element name="type">
                    # Assume negative staining
                    stain.set_type('NEGATIVE')
                    # element 2 - <xs:element name="staining" minOccurs="0">
                    # XSD: <xs:element name="material" type="xs:token">
                    # element 3 - <xs:element name="staining" minOccurs="0">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                    stain.set_details(stain_in)
                    prep.set_staining(stain)
                # element 4 - <xs:complexType name="base_preparation_type">
                # XSD: <xs:element name="sugar_embedding" minOccurs="0">
                # element 5 - <xs:complexType name="base_preparation_type">
                # XSD: <xs:element name="shadowing" minOccurs="0">
                # element 6 - <xs:complexType name="base_preparation_type">
                # XSD: <xs:element name="grid" type="grid_type" minOccurs="0"/>
                grid_in = spec_prep_in.get_specimenSupportDetails()
                if grid_in is not None:
                    # XSD: <xs:complexType name="grid_type"> has 6 elements
                    grid = emdb30.grid_type()
                    # element 1 - <xs:complexType name="grid_type">
                    # XSD: <xs:element name="model" type="xs:token" minOccurs="0">
                    # element 2 - <xs:complexType name="grid_type">
                    # XSD: <xs:element name="material" minOccurs="0">
                    # element 3 - <xs:complexType name="grid_type">
                    # XSD: <xs:element name="mesh" type="xs:positiveInteger" minOccurs="0">
                    # element 4 - <xs:complexType name="grid_type">
                    # XSD: <xs:element name="support_film" type="film_type" maxOccurs="unbounded" minOccurs="0">
                    # element 5 - <xs:complexType name="grid_type">
                    # XSD: <xs:element name="pretreatment" type="grid_pretreatment_type" minOccurs="0"/>
                    # element 6 - <xs:complexType name="grid_type">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                    grid.set_details(grid_in)
                    prep.set_grid(grid)
            # element 7 - <xs:complexType name="base_preparation_type">
            # XSD: <xs:element name="vitrification" type="vitrification_type" minOccurs="0">
            if vitr_in is not None:
                v_in = vitr_in[i]
                # XSD: <xs:complexType name="vitrification_type"> has 7 elements
                vitr = emdb30.vitrification_type()
                # element 1 - <xs:complexType name="vitrification_type">
                # XSD: <xs:element name="cryogen_name">
                cryo_name = v_in.get_cryogenName()
                if cryo_name is not None:
                    if self.roundtrip:
                        vitr.set_cryogen_name(cryo_name)
                    else:
                        allowed_cryo_names = ['ETHANE', 'ETHANE-PROPANE MIXTURE', 'METHANE', 'NITROGEN', 'HELIUM', 'PROPANE', 'FREON 12', 'FREON 22', 'NONE', 'OTHER']
                        if cryo_name in allowed_cryo_names:
                            vitr.set_cryogen_name(cryo_name)
                        else:
                            vitr.set_cryogen_name('OTHER')
                # element 2 - <xs:complexType name="vitrification_type">
                # XSD: <xs:element name="chamber_humidity" minOccurs="0">
                humidity_in = v_in.get_humidity()
                if humidity_in is not None:
                    vitr.set_chamber_humidity(emdb30.chamber_humidityType(valueOf_=humidity_in, units=const.U_PERCENTAGE))
                # element 3 - <xs:complexType name="vitrification_type">
                # XSD: <xs:element name="chamber_temperature" minOccurs="0">
                temperature_in = v_in.get_temperature()
                if temperature_in is not None:
                    vitr.set_chamber_temperature(emdb30.chamber_temperatureType(valueOf_=temperature_in.get_valueOf_(), units=const.U_KEL))
                # element 4 - <xs:complexType name="vitrification_type">
                # XSD: <xs:element name="instrument" minOccurs="0">
                vitr_instrument = v_in.get_instrument()
                # allowed_vitr_instruments = ['FEI VITROBOT MARK I', 'FEI VITROBOT MARK II', 'FEI VITROBOT MARK III', 'FEI VITROBOT MARK IV', 'GATAN CRYOPLUNGE 3', 'HOMEMADE PLUNGER', 'LEICA EM CPC', 'LEICA EM GP', 'LEICA KF80', 'LEICA PLUNGER', 'REICHERT-JUNG PLUNGER', 'OTHER']
                # if vitr_instrument in allowed_vitr_instruments:
                vitr.set_instrument(vitr_instrument)
                #else:
                # vitr.set_instrument('OTHER')
                # element 5 - <xs:complexType name="vitrification_type">
                # XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                self.check_set(v_in.get_details, vitr.set_details)
                # element 6 - <xs:complexType name="vitrification_type">
                # XSD: <xs:element name="timed_resolved_state" type="xs:token" minOccurs="0">
                self.check_set(v_in.get_timeResolvedState, vitr.set_timed_resolved_state)
                # element 7 - <xs:complexType name="vitrification_type">
                # XSD: <xs:element name="method" type="xs:string" minOccurs="0">
                self.check_set(v_in.get_method, vitr.set_method)
                prep.set_vitrification(vitr)
            # element 8 - <xs:complexType name="base_preparation_type">
            # XSD: <xs:element name="details" type="xs:string" minOccurs="0">
            # crystalGrowDetails added here
            if spec_prep_in is not None:
                cryst_grow_details = spec_prep_in.get_crystalGrowDetails()
                if cryst_grow_details is not None and cryst_grow_details != '':
                    prep.set_details('crystalGrowDetails: %s :crystalGrowDetails' % cryst_grow_details)

        def set_base_microscopy(mic, img, im_ac_in):
            """
            Method that sets the elements common to all microscopy types

            Parameters:
            @params: mic - microscopy object
            """
            # XSD: <xs:complexType name="base_microscopy_type"> has 26 elements and 1 attribute
            # attribute 1 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:attribute name="id" type="xs:positiveInteger" use="required">
            mic.set_microscopy_id(i)
            # element 1 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="specimen_preparations" minOccurs="0">
            # element 2 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="microscope">
            mic.set_microscope(img.get_microscope())
            # element 3 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="illumination_mode">
            mic.set_illumination_mode(img.get_illuminationMode())
            # element 4 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="imaging_mode">
            mic.set_imaging_mode(img.get_imagingMode())
            # element 5 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="electron_source">
            mic.set_electron_source(img.get_electronSource())
            # element 6 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="acceleration_voltage">
            acc_voltage_in = img.get_acceleratingVoltage()
            if acc_voltage_in is not None:
                mic.set_acceleration_voltage(emdb30.acceleration_voltageType(valueOf_=acc_voltage_in.valueOf_, units=const.U_KVOLT))
            # element 7 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="c2_aperture_diameter" minOccurs="0">
            # element 8 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="nominal_cs" minOccurs="0">
            nominal_cs_in = img.get_nominalCs()
            if nominal_cs_in is not None:
                mic.set_nominal_cs(emdb30.nominal_csType(valueOf_=nominal_cs_in.valueOf_, units=const.U_MM))
            # element 9 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="nominal_defocus_min" minOccurs="0">
            add_to_details_nom_defocus_min = ''
            nominal_defoc_min = img.get_nominalDefocusMin()
            if nominal_defoc_min is not None:
                units_def_min = nominal_defoc_min.get_units()
                if units_def_min is not None:
                    if units_def_min == const.U_NM:
                        nominal_defoc_min_val = nominal_defoc_min.valueOf_
                        if nominal_defoc_min_val.find('.') == -1:
                            # the number is integer; write this into details as a note for round trip
                            add_to_details_nom_defocus_min = '{nominal defocus min is int}'
                        mic.set_nominal_defocus_min(emdb30.nominal_defocus_minType(valueOf_=float(nominal_defoc_min_val) * 0.001, units=const.U_MICROM))
                    elif units_def_min == const.U_MICROM:
                        mic.set_nominal_defocus_min(emdb30.nominal_defocus_minType(valueOf_=nominal_defoc_min_val, units=const.U_MICROM))
            # element 10 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="calibrated_defocus_min" minOccurs="0">
            # element 11 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="nominal_defocus_max" minOccurs="0">
            add_to_details_nom_defocus_max = ''
            nominal_defoc_max = img.get_nominalDefocusMax()
            if nominal_defoc_max is not None:
                units_def_max = nominal_defoc_max.get_units()
                if units_def_max is not None:
                    if units_def_max == const.U_NM:
                        nominal_defoc_max_val = nominal_defoc_max.valueOf_
                        if nominal_defoc_max_val.find('.') == -1:
                            # the number is integer; write this into details as a note for round trip
                            add_to_details_nom_defocus_max = '{nominal defocus max is int}'
                        mic.set_nominal_defocus_max(emdb30.nominal_defocus_maxType(valueOf_=float(nominal_defoc_max_val) * 0.001, units=const.U_MICROM))
                    else:
                        mic.set_nominal_defocus_max(emdb30.nominal_defocus_maxType(valueOf_=nominal_defoc_max_val, units=const.U_MICROM))
            # element 12 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="calibrated_defocus_max" minOccurs="0">
            # element 13 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="nominal_magnification" type="allowed_magnification" minOccurs="0">
            self.check_set(img.get_nominalMagnification, mic.set_nominal_magnification)
            # element 14 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="calibrated_magnification" type="allowed_magnification" minOccurs="0">
            self.check_set(img.get_calibratedMagnification, mic.set_calibrated_magnification)
            # element 15 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="specimen_holder_model" minOccurs="0">
            self.check_set(img.get_specimenHolderModel, mic.set_specimen_holder_model)
            # element 16 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="cooling_holder_cryogen" minOccurs="0">
            # element 17 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="temperature" minOccurs="0"> has 3 elements
            temp = emdb30.temperatureType()
            # element 1 - <xs:element name="temperature" minOccurs="0">
            # XSD: <xs:element name="temperature_min" type="temperature_type" minOccurs="0">
            temp_min_in = img.get_temperatureMin()
            if temp_min_in is not None:
                temp.set_temperature_min(emdb30.temperature_type(valueOf_=temp_min_in.get_valueOf_(), units=const.U_KEL))
            # element 2 - <xs:element name="temperature" minOccurs="0">
            # XSD: <xs:element name="temperature_max" type="temperature_type" minOccurs="0">
            temp_max_in = img.get_temperatureMax()
            if temp_max_in is not None:
                temp.set_temperature_max(emdb30.temperature_type(valueOf_=temp_max_in.get_valueOf_(), units=const.U_KEL))
            # element 3 - <xs:element name="temperature" minOccurs="0">
            # XSD: <xs:element name="temperature_average" type="temperature_type" minOccurs="0">
            temp_av_in = img.get_temperature()
            if temp_av_in is not None:
                temp.set_temperature_average(emdb30.temperature_type(valueOf_=temp_av_in.get_valueOf_(), units=const.U_KEL))
            mic.set_temperature(temp)
            # element 18 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="alignment_procedure" minOccurs="0"> has 1 element of 6 choices
            align = emdb30.alignment_procedureType()
            # choice 1 - <xs:element name="alignment_procedure" minOccurs="0">
            # XSD: <xs:element name="none">
            # choice 2 - <xs:element name="alignment_procedure" minOccurs="0">
            # XSD: <xs:element name="basic">
            # choice 3 - <xs:element name="alignment_procedure" minOccurs="0">
            # XSD: <xs:element name="zemlin_tableau">
            # choice 4 - <xs:element name="alignment_procedure" minOccurs="0">
            # XSD: <xs:element name="coma_free">
            # choice 5 - <xs:element name="alignment_procedure" minOccurs="0">
            # XSD: <xs:element name="other">
            # choice 6 - <xs:element name="alignment_procedure" minOccurs="0">
            # XSD: <xs:element name="legacy"> has 2 elements
            leg = emdb30.legacyType()
            # element 1 - <xs:element name="legacy">
            # XSD: <xs:element name="astigmatism" type="xs:string"/>
            ast_in = img.get_astigmatism()
            if ast_in is not None:
                leg.set_astigmatism(ast_in)
            # element 2 - <xs:element name="legacy">
            # XSD: <xs:element name="electron_beam_tilt_params" type="xs:string"/>
            tilt_in = img.get_electronBeamTiltParams()
            if tilt_in is not None:
                leg.set_electron_beam_tilt_params(tilt_in)
            if leg.hasContent_():
                align.set_legacy(leg)
            if align.hasContent_():
                mic.set_alignment_procedure(align)
            # element 19 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="specialist_optics" type="specialist_optics_type" minOccurs="0"/>
            # XSD: <xs:complexType name="specialist_optics_type"> has 4 elements
            spop = emdb30.specialist_optics_type()
            # element 1 - <xs:complexType name="specialist_optics_type">
            # XSD: <xs:element name="phase_plate" type="xs:token" minOccurs="0"/>
            # element 2 - <xs:complexType name="specialist_optics_type">
            # XSD: <xs:element name="sph_aberration_corrector" type="xs:token" minOccurs="0"/>
            # element 3 - <xs:complexType name="specialist_optics_type">
            # XSD: <xs:element name="chr_aberration_corrector" type="xs:token" minOccurs="0"/>
            # element 4 - <xs:complexType name="specialist_optics_type">
            # XSD: <xs:element name="energy_filter" minOccurs="0"> has 3 elements
            egf = emdb30.energy_filterType()
            # element 1 - <xs:element name="energy_filter" minOccurs="0">
            # XSD: <xs:element name="name" type="xs:token" minOccurs="0">
            egf_in = img.get_energyFilter()
            if egf_in is not None:
                egf.set_name(egf_in)
            # elements 2 and 3
            energy_window_in = img.get_energyWindow()
            is_negative = False
            if energy_window_in is not None:
                e_val_in = energy_window_in.get_valueOf_()
                if len(e_val_in) > 0:
                    # check if the first value is a negative number
                    if e_val_in[0] == '-':
                        is_negative = True
                e_units_in = energy_window_in.get_units()
                energy_window_values = e_val_in.split('-')
                len_energy_window_values = len(energy_window_values)
                if len_energy_window_values >= 2:
                    # element 2 - <xs:element name="energy_filter" minOccurs="0">
                    # XSD: <xs:element name="lower_energy_threshold" minOccurs="0">
                    # element 3 - <xs:element name="energy_filter" minOccurs="0">
                    # XSD: <xs:element name="upper_energy_threshold" minOccurs="0">
                    if is_negative:
                        egf.set_lower_energy_threshold(emdb30.lower_energy_thresholdType(valueOf_=float(energy_window_values[len_energy_window_values - 2]) * (-1), units=e_units_in))
                    else:
                        egf.set_lower_energy_threshold(emdb30.lower_energy_thresholdType(valueOf_=float(energy_window_values[len_energy_window_values - 2]), units=e_units_in))
                    egf.set_upper_energy_threshold(emdb30.upper_energy_thresholdType(valueOf_=float(energy_window_values[len_energy_window_values - 1]), units=e_units_in))
                else:
                    # element 3 - <xs:element name="energy_filter" minOccurs="0">
                    # XSD: <xs:element name="upper_energy_threshold" minOccurs="0">
                    egf.set_upper_energy_threshold(emdb30.upper_energy_thresholdType(valueOf_=-1.0, units=e_units_in))
            if egf.hasContent_():
                spop.set_energy_filter(egf)
            if spop.hasContent_():
                mic.set_specialist_optics(spop)
            # element 20 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
            # element 21 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="details" type="xs:string" minOccurs="0">
            im_details = img.get_details()
            all_details = ''
            if im_details is not None:
                all_details = im_details + add_to_details_nom_defocus_min + add_to_details_nom_defocus_max
            elif add_to_details_nom_defocus_min != '' or add_to_details_nom_defocus_max != '':
                all_details = add_to_details_nom_defocus_min + add_to_details_nom_defocus_max
            else:
                all_details = None
            mic.set_details(all_details)
            # element 22 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="date" type="xs:date" minOccurs="0">
            date_in = img.get_date()
            if date_in is not None:
                if self.roundtrip:
                    known_date_issues = {'2001-19-12': '2001-12-19'}
                    if date_in in known_date_issues:
                        date_corr = known_date_issues.get(date_in)
                        date_corr_prs = dtp.parse(date_corr)
                        mic.set_date(date_corr_prs)
                    else:
                        try:
                            date_prs = dtp.parse(date_in)
                            mic.set_date(date_prs)
                        except:
                            self.warn(1, "Unrecognized date format: %s" % date_in)
                else:
                    try:
                        date_prs = dtp.parse(date_in)
                        mic.set_date(date_prs)
                    except:
                        self.warn(1, "Unrecognized date format: %s" % date_in)
            # element 23 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="image_recording_list"> has 1 element
            im_rec_list = emdb30.image_recording_listType()
            if img is not None and im_ac_in == []:
                im_rec = emdb30.image_recordingType()
                detector_in = img.get_detector()
                if detector_in is not None:
                    # can happen
                    # element 1 - <xs:element name="image_recording_list">
                    # XSD: <xs:element name="image_recording"> has 12 elements and 1 attribute
                    # attribute 1 - <xs:element name="image_recording">
                    # XSD: <xs:attribute name="id" type="xs:positiveInteger"/>
                    # element 1 - <xs:element name="image_recording">
                    # XSD: <xs:element name="film_or_detector_model"> extends <xs:extension base="allowed_film_or_detector_model"> + 1 attribute
                    fod = emdb30.film_or_detector_modelType()
                    detector_in = img.get_detector()
                    if detector_in is not None:
                        if self.roundtrip:
                            fod.set_valueOf_(detector_in)
                        else:
                            detect_in = detector_in.upper()
                            allowed_detectors = ['AGFA SCIENTA FILM', 'DIRECT ELECTRON DE-10 (5k x 4k)', 'DIRECT ELECTRON DE-12 (4k x 3k)', 'DIRECT ELECTRON DE-16 (4k x 4k)',
                                                 'DIRECT ELECTRON DE-20 (5k x 3k)', 'DIRECT ELECTRON DE-64 (8k x 8k)', 'FEI CETA (4k x 4k)', 'FEI EAGLE (2k x 2k)',
                                                 'FEI EAGLE (4k x 4k)', 'FEI FALCON I (4k x 4k)', 'FEI FALCON II (4k x 4k)', 'FEI FALCON III (4k x 4k)', 'GATAN K2 (4k x 4k)',
                                                 'GATAN K2 BASE (4k x 4k)', 'GATAN K2 IS (4k x 4k)', 'GATAN K2 QUANTUM (4k x 4k)', 'GATAN K2 SUMMIT (4k x 4k)', 'GATAN MULTISCAN',
                                                 'GATAN ORIUS SC1000 (4k x 2.7k)', 'GATAN ORIUS SC200 (2k x 2k)', 'GATAN ORIUS SC600 (2.7k x 2.7k)',
                                                 'GATAN ULTRASCAN 1000 (2k x 2k)', 'GATAN ULTRASCAN 10000 (10k x 10k)', 'GATAN ULTRASCAN 4000 (4k x 4k)', 'GENERIC CCD',
                                                 'GENERIC CCD (2k x 2k)', 'GENERIC CCD (4k x 4k)', 'GENERIC FILM', 'GENERIC GATAN', 'GENERIC GATAN (2k x 2k)',
                                                 'GENERIC GATAN (4k x 4k)', 'GENERIC IMAGE PLATES', 'GENERIC TVIPS', 'GENERIC TVIPS (2k x 2k)', 'GENERIC TVIPS (4k x 4k)',
                                                 'KODAK 4489 FILM', 'KODAK SO-163 FILM', 'OTHER', 'PROSCAN TEM-PIV (2k x 2k)', 'SIA 15C (3k x 3k)', 'TVIPS TEMCAM-F216 (2k x 2k)',
                                                 'TVIPS TEMCAM-F224 (2k x 2k)', 'TVIPS TEMCAM-F415 (4k x 4k)', 'TVIPS TEMCAM-F416 (4k x 4k)', 'TVIPS TEMCAM-F816 (8k x 8k)']
                            if detect_in in allowed_detectors:
                                fod.set_valueOf_(detect_in)
                            else:
                                fod.set_valueOf_('OTHER')
                    if fod.hasContent_():
                        im_rec.set_film_or_detector_model(fod)
                # element 8 - <xs:element name="image_recording">
                # XSD: <xs:element name="average_electron_dose_per_image" minOccurs="0">
                dose_in = img.get_electronDose()
                if dose_in is not None:
                    im_rec.set_average_electron_dose_per_image(emdb30.average_electron_dose_per_imageType(valueOf_=dose_in.get_valueOf_(), units=const.U_EOVERANGSQR))
                # element 9 - <xs:element name="image_recording">
                # XSD: <xs:element name="detector_distance" type="xs:string" minOccurs="0"/>
                self.check_set(img.get_detectorDistance, im_rec.set_detector_distance)
                if im_rec.hasContent_():
                    im_rec_list.add_image_recording(im_rec)
                if im_rec_list.hasContent_():
                    mic.set_image_recording_list(im_rec_list)
            # if number of detectors > number of microscopes
            # mic: 1 2 3 4
            # det: 1 2 3 4,5,6
            # else
            # mic: 1 2 3 4
            # det: 1 2 2 2
            if num_det > 0:
                if num_det >= num_imaging_in:
                    min_idx = i - 1
                    if i == num_imaging_in:
                        max_idx = num_det
                    else:
                        max_idx = min_idx + 1
                else:
                    if i < num_det:
                        min_idx = i - 1
                    else:
                        min_idx = num_det - 1
                    max_idx = min_idx + 1

                for im_ac in im_ac_in[min_idx:max_idx]:
                    # element 1 - <xs:element name="image_recording_list">
                    # XSD: <xs:element name="image_recording"> has 12 elements and 1 attribute
                    im_rec = emdb30.image_recordingType()
                    # attribute 1 - <xs:element name="image_recording">
                    # XSD: <xs:attribute name="id" type="xs:positiveInteger"/>
                    # element 1 - <xs:element name="image_recording">
                    # XSD: <xs:element name="film_or_detector_model"> extends <xs:extension base="allowed_film_or_detector_model"> + 1 attribute
                    fod = emdb30.film_or_detector_modelType()
                    detector_in = img.get_detector()
                    if detector_in is not None:
                        if self.roundtrip:
                            fod.set_valueOf_(detector_in)
                        else:
                            detect_in = detector_in.upper()
                            allowed_detectors = ['AGFA SCIENTA FILM', 'DIRECT ELECTRON DE-10 (5k x 4k)', 'DIRECT ELECTRON DE-12 (4k x 3k)', 'DIRECT ELECTRON DE-16 (4k x 4k)',
                                                 'DIRECT ELECTRON DE-20 (5k x 3k)', 'DIRECT ELECTRON DE-64 (8k x 8k)', 'FEI CETA (4k x 4k)', 'FEI EAGLE (2k x 2k)',
                                                 'FEI EAGLE (4k x 4k)', 'FEI FALCON I (4k x 4k)', 'FEI FALCON II (4k x 4k)', 'FEI FALCON III (4k x 4k)', 'GATAN K2 (4k x 4k)',
                                                 'GATAN K2 BASE (4k x 4k)', 'GATAN K2 IS (4k x 4k)', 'GATAN K2 QUANTUM (4k x 4k)', 'GATAN K2 SUMMIT (4k x 4k)', 'GATAN MULTISCAN',
                                                 'GATAN ORIUS SC1000 (4k x 2.7k)', 'GATAN ORIUS SC200 (2k x 2k)', 'GATAN ORIUS SC600 (2.7k x 2.7k)',
                                                 'GATAN ULTRASCAN 1000 (2k x 2k)', 'GATAN ULTRASCAN 10000 (10k x 10k)', 'GATAN ULTRASCAN 4000 (4k x 4k)', 'GENERIC CCD',
                                                 'GENERIC CCD (2k x 2k)', 'GENERIC CCD (4k x 4k)', 'GENERIC FILM', 'GENERIC GATAN', 'GENERIC GATAN (2k x 2k)',
                                                 'GENERIC GATAN (4k x 4k)', 'GENERIC IMAGE PLATES', 'GENERIC TVIPS', 'GENERIC TVIPS (2k x 2k)', 'GENERIC TVIPS (4k x 4k)',
                                                 'KODAK 4489 FILM', 'KODAK SO-163 FILM', 'OTHER', 'PROSCAN TEM-PIV (2k x 2k)', 'SIA 15C (3k x 3k)', 'TVIPS TEMCAM-F216 (2k x 2k)',
                                                 'TVIPS TEMCAM-F224 (2k x 2k)', 'TVIPS TEMCAM-F415 (4k x 4k)', 'TVIPS TEMCAM-F416 (4k x 4k)', 'TVIPS TEMCAM-F816 (8k x 8k)']
                            if detect_in in allowed_detectors:
                                fod.set_valueOf_(detect_in)
                            else:
                                fod.set_valueOf_('OTHER')
                    # attribute 1 - <xs:element name="film_or_detector_model">
                    # XSD: <xs:attribute name="category">
                    scanner_in = im_ac.get_scanner()
                    if scanner_in is not None:
                        fod.set_category('FILM')
                    else:
                        # For now classify all as CCD - this may need remediation
                        fod.set_category('CCD')

                    if fod.hasContent_():
                        im_rec.set_film_or_detector_model(fod)
                    # element 2 - <xs:element name="image_recording">
                    # XSD: <xs:element name="detector_mode" minOccurs="0">
                    # element 3 - <xs:element name="image_recording">
                    # XSD: <xs:element name="digitization_details" minOccurs="0"> has 4 elements
                    dig = emdb30.digitization_detailsType()
                    # element 1 - <xs:element name="digitization_details" minOccurs="0">
                    # XSD: <xs:element name="scanner" minOccurs="0">
                    if scanner_in is not None:
                        dig.set_scanner(scanner_in)
                    # element 2 - <xs:element name="digitization_details" minOccurs="0">
                    # XSD: <xs:element name="dimensions" minOccurs="0">
                    # element 3 - <xs:element name="digitization_details" minOccurs="0">
                    # XSD: <xs:element name="sampling_interval" minOccurs="0">
                    sampling_size_in = im_ac.get_samplingSize()
                    if sampling_size_in is not None:
                        dig.set_sampling_interval(emdb30.sampling_intervalType(valueOf_=sampling_size_in.get_valueOf_(), units=const.U_MICROM))
                    # element 4 - <xs:element name="digitization_details" minOccurs="0">
                    # XSD: <xs:element name="frames_per_image" type="xs:token" minOccurs="0">
                    if dig.hasContent_():
                        im_rec.set_digitization_details(dig)
                    # element 4 - <xs:element name="image_recording">
                    # XSD: <xs:element name="number_grids_imaged" type="xs:positiveInteger" minOccurs="0"/>
                    # element 5 - <xs:element name="image_recording">
                    # XSD: <xs:element name="number_real_images" type="xs:positiveInteger" minOccurs="0"/>
                    self.check_set(im_ac.get_numDigitalImages, im_rec.set_number_real_images)
                    # element 6 - <xs:element name="image_recording">
                    # XSD: <xs:element name="number_diffraction_images" type="xs:positiveInteger" minOccurs="0"/>
                    # element 7 - <xs:element name="image_recording">
                    # XSD: <xs:element name="average_exposure_time" minOccurs="0">
                    # element 8 - <xs:element name="image_recording">
                    # XSD: <xs:element name="average_electron_dose_per_image" minOccurs="0">
                    dose_in = img.get_electronDose()
                    if dose_in is not None:
                        im_rec.set_average_electron_dose_per_image(emdb30.average_electron_dose_per_imageType(valueOf_=dose_in.get_valueOf_(), units=const.U_EOVERANGSQR))
                    # element 9 - <xs:element name="image_recording">
                    # XSD: <xs:element name="detector_distance" type="xs:string" minOccurs="0"/>
                    self.check_set(img.get_detectorDistance, im_rec.set_detector_distance)
                    # element 10 - <xs:element name="image_recording">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    self.check_set(im_ac.get_details, im_rec.set_details)
                    # element 11 - <xs:element name="image_recording">
                    # XSD: <xs:element name="od_range" type="xs:float" minOccurs="0">
                    self.check_set(im_ac.get_odRange, im_rec.set_od_range)
                    # element 12 - <xs:element name="image_recording">
                    # XSD: <xs:element name="bits_per_pixel" type="xs:float" minOccurs="0">
                    self.check_set(im_ac.get_quantBitNumber, im_rec.set_bits_per_pixel)

                    if im_rec.hasContent_():
                        im_rec_list.add_image_recording(im_rec)

                if im_rec_list.hasContent_():
                    mic.set_image_recording_list(im_rec_list)
            # element 24 - <xs:complexType name="base_microscopy_type">
            # XSD: <xs:element name="specimen_holder" type="xs:string" minOccurs="0">
            self.check_set(img.get_specimenHolder, mic.set_specimen_holder)
            # elements 25 and 26
            tilt_min_in = img.get_tiltAngleMin()
            tilt_max_in = img.get_tiltAngleMax()
            if not isinstance(mic, emdb30.tomography_microscopy_type):
                if tilt_min_in is not None:
                    # element 25 - <xs:complexType name="base_microscopy_type">
                    # XSD: <xs:element name="tilt_angle_min" minOccurs="0">
                    tilt_min = tilt_min_in.valueOf_
                    mic.set_tilt_angle_min(tilt_min)
                if tilt_max_in is not None:
                    # element 26 - <xs:complexType name="base_microscopy_type">
                    # XSD:<xs:element name="tilt_angle_max" minOccurs="0">
                    tilt_max = tilt_max_in.valueOf_
                    mic.set_tilt_angle_max(tilt_max)

        def set_tilt_series(mic, img, tom_proc=None):
            """
            Method that sets min, max and increment of an angle in the tilt series

            @param: mic - microscopy object v3.0
            @param: img - imaging from v1.9
            @params: tom_proc - tomography processing v1.9
            """
            tilt_series = emdb30.tilt_series_type()
            # element 1 - <xs:complexType name="tilt_series_type">
            # XSD: <xs:element name="axis1" type="axis_type"/>
            axis1 = emdb30.axis_type()
            # XSD: <xs:complexType name="axis_type"> has 3 elements
            # element 1 - <xs:complexType name="axis_type">
            # XSD: <xs:element name="min_angle" minOccurs="0">
            tilt_min_in = img.get_tiltAngleMin()
            if tilt_min_in is not None:
                tilt_min_value = tilt_min_in.get_valueOf_()
                tilt_min = emdb30.min_angleType(valueOf_=tilt_min_value, units=const.U_DEG)
                axis1.set_min_angle(tilt_min)
            # element 2 - <xs:complexType name="axis_type">
            # XSD: <xs:element name="max_angle" minOccurs="0">
            tilt_max_in = img.get_tiltAngleMax()
            if tilt_max_in is not None:
                tilt_max_value = tilt_max_in.get_valueOf_()
                tilt_max = emdb30.max_angleType(valueOf_=tilt_max_value, units=const.U_DEG)
                axis1.set_max_angle(tilt_max)
            # element 3 - <xs:complexType name="axis_type">
            # XSD: <xs:element name="angle_increment" minOccurs="0">
            if tom_proc is not None:
                tilt_inc_in = tom_proc.get_tiltAngleIncrement()
                if tilt_inc_in is not None:
                    tilt_inc = emdb30.angle_incrementType(valueOf_=tilt_inc_in, units=const.U_DEG)
                    axis1.set_angle_increment(tilt_inc)

            tilt_series.set_axis1(axis1)
            # element 2 - <xs:complexType name="tilt_series_type">
            # XSD: <xs:element name="axis2" minOccurs="0">
            # element 3 - <xs:complexType name="tilt_series_type">
            # XSD: <xs:element name="axis_rotation" fixed="90" minOccurs="0"
            mic.add_tilt_series(tilt_series)

        xml_in = emdb_19.parse(input_file, silence=True)
        # XSD: <xs:complexType name="entry_type"> has 7 elements and 2 attributes
        xml_out = emdb30.entry_type()
        # attribute 1 - <xs:complexType name="entry_type">
        # XSD: <xs:attribute name="emdb_id" type="emdb_id_type" use="required">
        code_in = xml_in.get_accessCode()
        fmt_code_in = self.format_emdb_code(code_in)
        xml_out.set_emdb_id(fmt_code_in)
        # attribute 2 - <xs:complexType name="entry_type">
        # XSD: <xs:attribute name="version" use="required">
        xml_out.set_version(const.LATEST_XML_VERSION)
        # element 1 - <xs:complexType name="entry_type">
        # XSD: <xs:element name="admin" type="admin_type">
        # XSD: <<xs:complexType name="admin_type"> has 13 elements
        admin = emdb30.admin_type()
        # element 1 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="status_history_list" type="version_list_type" minOccurs="0">
        # XSD: <xs:complexType name="version_list_type"> has 1 element
        status_history_list = emdb30.version_list_type()
        # element 1 - <xs:complexType name="version_list_type">
        # XSD: <xs:element name="status" maxOccurs="unbounded"> extension of base="version_type" + 1 attribute
        dep_in = xml_in.get_deposition()
        status_prior_in = dep_in.get_status().get_prior()
        if status_prior_in is not None:
            prior_status = emdb30.statusType()
            # <xs:complexType name="version_type"> has 5 elements
            # element 1 - <xs:complexType name="version_type">
            # XSD: <xs:element name="date" minOccurs="0">
            # element 2 - <xs:complexType name="version_type">
            # XSD: <xs:element name="code" type="code_type">
            prior_status.set_code(emdb30.code_type(valueOf_=status_prior_in))
            # element 3 - <xs:complexType name="version_type">
            # XSD: <xs:element name="processing_site" minOccurs="0">
            # element 4 - <xs:complexType name="version_type">
            # XSD: <xs:element name="annotator" minOccurs="0">
            # element 5 - <xs:complexType name="version_type">
            # XSD: <xs:element name="details" type="xs:string" minOccurs="0">
            # attribute 1 - <xs:element name="status" maxOccurs="unbounded">
            # XSD: <xs:attribute name="id" type="xs:positiveInteger" use="required"/>
            prior_status.set_id(1)
            status_history_list.add_status(prior_status)
        if status_history_list.hasContent_():
            admin.set_status_history_list(status_history_list)
        # element 2 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="current_status" type="version_type">
        # XSD: <xs:complexType name="version_type"> has 5 elements
        current_status = emdb30.version_type()
        # element 1 - <xs:complexType name="version_type">
        # XSD: <xs:element name="date" minOccurs="0">
        # element 2 - <xs:complexType name="version_type">
        # XSD: <xs:element name="code" type="code_type">
        # XSD: <xs:complexType name="code_type"> has base base="status_code_type" + 2 attributes
        code = emdb30.code_type()
        # XSD: <xs:simpleType name="status_code_type"> <xs:restriction base="xs:token">
        code.set_valueOf_(dep_in.get_status().get_valueOf_())
        # attribute 1 - <xs:complexType name="code_type">
        # XSD: <xs:attribute name="superseded" type="xs:boolean">
        superseded_list_in = dep_in.get_supersededByList()
        if superseded_list_in is not None:
            code.set_superseded(True)
        # attribute 2 - <xs:complexType name="code_type">
        # XSD: xs:attribute name="supersedes" type="xs:boolean">
        obs_list_in = dep_in.get_obsoleteList()
        if obs_list_in is not None:
            code.set_supersedes(True)

        current_status.set_code(code)
        # element 3 - <xs:complexType name="version_type">
        # XSD: <xs:element name="processing_site" minOccurs="0">
        current_status.set_processing_site(dep_in.get_processingSite())
        # element 4 - <xs:complexType name="version_type">
        # XSD: <xs:element name="annotator" minOccurs="0">
        # element 5 - <xs:complexType name="version_type">
        # XSD: <xs:element name="details" type="xs:string" minOccurs="0">
        if current_status.hasContent_():
            admin.set_current_status(current_status)
        # element 3 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="sites"> has 2 elements
        sites = emdb30.sitesType()
        # element 1 - <xs:element name="sites">
        # XSD: <xs:element name="deposition">
        sites.set_deposition(dep_in.get_depositionSite())
        # element 2 - <xs:element name="sites">
        # XSD: <xs:element name="last_processing">
        sites.set_last_processing(dep_in.get_processingSite())
        admin.set_sites(sites)
        # element 4 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="key_dates"> has 5 elements
        key_dates = emdb30.key_datesType()
        # element 1 - <xs:element name="key_dates">
        # XSD: <xs:element name="deposition" type="xs:date">
        key_dates.set_deposition(dep_in.get_depositionDate())
        # element 2 - <xs:element name="key_dates">
        # XSD: <xs:element name="header_release" type="xs:date" minOccurs="0">
        key_dates.set_header_release(dep_in.get_headerReleaseDate())
        # element 3 - <xs:element name="key_dates">
        # XSD: <xs:element name="map_release" type="xs:date" minOccurs="0">
        self.check_set(dep_in.get_mapReleaseDate, key_dates.set_map_release)
        # element 4 - <xs:element name="key_dates">
        # XSD: <xs:element name="obsolete" type="xs:date" minOccurs="0">
        self.check_set(dep_in.get_obsoletedDate, key_dates.set_obsolete)
        # element 5 - <xs:element name="key_dates">
        # XSD: <xs:element name="update" type="xs:date">
        key_dates.set_update(xml_in.get_admin().get_lastUpdate())
        admin.set_key_dates(key_dates)
        # element 5 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="obsolete_list" minOccurs="0"> has 1 element
        if obs_list_in is not None:
            obs_list = emdb30.obsolete_listType()
            obs_entries_in = obs_list_in.get_entry()
            # element 1 - <xs:element name="obsolete_list" minOccurs="0">
            # XSD: <xs:element name="entry" type="supersedes_type" maxOccurs="unbounded">
            for obs_in in obs_entries_in:
                obs = emdb30.supersedes_type()
                obs.set_entry(obs_in)
                obs_list.add_entry(obs)
            if obs_list.hasContent_():
                admin.set_obsolete_list(obs_list)
        # element 6 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="superseded_by_list" minOccurs="0"> has 1 element
        if superseded_list_in is not None:
            supersede_list = emdb30.superseded_by_listType()
            supersede_entries_in = superseded_list_in.get_entry()
            # element 1 - <xs:element name="superseded_by_list" minOccurs="0">
            # XSD: <xs:element name="entry" type="supersedes_type" maxOccurs="unbounded">
            for supersede_in in supersede_entries_in:
                supersede = emdb30.supersedes_type()
                supersede.set_entry(supersede_in)
                supersede_list.add_entry(supersede)
            if supersede_list.hasContent_():
                admin.set_superseded_by_list(supersede_list)
        # element 7 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="grant_support" minOccurs="0">
        # element 8 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="contact_author" maxOccurs="unbounded" minOccurs="0">
        # element 9 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="title" type="xs:token">
        admin.set_title(dep_in.get_title())
        # element 10 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="authors_list"> has 1 element
        author_list = emdb30.authors_listType()
        # element 1 - <xs:element name="authors_list">
        # XSD: <xs:element name="author" type="author_type" maxOccurs="unbounded">
        copy_authors(dep_in.get_authors, author_list.add_author, simple=True)
        # if author_list.hasContent_():
        admin.set_authors_list(author_list)
        # element 11 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="details" type="xs:token" minOccurs="0">
        # element 12 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="keywords" type="xs:string" minOccurs="0">
        kwrds = dep_in.get_keywords()
        if kwrds is not None:
            admin.set_keywords(kwrds)
        # element 13 - <xs:complexType name="admin_type">
        # XSD: <xs:element name="replace_existing_entry" type="xs:boolean" minOccurs="0">
        self.check_set(dep_in.get_replaceExistingEntry, admin.set_replace_existing_entry)

        xml_out.set_admin(admin)
        # element 2 - <xs:complexType name="entry_type">
        # XSD: <xs:element name="crossreferences" type="crossreferences_type"/>
        # XSD: <xs:complexType name="crossreferences_type"> has 4 elements
        cref = emdb30.crossreferences_type()
        # element 1 - <xs:complexType name="crossreferences_type">
        # XSD: <xs:element name="citation_list"> has 2 elements
        cite_list = emdb30.citation_listType()
        # element 1 - <xs:element name="citation_list">
        # XSD: <xs:element name="primary_citation"> has 1 element
        # element 1 - <xs:element name="primary_citation">
        ref_out = emdb30.primary_citationType()
        # XSD: <xs:element ref="citation_type"/>
        ref_in = dep_in.get_primaryReference()
        copy_citation(ref_in, ref_out)
        cite_list.set_primary_citation(ref_out)
        # element 2 - <xs:element name="citation_list">
        # XSD: <xs:element name="secondary_citation" maxOccurs="unbounded" minOccurs="0"> has 1 element
        refs_in = dep_in.get_secondaryReference()
        for ref_in in refs_in:
            ref_out = emdb30.secondary_citationType()
            # element 1 - <xs:element name="secondary_citation" maxOccurs="unbounded" minOccurs="0">
            # XSD: <xs:element ref="citation_type"/>
            cite_list.add_secondary_citation(ref_out)
            copy_citation(ref_in, ref_out)
        cref.set_citation_list(cite_list)
        # element 2 - <xs:complexType name="crossreferences_type">
        # XSD: <xs:element name="emdb_list" type="emdb_cross_reference_list_type" minOccurs="0"> has 1 element
        emdb_in = dep_in.get_inFrameEMDBId()
        if emdb_in is not None and emdb_in != '':
            emdb_list_in = emdb_in.replace(' ', '').split(',')
            if len(emdb_list_in) > 0:
                # element 1 - <xs:element name="emdb_list" type="emdb_cross_reference_list_type" minOccurs="0">
                # XSD: <xs:element name="emdb_reference" type="emdb_cross_reference_type" maxOccurs="unbounded"/>
                emdb_list = emdb30.emdb_cross_reference_list_type()
                map_in_file = 1
                for e_value in emdb_list_in:
                    # XSD: <xs:complexType name="emdb_cross_reference_type"> has 3 elements
                    emdb_elem = emdb30.emdb_cross_reference_type()
                    # element 1 - <xs:complexType name="emdb_cross_reference_type">
                    # XSD: <xs:element name="emdb_id" type="emdb_id_type"/>
                    emdb_elem.set_emdb_id(e_value)
                    # element 2 - <xs:complexType name="emdb_cross_reference_type">
                    # XSD: <xs:element name="relationship" minOccurs="0">
                    emdb_elem.set_relationship(emdb30.relationshipType('FULLOVERLAP'))
                    # element 3 - <xs:complexType name="emdb_cross_reference_type">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    # attribute 1 - not used any more
                    # XSD: <!--  <xs:attribute name="id" type="xs:positiveInteger"/>  -->
                    # emdb_elem.set_id(map_in_file)
                    emdb_list.add_emdb_reference(emdb_elem)
                    map_in_file += 1
                cref.set_emdb_list(emdb_list)
        # element 3 - <xs:complexType name="crossreferences_type">
        # XSD: <xs:element name="pdb_list" type="pdb_cross_reference_list_type" minOccurs="0">
        pdb_list_in = dep_in.get_fittedPDBEntryIdList()
        if pdb_list_in:
            pdbs_in = pdb_list_in.get_fittedPDBEntryId()
            if pdbs_in is not None and len(pdbs_in) > 0:
                # XSD: <xs:complexType name="pdb_cross_reference_list_type"> has 1 element
                pdb_list = emdb30.pdb_cross_reference_list_type()
                cref.set_pdb_list(pdb_list)
                # map_in_file = 1
                for p_in in pdbs_in:
                    # element 1 - <xs:element name="pdb_list" type="pdb_cross_reference_list_type" minOccurs="0">
                    # XSD: <xs:element name="pdb_reference" type="pdb_cross_reference_type" maxOccurs="unbounded">
                    # XSD: <xs:complexType name="pdb_cross_reference_type"> has 3 elements
                    pdb_elem = emdb30.pdb_cross_reference_type()
                    # element 1 - <xs:complexType name="pdb_cross_reference_type">
                    # XSD: <xs:element name="pdb_id" type="pdb_code_type"/>
                    pdb_elem.set_pdb_id(p_in)
                    # element 2 - <xs:complexType name="pdb_cross_reference_type">
                    # XSD: <xs:element name="relationship" minOccurs="0">
                    pdb_elem.set_relationship(emdb30.relationshipType('FULLOVERLAP'))
                    # element 3 - <xs:complexType name="pdb_cross_reference_type">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    # pdb_elem.set_details(???)
                    # pdb_elem.set_id(map_in_file)
                    pdb_list.add_pdb_reference(pdb_elem)
                    # map_in_file += 1
        # element 4 - <xs:complexType name="crossreferences_type">
        # XSD: <xs:element name="auxiliary_link_list" minOccurs="0">
        exp_in = xml_in.get_experiment()
        im_ac_in = exp_in.get_imageAcquisition()
        num_det = len(im_ac_in)
        if num_det == 0 and not self.roundtrip:
            self.warn(1, "No image acquisition elements found!")
        if num_det > 0:
            min_idx = 0
            max_idx = 0
            i = 1
            imaging_list_in = exp_in.get_imaging()
            num_imaging_in = len(imaging_list_in)
            if num_det >= num_imaging_in:
                min_idx = i - 1
                if i == num_imaging_in:
                    max_idx = num_det
                else:
                    max_idx = min_idx + 1
            else:
                if i < num_det:
                    min_idx = i - 1
                else:
                    min_idx = num_det - 1
                max_idx = min_idx + 1
            for im_ac in im_ac_in[min_idx:max_idx]:
                url_raw_data_in = im_ac.get_URLRawData()
                if url_raw_data_in is not None:
                    url_ref_list = emdb30.auxiliary_link_listType()
                    cref.set_auxiliary_link_list(url_ref_list)
                    url_ref = emdb30.auxiliary_link_type()
                    url_ref_list.add_auxiliary_link(url_ref)
                    url_ref.set_link(url_raw_data_in)

        xml_out.set_crossreferences(cref)

        def set_base_supramolecule(supmol, supmol_id, supmol_in, comp_in=None):
            """
            1.9 -> 3.0: Method that sets the elements that are common to all supramolecules.

            Parameters:
            @params: supmol - sample_supmol/virus_smol/comp_smol/complex_smol - object of a supramolecule type
            @params: supmol_in - sample_in/virus_in/cell_comp_in/complex_smol_in - v.19
            @params: comp_in - protein/cell-component etc
            @params: rib_cat - 'ribosome-eukaryote' or 'ribosome-prokaryote' used in setting the category element
            """
            # XSD: <xs:complexType name="base_supramolecule_type"> has 9 elements and 1 attribute
            if comp_in is not None:
                # all but sample
                # attribute 1 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:attribute name="id" type="xs:positiveInteger" use="required"/>
                supmol.set_supramolecule_id(supmol_id)
                # element 1 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="name" type="sci_name_type">
                set_sci_name(comp_in, supmol)
                # element 2 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="category" minOccurs="0">
                # element 3 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="parent" type="xs:nonNegativeInteger">
                # element 4 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="macromolecule_list" minOccurs="0">
                if supmol.original_tagname_ != 'virus_supramolecule':
                    if supmol_in is not None:
                        # element 5 - <xs:complexType name="base_supramolecule_type">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        self.check_set(supmol_in.get_details, supmol.set_details)
                        # element 6 - <xs:complexType name="base_supramolecule_type">
                        # XSD: <xs:element name="number_of_copies" type="pos_int_or_string_type" minOccurs="0">
                        self.check_set(supmol_in.get_numCopies, supmol.set_number_of_copies)
                        # element 7 - <xs:complexType name="base_supramolecule_type">
                        # XSD: <xs:element name="oligomeric_state" type="pos_int_or_string_type" minOccurs="0">
                        self.check_set(supmol_in.get_oligomericDetails, supmol.set_oligomeric_state)
                else:
                    # element 5 - <xs:complexType name="base_supramolecule_type">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    self.check_set(comp_in.get_details, supmol.set_details)
                # element 8 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                ext_refs = supmol_in.get_externalReferences()
                if ext_refs is not None and ext_refs != []:
                    add_mol_references(supmol.add_external_references, ext_refs)
                if supmol.original_tagname_ != 'virus_supramolecule':
                    # element 9 - <xs:complexType name="base_supramolecule_type">
                    # XSD: <xs:element name="recombinant_exp_flag" type="xs:boolean" maxOccurs="1" minOccurs="0">
                    self.check_set(supmol_in.get_recombinantExpFlag, supmol.set_recombinant_exp_flag)
            else:
                # sample - the sample supramolecule doesn't have components
                # attribute 1 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:attribute name="id" type="xs:positiveInteger" use="required"/>
                supmol.set_supramolecule_id(const.EM_SAMPLE_ID)
                # element 1 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="name" type="sci_name_type">
                supmol.set_name(emdb30.sci_name_type(valueOf_=supmol_in.get_name()))
                # element 2 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="category" minOccurs="0">
                # element 3 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="parent" type="xs:nonNegativeInteger">
                # element 4 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="macromolecule_list" minOccurs="0">
                # element 5 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                self.check_set(supmol_in.get_details, supmol.set_details)
                # element 6 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="number_of_copies" type="pos_int_or_string_type" minOccurs="0">
                # element 7 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="oligomeric_state" type="pos_int_or_string_type" minOccurs="0">
                self.check_set(supmol_in.get_compDegree, supmol.set_oligomeric_state)
                # element 8 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                # element 9 - <xs:complexType name="base_supramolecule_type">
                # XSD: <xs:element name="recombinant_exp_flag" type="xs:boolean" maxOccurs="1" minOccurs="0">

        # element 3 - <xs:complexType name="entry_type">
        # XSD: <xs:element name="sample" type="sample_type">

        def set_base_source(mol_nat_source, spec_component_in):
            """
            Method called from: def set_mol_natural_source(nat_source, spec_component_in, tissue=True, cell=True, organelle=True, cell_loc=True)
            Sets the elements that feature in all molecules needed to set natural source
            """
            species_in = spec_component_in.get_sciSpeciesName()
            if species_in is not None:
                # XSD: <xs:complexType name="base_source_type"> has 3 elements and 1 attribute
                # attribute 1 - <xs:complexType name="base_source_type">
                # XSD: <xs:attribute name="database" use="required">
                mol_nat_source.set_database('NCBI')
                # element 1 - <xs:complexType name="base_source_type">
                # XSD: <xs:element name="organism" type="organism_type">
                org = emdb30.organism_type()
                # XSD: <xs:complexType name="organism_type"> has 1 attribute and is ext of token
                self.check_set(species_in.get_valueOf_, org.set_valueOf_)
                # attribute 1 - <xs:complexType name="organism_type">
                # XSD: <xs:attribute name="ncbi" type="xs:positiveInteger"/>
                self.check_set(species_in.get_ncbiTaxId, org.set_ncbi)
                mol_nat_source.set_organism(org)
                # element 2 - <xs:complexType name="base_source_type">
                # <xs:element name="strain" type="xs:token" minOccurs="0"/>
                strain_in = spec_component_in.get_sciSpeciesStrain()
                if strain_in is not None:
                    # strain = emdb30.organism_type()
                    ## XSD: <xs:complexType name="organism_type"> has 1 attribute and is ext of token
                    #strain.set_valueOf_(strain_in)
                    ## attribute 1 - <xs:complexType name="organism_type">
                    ## XSD: <xs:attribute name="ncbi" type="xs:positiveInteger"/>
                    mol_nat_source.set_strain(strain_in)
                # element 3 - <xs:complexType name="base_source_type">
                # XSD: <xs:element name="synonym_organism" type="xs:token" minOccurs="0">
                self.check_set(spec_component_in.get_synSpeciesName, mol_nat_source.set_synonym_organism)

        def set_mol_natural_source(nat_source, component_in, spec_component_in, tissue=True, cell=True, organelle=True, cell_loc=True):
            """
            v19 to v3.0: Method that sets natural source for molecules from proteinType/cellCompType/virusType/nuclAcidType/ligandType/riboTypeEu/riboTypePro
            """
            # <xs:complexType name="macromolecule_source_type"> has a base and 5 elements
            # base - <xs:complexType name="macromolecule_source_type">
            # XSD: <xs:extension base="base_source_type">
            if spec_component_in is not None:
                set_base_source(nat_source, spec_component_in)
                # nuclAcidType and labelType (this method is never called for labels) do not have natSource
                if component_in.get_entry() != 'nucleic-acid':
                    ns_in = spec_component_in.get_natSource()
                    if ns_in is not None and ns_in.hasContent_():
                        # element 1 - <xs:complexType name="macromolecule_source_type">
                        # XSD: <xs:element name="organ" type="xs:token" minOccurs="0"/>
                        # element 2 - <xs:complexType name="macromolecule_source_type">
                        # XSD: <xs:element name="tissue" type="xs:token" minOccurs="0">
                        if tissue:
                            self.check_set(ns_in.get_organOrTissue, nat_source.set_tissue)
                        # element 3 - <xs:complexType name="macromolecule_source_type">
                        # XSD:<xs:element name="cell" type="xs:token" minOccurs="0">
                        if cell:
                            self.check_set(ns_in.get_cell, nat_source.set_cell)
                        # element 4 - <xs:complexType name="macromolecule_source_type">
                        # XSD: <xs:element name="organelle" type="xs:token" minOccurs="0">
                        if organelle:
                            self.check_set(ns_in.get_organelle, nat_source.set_organelle)
                        # element 5 - <xs:complexType name="macromolecule_source_type">
                        # XSD: <xs:element name="cellular_location" type="xs:token" minOccurs="0">
                        if cell_loc:
                            self.check_set(ns_in.get_cellLocation, nat_source.set_cellular_location)

        def set_base_macromolecule(mol, c_id, component_in, spec_component_in, nucleic_acid=False, label=False):
            """
            1.9 -> 3.0
            Parameters:

            @param mol:
            @param c_id: componenet ID. Count kept separately for macromolecules
            """
            # XSD: <xs:complexType name="base_macromolecule_type"> has 7 elements and 3 attributes
            # attribute 1 - <xs:complexType name="base_macromolecule_type">
            # XSD: <xs:attribute name="id" type="xs:positiveInteger" use="required">
            mol.set_macromolecule_id(c_id)
            # attribute 2 - <xs:complexType name="base_macromolecule_type">
            # XSD: <xs:attribute name="mutant" type="xs:boolean"/>
            # attribute 3 - <xs:complexType name="base_macromolecule_type">
            # XSD: <xs:attribute name="chimera" type="xs:boolean"/>
            # element 1 - <xs:complexType name="base_macromolecule_type">
            # XSD: <xs:element name="name" type="sci_name_type"/>
            if not label:
                set_sci_name(component_in, mol)
            else:
                name = emdb30.sci_name_type()
                sci_name = component_in.get_sciName()
                if sci_name is not None:
                    name.set_valueOf_(sci_name)
                mol.set_name(name)
            # element 2 - <xs:complexType name="base_macromolecule_type">
            # XSD: <xs:element name="natural_source" type="macromolecule_source_type" minOccurs="0"/>
            if not label:
                mol_nat_source = emdb30.macromolecule_source_type()
                set_mol_natural_source(mol_nat_source, component_in, spec_component_in)
                if mol_nat_source.hasContent_():
                    mol.set_natural_source(mol_nat_source)
            # element 3 - <xs:complexType name="base_macromolecule_type">
            # XSD: <xs:element name="molecular_weight" type="molecular_weight_type" minOccurs="0">
            set_mol_weight(mol.set_molecular_weight, component_in.get_molWtTheo(), component_in.get_molWtExp())
            # element 4 - <xs:complexType name="base_macromolecule_type">
            # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
            self.check_set(component_in.get_details, mol.set_details)

            if not nucleic_acid and spec_component_in is not None:
                # element 5 - <xs:complexType name="base_macromolecule_type">
                # XSD: <xs:element name="number_of_copies" type="pos_int_or_string_type" minOccurs="0"/>
                self.check_set(spec_component_in.get_numCopies, mol.set_number_of_copies)
                # element 6 - <xs:complexType name="base_macromolecule_type">
                # XSD: <xs:element name="oligomeric_state" type="pos_int_or_string_type" minOccurs="0"/>
                self.check_set(spec_component_in.get_oligomericDetails, mol.set_oligomeric_state)
                # element 7 - <xs:complexType name="base_macromolecule_type">
                # XSD: <xs:element name="recombinant_exp_flag" type="xs:boolean" maxOccurs="1" minOccurs="0">
                if not label:
                    self.check_set(spec_component_in.get_recombinantExpFlag, mol.set_recombinant_exp_flag)

        def set_oddity_details(smol_or_mol_in, comp_in, smol_or_mol):
            """
            """
            all_details = ''
            base_details = smol_or_mol_in.get_details()
            if base_details is not None and base_details != '':
                all_details = base_details
            comp_details = comp_in.get_details()
            if comp_details is not None and comp_details != '':
                # add this to the base details
                all_details += ' {comp details}: ' + comp_details
            if all_details != '':
                smol_or_mol.set_details(all_details)

        sample_in = xml_in.get_sample()
        smol_index = 0
        if sample_in is not None:
            # XSD: <xs:complexType name="sample_type"> has 3 elements
            sample = emdb30.sample_type()
            # element 1 - <xs:complexType name="sample_type">
            # XSD: <xs:element name="name" type="sci_name_type">
            sample.set_name(emdb30.sci_name_type(valueOf_=sample_in.get_name()))
            # element 2 - <xs:complexType name="sample_type">
            # XSD: <xs:element name="supramolecule_list" maxOccurs="1"> has 1 element
            sup_mol_list = emdb30.supramolecule_listType()
            # element 1 - <xs:element name="supramolecule_list" maxOccurs="1">
            # XSD: <xs:element ref="supramolecule" maxOccurs="unbounded"/> head element for 6 substitution groups
            # substitution group 1 - <xs:element ref="supramolecule" maxOccurs="unbounded"/>
            # XSD: <xs:element name="sample_supramolecule" substitutionGroup="supramolecule" type="sample_supramolecule_type">
            # XSD: <xs:complexType name="sample_supramolecule_type"> has a base and 3 elements
            sample_supmol = emdb30.sample_supramolecule_type()
            sample_supmol.original_tagname_ = 'sample_supramolecule'
            # base - <xs:complexType name="sample_supramolecule_type">
            # XSD: <xs:extension base="base_supramolecule_type">
            set_base_supramolecule(sample_supmol, const.EM_SAMPLE_ID, sample_in)
            # element 1 - <xs:complexType name="sample_supramolecule_type">
            # XSD: <xs:element name="natural_source" minOccurs="0" type="sample_natural_source_type" maxOccurs="unbounded"/>
            # element 2 - <xs:complexType name="sample_supramolecule_type">
            # XSD: <xs:element name="number_unique_components" type="xs:positiveInteger" minOccurs="0"/>
            self.check_set(sample_in.get_numComponents, sample_supmol.set_number_unique_components)
            # element 3 - <xs:complexType name="sample_supramolecule_type">
            # XSD: <xs:element name="molecular_weight" type="molecular_weight_type" minOccurs="0"/>
            set_mol_weight(sample_supmol.set_molecular_weight, sample_in.get_molWtTheo(), sample_in.get_molWtExp(), sample_in.get_molWtMethod())

            sup_mol_list.add_supramolecule(sample_supmol)
            # element 3 - <xs:complexType name="sample_type">
            # XSD: <xs:element name="macromolecule_list" type="macromolecule_list_type" minOccurs="0"> has 1 element
            mol_list = emdb30.macromolecule_list_type()
            # element 1 - <xs:element name="macromolecule_list" type="macromolecule_list_type" minOccurs="0">
            # XSD: <xs:element ref="macromolecule" maxOccurs="unbounded"/>
            # XSD: <xs:element name="macromolecule" type="base_macromolecule_type" abstract="true"> head for 7 substitute groups
            comp_in = sample_in.get_sampleComponentList()
            mol_index = 0
            if comp_in is not None:
                comp_list_in = comp_in.get_sampleComponent()
                for component_in in comp_list_in:
                    c_type = component_in.get_entry()
                    if c_type == 'protein':
                        # substitution group 1 - <xs:element ref="macromolecule" maxOccurs="unbounded"/>
                        # XSD: <xs:element name="protein_or_peptide" substitutionGroup="macromolecule" type="protein_or_peptide_macromolecule_type"/> has a base and 4 elements
                        protein_mol = emdb30.protein_or_peptide_macromolecule_type()
                        protein_mol.original_tagname_ = 'protein_or_peptide'
                        # base - <xs:element name="protein_or_peptide" substitutionGroup="macromolecule" type="protein_or_peptide_macromolecule_type"/>
                        # XSD: <xs:extension base="base_macromolecule_type">
                        p_in = component_in.get_protein()
                        mol_index = mol_index + 1
                        set_base_macromolecule(protein_mol, mol_index, component_in, p_in)
                        if p_in is not None:
                            # element 1 - <xs:element name="protein_or_peptide" substitutionGroup="macromolecule" type="protein_or_peptide_macromolecule_type"/>
                            # XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" minOccurs="0"/>
                            eng_source = p_in.get_engSource()
                            copy_recombinant_source(eng_source, protein_mol.set_recombinant_expression)
                            # element 2 - <xs:element name="protein_or_peptide" substitutionGroup="macromolecule" type="protein_or_peptide_macromolecule_type"/>
                            # XSD: <xs:element name="enantiomer">
                            # element 3 - <xs:element name="protein_or_peptide" substitutionGroup="macromolecule" type="protein_or_peptide_macromolecule_type"/>
                            # XSD: <xs:element name="sequence">
                            seq = emdb30.sequenceType()
                            p_ext_refs = p_in.get_externalReferences()
                            if p_ext_refs is not None and p_ext_refs != []:
                                add_mol_references(seq.add_external_references, p_ext_refs)
                            protein_mol.set_sequence(seq)
                            # element 4 - <xs:element name="protein_or_peptide" substitutionGroup="macromolecule" type="protein_or_peptide_macromolecule_type"/>
                            # XSD: <xs:element name="ec_number" maxOccurs="unbounded" minOccurs="0">

                            # AN ODDITY: proteins in 1.9 can have details given at two elements and only one element to write them in
                            # if there are details in p_in they should be added with a flag
                            set_oddity_details(component_in, p_in, protein_mol)

                        mol_list.add_macromolecule(protein_mol)
                    elif c_type == 'ligand':
                        # substitution group 2 - <xs:element ref="macromolecule" maxOccurs="unbounded"/>
                        # XSD: <xs:element name="ligand" substitutionGroup="macromolecule" type="ligand_macromolecule_type"> has a base and 3 elements
                        ligand_mol = emdb30.ligand_macromolecule_type()
                        ligand_mol.original_tagname_ = 'ligand'
                        # base - <xs:element name="ligand" substitutionGroup="macromolecule" type="ligand_macromolecule_type">
                        # XSD: <xs:extension base="base_macromolecule_type">
                        l_in = component_in.get_ligand()
                        mol_index = mol_index + 1
                        set_base_macromolecule(ligand_mol, mol_index, component_in, l_in)
                        if l_in is not None:
                            # element 1 - <xs:element name="ligand" substitutionGroup="macromolecule" type="ligand_macromolecule_type">
                            # XSD: <xs:element name="formula" type="formula_type" minOccurs="0"/>
                            # element 2 - <xs:element name="ligand" substitutionGroup="macromolecule" type="ligand_macromolecule_type">
                            # XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                            l_ext_refs = l_in.get_externalReferences()
                            if l_ext_refs is not None and l_ext_refs.hasContent_():
                                add_mol_references(ligand_mol.add_external_references, l_ext_refs)
#                                 for ext_ref in ext_refs:
#                                     if ext_ref.original_tagname_ == 'refUniProt':
#                                         ligand_mol.add_external_references(valueOf_=ext_ref.valueOf_, type='')
#                                     if ext_ref.original_tagname_ == 'refGo':
#                                         ligand_mol.add_external_references(valueOf_=ext_ref.valueOf_, type='')
#                                     if ext_ref.original_tagname_ == 'refInterpro':
#                                         ligand_mol.add_external_references(valueOf_=ext_ref.valueOf_, type='')
                            # element 3 - <xs:element name="ligand" substitutionGroup="macromolecule" type="ligand_macromolecule_type">
                            # XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" minOccurs="0"/>
                            eng_source = l_in.get_engSource()
                            copy_recombinant_source(eng_source, ligand_mol.set_recombinant_expression)

                            # AN ODDITY: proteins in 1.9 can have details given at two elements and only one element to write them in
                            # if there are details in p_in they should be added with a flag
                            set_oddity_details(component_in, l_in, ligand_mol)

                        mol_list.add_macromolecule(ligand_mol)
                    elif c_type == 'label':
                        # substitution group 3 - <xs:element ref="macromolecule" maxOccurs="unbounded"/>
                        # XSD: <xs:element name="em_label" substitutionGroup="macromolecule" type="em_label_macromolecule_type"/> has a base and 1 element
                        em_label_mol = emdb30.em_label_macromolecule_type()
                        em_label_mol.original_tagname_ = 'em_label'
                        # base - <xs:element name="em_label" substitutionGroup="macromolecule" type="em_label_macromolecule_type"/>
                        # XSD: <xs:extension base="base_macromolecule_type">
                        l_in = component_in.get_label()
                        mol_index = mol_index + 1
                        set_base_macromolecule(em_label_mol, mol_index, component_in, l_in, label=True)
                        if l_in is not None:
                            # element 1 - <xs:element name="em_label" substitutionGroup="macromolecule" type="em_label_macromolecule_type"/>
                            # XSD: <xs:element name="formula" type="formula_type" minOccurs="0"/>
                            self.check_set(l_in.get_formula, em_label_mol.set_formula)

                        mol_list.add_macromolecule(em_label_mol)
                    elif c_type == 'nucleic-acid':
                        na_in = component_in.get_nucleic_acid()
                        if na_in is not None:
                            na_class_in = na_in.get_class()
                            if na_class_in == 'DNA':
                                # substitution group 4 - <xs:element ref="macromolecule" maxOccurs="unbounded"/>
                                # XSD: <xs:element name="dna" substitutionGroup="macromolecule" type="dna_macromolecule_type"> has a base and 4 elements
                                dna_mol = emdb30.dna_macromolecule_type()
                                dna_mol.original_tagname_ = 'dna'
                                # base - <xs:element name="dna" substitutionGroup="macromolecule" type="dna_macromolecule_type">
                                # XSD: <xs:extension base="base_macromolecule_type">
                                mol_index = mol_index + 1
                                set_base_macromolecule(dna_mol, mol_index, component_in, na_in, nucleic_acid=True)
                                # element 1 - <xs:element name="dna" substitutionGroup="macromolecule" type="dna_macromolecule_type">
                                # XSD: <xs:element name="sequence"> has 3 elements
                                seq = emdb30.sequenceType()
                                # element 1 - <xs:element name="sequence">
                                # XSD: <xs:element name="string">
                                seq_in = na_in.get_sequence()
                                if seq_in is not None:
                                    seq.set_string(seq_in)
                                # element 2 - <xs:element name="sequence">
                                # XSD: <xs:element name="discrepancy_list" minOccurs="0">
                                # element 3 - <xs:element name="sequence">
                                # XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                                if seq.hasContent_():
                                    dna_mol.set_sequence(seq)
                                # element 2 - <xs:element name="dna" substitutionGroup="macromolecule" type="dna_macromolecule_type">
                                # XSD: <xs:element name="classification" minOccurs="0">
                                dna_mol.set_classification('DNA')
                                # element 3 - <xs:element name="dna" substitutionGroup="macromolecule" type="dna_macromolecule_type">
                                # XSD: <xs:element name="structure" type="xs:token" minOccurs="0">
                                self.check_set(na_in.get_structure, dna_mol.set_structure)
                                # element 4 - <xs:element name="dna" substitutionGroup="macromolecule" type="dna_macromolecule_type">
                                # XSD: <xs:element name="synthetic_flag" type="xs:boolean" minOccurs="0">
                                self.check_set(na_in.get_syntheticFlag, dna_mol.set_synthetic_flag)

                                mol_list.add_macromolecule(dna_mol)
                            elif na_class_in == 'RNA' or na_class_in == 'T-RNA':
                                # substitution group 5 - <xs:element ref="macromolecule" maxOccurs="unbounded"/>
                                # XSD: <xs:element name="rna" substitutionGroup="macromolecule" type="rna_macromolecule_type"> has a base and 5 elements
                                rna_mol = emdb30.rna_macromolecule_type()
                                rna_mol.original_tagname_ = 'rna'
                                # base - <xs:element name="rna" substitutionGroup="macromolecule" type="rna_macromolecule_type">
                                # XSD: <xs:extension base="base_macromolecule_type">
                                mol_index = mol_index + 1
                                set_base_macromolecule(rna_mol, mol_index, component_in, na_in, nucleic_acid=True)
                                # element 1 - <xs:element name="rna" substitutionGroup="macromolecule" type="rna_macromolecule_type">
                                # XSD: <xs:element name="sequence"> has 3 elements
                                seq = emdb30.sequenceType()
                                seq_in = na_in.get_sequence()
                                if seq_in is not None:
                                    # element 1 - <xs:element name="sequence">
                                    # XSD: <xs:element name="string">
                                    seq.set_string(seq_in)
                                    # element 2 - <xs:element name="sequence">
                                    # XSD: <xs:element name="discrepancy_list" minOccurs="0">
                                    # element 3 - <xs:element name="sequence">
                                    # XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                                if seq.hasContent_():
                                    rna_mol.set_sequence(seq)
                                # element 2 - <xs:element name="rna" substitutionGroup="macromolecule" type="rna_macromolecule_type">
                                # XSD: <xs:element name="classification" minOccurs="0">
                                na_class = 'OTHER'
                                if na_class_in == 'T-RNA':
                                    na_class = 'TRANSFER'
                                rna_mol.set_classification(na_class)
                                # element 3 - <xs:element name="rna" substitutionGroup="macromolecule" type="rna_macromolecule_type">
                                # XSD: <xs:element name="structure" type="xs:token" minOccurs="0">
                                self.check_set(na_in.get_structure, rna_mol.set_structure)
                                # element 4 - <xs:element name="rna" substitutionGroup="macromolecule" type="rna_macromolecule_type">
                                # XSD: <xs:element name="synthetic_flag" type="xs:boolean" minOccurs="0">
                                self.check_set(na_in.get_syntheticFlag, rna_mol.set_synthetic_flag)
                                # element 5 - <xs:element name="rna" substitutionGroup="macromolecule" type="rna_macromolecule_type">
                                # XSD: <xs:element name="ec_number" maxOccurs="unbounded" minOccurs="0">
                                if rna_mol.hasContent_():
                                    mol_list.add_macromolecule(rna_mol)
                            elif na_class_in in ['DNA/RNA', 'OTHER']:
                                # substitution group 6 - <xs:element ref="macromolecule" maxOccurs="unbounded"/>
                                # XSD: <xs:element name="other_macromolecule" substitutionGroup="macromolecule" type="other_macromolecule_type"> has a base and 5 elements
                                other_mol = emdb30.other_macromolecule_type()
                                other_mol.original_tagname_ = 'other_macromolecule'
                                # base - <xs:element name="other_macromolecule" substitutionGroup="macromolecule" type="other_macromolecule_type">
                                # XSD: <xs:extension base="base_macromolecule_type">
                                mol_index = mol_index + 1
                                set_base_macromolecule(other_mol, mol_index, component_in, na_in, nucleic_acid=True)
                                # element 1 - <xs:element name="other_macromolecule" substitutionGroup="macromolecule" type="other_macromolecule_type">
                                # XSD: <xs:element name="sequence" minOccurs="0"> has 3 elements
                                seq = emdb30.sequenceType()
                                seq_in = na_in.get_sequence()
                                if seq_in is not None:
                                    # element 1 - <xs:element name="sequence" minOccurs="0">
                                    # XSD: <xs:element name="string">
                                    seq.set_string(seq_in)
                                    # element 2 - <xs:element name="sequence" minOccurs="0">
                                    # XSD: <xs:element name="discrepancy_list" minOccurs="0">
                                    # element 3 - <xs:element name="sequence" minOccurs="0">
                                    # XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                                if seq.hasContent_():
                                    other_mol.set_sequence(seq)
                                # element 2 - <xs:element name="other_macromolecule" substitutionGroup="macromolecule" type="other_macromolecule_type">
                                # XSD: <xs:element name="classification" type="xs:token">
                                if na_class_in == 'OTHER':
                                    other_mol.set_classification('OTHER_NA')
                                else:
                                    other_mol.set_classification(na_class_in)
                                # element 3 - <xs:element name="other_macromolecule" substitutionGroup="macromolecule" type="other_macromolecule_type">
                                # XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" minOccurs="0"/>
                                # element 4 - <xs:element name="other_macromolecule" substitutionGroup="macromolecule" type="other_macromolecule_type">
                                # XSD: <xs:element name="structure" type="xs:token" minOccurs="0">
                                self.check_set(na_in.get_structure, other_mol.set_structure)
                                # element 5 - <xs:element name="other_macromolecule" substitutionGroup="macromolecule" type="other_macromolecule_type">
                                # XSD: <xs:element name="synthetic_flag" type="xs:boolean" minOccurs="0">
                                self.check_set(na_in.get_syntheticFlag, other_mol.set_synthetic_flag)

                                if other_mol.hasContent_():
                                    mol_list.add_macromolecule(other_mol)
                    elif c_type == 'virus':
                        # substitution group 2 - <xs:element ref="supramolecule" maxOccurs="unbounded"/>
                        # XSD: <xs:element name="virus_supramolecule" substitutionGroup="supramolecule" type="virus_supramolecule_type"/>
                        virus_smol = emdb30.virus_supramolecule_type()
                        virus_smol.original_tagname_ = 'virus_supramolecule'
                        # XSD: <xs:complexType name="virus_supramolecule_type"> has a base and 14 elements
                        virus_in = component_in.get_virus()
                        if virus_in is not None:
                            # base - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:extension base="base_supramolecule_type">
                            smol_index = smol_index + 1
                            set_base_supramolecule(virus_smol, smol_index, virus_in, component_in)
                            # element 1 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="sci_species_name" type="virus_species_name_type" minOccurs="0"/>
                            sci_species_name = virus_in.get_sciSpeciesName()
                            if sci_species_name is not None:
                                virus_smol.set_sci_species_name(emdb30.virus_species_name_type(valueOf_=sci_species_name.get_valueOf_(), ncbi=sci_species_name.get_ncbiTaxId()))
                            # element 2 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="sci_species_strain" type="xs:string" maxOccurs="1" minOccurs="0"/>
                            self.check_set(virus_in.get_sciSpeciesStrain, virus_smol.set_sci_species_strain)
                            # element 3 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="natural_host" type="virus_natural_host_type" minOccurs="0" maxOccurs="unbounded"/>
                            all_vir_ns_in = virus_in.get_natSource()
                            if all_vir_ns_in is not None:
                                for vir_ns_in in all_vir_ns_in:
                                    vir_nat_source = emdb30.virus_natural_host_type()
                                    # XSD: <xs:complexType name="virus_natural_host_type">, xs:extension base="base_source_type"/> has 3 elements and 1 attribute
                                    # attribute 1 - <xs:complexType name="base_source_type">
                                    # XSD: <xs:attribute name="database" use="required">
                                    vir_nat_source.set_database('NCBI')
                                    # element 1 - <xs:complexType name="base_source_type">
                                    # XSD: <xs:element name="organism" type="organism_type">
                                    org = emdb30.organism_type()
                                    # XSD: <xs:complexType name="organism_type"> has 1 attribute and is ext of token
                                    virus_host_species = vir_ns_in.get_hostSpecies()
                                    if virus_host_species is not None:
                                        virus_hst_spec = virus_host_species.valueOf_
                                        org.set_valueOf_(virus_hst_spec)
                                        # attribute 1 - <xs:complexType name="organism_type">
                                        # XSD: <xs:attribute name="ncbi" type="xs:positiveInteger"/>
                                        self.check_set(virus_host_species.get_ncbiTaxId, org.set_ncbi)
                                    vir_nat_source.set_organism(org)
                                    # element 2 - <xs:complexType name="base_source_type">
                                    # XSD: <xs:element name="strain" type="organism_type" minOccurs="0"/>
                                    strain_in = vir_ns_in.get_hostSpeciesStrain()
                                    if strain_in is not None:
                                        strain = emdb30.organism_type()
                                        # XSD: <xs:complexType name="organism_type"> has 1 attribute and is ext of token
                                        strain.set_valueOf_(strain_in)
                                        # attribute 1 - <xs:complexType name="organism_type">
                                        # XSD: <xs:attribute name="ncbi" type="xs:positiveInteger"/>
                                        vir_nat_source.set_strain(strain)
                                    # element 3 - <xs:complexType name="base_source_type">
                                    # XSD: <xs:element name="synonym_organism" type="xs:token" minOccurs="0">
                                    self.check_set(vir_ns_in.get_hostCategory, vir_nat_source.set_synonym_organism)
                                    virus_smol.add_natural_host(vir_nat_source)
                            # element 4 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="host_system" type="recombinant_source_type" minOccurs="0"/>
                            es_in = virus_in.get_engSource()
                            if es_in is not None and es_in != []:
                                e_value = es_in[0]
                                if e_value.hasContent_:
                                    # XSD: <xs:complexType name="recombinant_source_type"> has 1 attribute and 5 elements
                                    rec_source = emdb30.recombinant_source_type()
                                    # attribute 1 - <xs:complexType name="recombinant_source_type">
                                    # XSD: <xs:attribute name="database" use="required">
                                    rec_source.set_database('NCBI')
                                    # element 1 - <xs:complexType name="recombinant_source_type">
                                    # XSD: <xs:element name="organism" type="organism_type">
                                    exp_sys_in = e_value.get_expSystem()
                                    if exp_sys_in is not None:
                                        # XSD: <xs:complexType name="organism_type"> is a token and has 1 attribute
                                        org = emdb30.organism_type()
                                        self.check_set(exp_sys_in.get_valueOf_, org.set_valueOf_)
                                        # attribute 1 - <xs:complexType name="organism_type">
                                        # XSD: <xs:attribute name="ncbi" type="xs:positiveInteger"/>
                                        self.check_set(exp_sys_in.get_ncbiTaxId, org.set_ncbi)
                                        rec_source.set_recombinant_organism(org)
                                    # element 2 - <xs:complexType name="recombinant_source_type">
                                    # XSD: <xs:element name="strain" type="xs:token" minOccurs="0"/>
                                    self.check_set(e_value.get_expSystemStrain, rec_source.set_recombinant_strain)
                                    # element 3 - <xs:complexType name="recombinant_source_type">
                                    # XSD: <xs:element name="cell" type="xs:token" minOccurs="0">
                                    self.check_set(e_value.get_expSystemCell, rec_source.set_recombinant_cell)
                                    # element 4 - <xs:complexType name="recombinant_source_type">
                                    # XSD: <xs:element name="plasmid" type="xs:token" minOccurs="0"/>
                                    self.check_set(e_value.get_vector, rec_source.set_recombinant_plasmid)
                                    # element 5 - <xs:complexType name="recombinant_source_type">
                                    # XSD: <xs:element name="synonym_organism" type="xs:token" minOccurs="0">
                                    # if rec_source.hasContent_():
                                    virus_smol.set_host_system(rec_source)
                            # element 5 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="molecular_weight" type="molecular_weight_type" minOccurs="0"/>
                            set_mol_weight(virus_smol.set_molecular_weight, component_in.get_molWtTheo(), component_in.get_molWtExp())
                            # element 6 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="virus_shell" maxOccurs="unbounded" minOccurs="0"> has 1 attribute and 3 elements
                            shell_list_in = virus_in.get_shell()
                            for shell_in in shell_list_in:
                                shell = emdb30.virus_shellType()
                                # attribute 1 - <xs:element name="virus_shell" maxOccurs="unbounded" minOccurs="0">
                                # XSD: <xs:attribute name="id" type="xs:positiveInteger"/>
                                self.check_set(shell_in.get_id, shell.set_shell_id)
                                # element 1 - <xs:element name="virus_shell" maxOccurs="unbounded" minOccurs="0">
                                # XSD: <xs:element name="name" type="xs:token" nillable="false" minOccurs="0"/>
                                self.check_set(shell_in.get_nameElement, shell.set_name)
                                # element 2 - <xs:element name="virus_shell" maxOccurs="unbounded" minOccurs="0">
                                # XSD: <xs:element name="diameter" minOccurs="0">
                                shell_diam = shell_in.get_diameter()
                                if shell_diam is not None:
                                    shell.set_diameter(emdb30.diameterType(valueOf_=shell_diam.valueOf_, units=const.U_ANG))
                                # element 3 - <xs:element name="virus_shell" maxOccurs="unbounded" minOccurs="0">
                                # XSD: <xs:element name="triangulation" type="xs:positiveInteger" minOccurs="0"/>
                                self.check_set(shell_in.get_tNumber, shell.set_triangulation, int)
                                virus_smol.add_virus_shell(shell)
                            # element 7 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="virus_type">
                            vir_type = virus_in.get_class()
                            if vir_type is not None:
                                virus_smol.set_virus_type(vir_type)
                            else:
                                if not self.roundtrip:
                                    # virus_type  and class are mandatory
                                    virus_smol.set_virus_type('OTHER')
                            # element 8 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="virus_isolate">
                            virus_smol.set_virus_isolate(virus_in.get_isolate())
                            # element 9 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="virus_enveloped" type="xs:boolean"/>
                            virus_smol.set_virus_enveloped(virus_in.get_enveloped())
                            # element 10 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="virus_empty" type="xs:boolean"/>
                            virus_smol.set_virus_empty(virus_in.get_empty())
                            # element 11 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="syn_species_name" type="xs:string" maxOccurs="1" minOccurs="0">
                            self.check_set(virus_in.get_synSpeciesName, virus_smol.set_syn_species_name)
                            # element 12 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="sci_species_serotype" type="xs:string" maxOccurs="1" minOccurs="0">
                            self.check_set(virus_in.get_sciSpeciesSerotype, virus_smol.set_sci_species_serotype)
                            # element 13 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="sci_species_serocomplex" type="xs:string" maxOccurs="1" minOccurs="0">
                            self.check_set(virus_in.get_sciSpeciesSerocomplex, virus_smol.set_sci_species_serocomplex)
                            # element 14 - <xs:complexType name="virus_supramolecule_type">
                            # XSD: <xs:element name="sci_species_subspecies" type="xs:string" maxOccurs="1" minOccurs="0">
                            self.check_set(virus_in.get_sciSpeciesSubspecies, virus_smol.set_sci_species_subspecies)

                        sup_mol_list.add_supramolecule(virus_smol)
                    elif c_type == 'cellular-component':
                        # substitution group 3 - <xs:element ref="supramolecule" maxOccurs="unbounded"/>
                        # XSD: <xs:element name="organelle_or_cellular_component_supramolecule" substitutionGroup="supramolecule" type="organelle_or_cellular_component_supramolecule_type"/>
                        comp_smol = emdb30.organelle_or_cellular_component_supramolecule_type()
                        comp_smol.original_tagname_ = 'organelle_or_cellular_component_supramolecule'
                        cell_comp_in = component_in.get_cellular_component()
                        if cell_comp_in is not None:
                            # base - <xs:element name="organelle_or_cellular_component_supramolecule" substitutionGroup="supramolecule" type="organelle_or_cellular_component_supramolecule_type"/>
                            # has a base and 3 elements
                            # XSD: <xs:extension base="base_macromolecule_type">
                            smol_index = smol_index + 1
                            set_base_supramolecule(comp_smol, smol_index, cell_comp_in, component_in)
                            # element 1 - <xs:element name="organelle_or_cellular_component_supramolecule" substitutionGroup="supramolecule" type="organelle_or_cellular_component_supramolecule_type"/>
                            # XSD: <xs:element name="natural_source" minOccurs="0" type="organelle_natural_source_type" maxOccurs="unbounded"/>
                            # XSD: <xs:complexType name="organelle_natural_source_type"> has base and 5 elements
                            org_nat_source = emdb30.organelle_natural_source_type()
                            set_mol_natural_source(org_nat_source, component_in, cell_comp_in)#, organelle=False)
                            if org_nat_source.hasContent_():
                                comp_smol.add_natural_source(org_nat_source)
                        # element 2 - <xs:element name="organelle_or_cellular_component_supramolecule" substitutionGroup="supramolecule" type="organelle_or_cellular_component_supramolecule_type"/>
                        # XSD: <xs:element name="molecular_weight" type="molecular_weight_type" minOccurs="0"/>
                        set_mol_weight(comp_smol.set_molecular_weight, component_in.get_molWtTheo(), component_in.get_molWtExp())
                        # element 3 - <xs:element name="organelle_or_cellular_component_supramolecule" substitutionGroup="supramolecule" type="organelle_or_cellular_component_supramolecule_type"/>
                        # XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" maxOccurs="1" minOccurs="0">
                        if cell_comp_in is not None:
                            eng_src = cell_comp_in.get_engSource()
                            copy_recombinant_source(eng_src, comp_smol.set_recombinant_expression)

                        # AN ODDITY: proteins in 1.9 can have details given at two elements and only one element to write them in
                        # if there are details in p_in they should be added with a flag
                        set_oddity_details(component_in, cell_comp_in, comp_smol)

                        sup_mol_list.add_supramolecule(comp_smol)
                    elif c_type in ['ribosome-eukaryote', 'ribosome-prokaryote']:
                        # substitution group 4 - <xs:element ref="supramolecule" maxOccurs="unbounded"/>
                        # XSD: <xs:element name="complex_supramolecule" substitutionGroup="supramolecule" type="complex_supramolecule_type"/>
                        complex_smol = emdb30.complex_supramolecule_type()
                        complex_smol.original_tagname_ = 'complex_supramolecule'
                        complex_smol_in = None
                        if c_type == 'ribosome-eukaryote':
                            complex_smol_in = component_in.get_ribosome_eukaryote()
                        elif c_type == 'ribosome-prokaryote':
                            complex_smol_in = component_in.get_ribosome_prokaryote()
                        if complex_smol_in is not None:
                            # XSD: <xs:complexType name="complex_supramolecule_type"> has a base and 4 elements and 1 attribute
                            # base - <xs:complexType name="complex_supramolecule_type">
                            # XSD: <xs:extension base="base_supramolecule_type">
                            smol_index = smol_index + 1
                            set_base_supramolecule(complex_smol, smol_index, complex_smol_in, component_in, rib_cat=c_type)
                            # attribute 1 - <xs:complexType name="complex_supramolecule_type">
                            # XSD: <xs:attribute name="chimera" type="xs:boolean" fixed="true"/>
                            # element 1 - <xs:complexType name="complex_supramolecule_type">
                            # XSD: <xs:element name="natural_source" type="complex_natural_source_type" minOccurs="0" maxOccurs="unbounded"/>
                            complex_smol_nat_source = emdb30.complex_natural_source_type()
                            set_mol_natural_source(complex_smol_nat_source, component_in, complex_smol_in, tissue=False, cell=False, organelle=False, cell_loc=False)
                            if complex_smol_nat_source.hasContent_():
                                complex_smol.add_natural_source(complex_smol_nat_source)
                            # element 2 - <xs:complexType name="complex_supramolecule_type">
                            # XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" maxOccurs="unbounded" minOccurs="0">
                            copy_recombinant_source(complex_smol_in.get_engSource(), complex_smol.add_recombinant_expression)
                            # element 3 - <xs:complexType name="complex_supramolecule_type">
                            # XSD: <xs:element name="molecular_weight" type="molecular_weight_type" minOccurs="0"/>
                            set_mol_weight(complex_smol.set_molecular_weight, component_in.get_molWtTheo(), component_in.get_molWtExp())
                            # element 4 - <xs:complexType name="complex_supramolecule_type">
                            # XSD: <xs:element name="ribosome-details" type="xs:string" minOccurs="0">
                            add_datails = c_type + ': '
                            if c_type == 'ribosome-eukaryote':
                                complex_smol.set_ribosome_details(add_datails + complex_smol_in.get_eukaryote())
                            elif c_type == 'ribosome-prokaryote':
                                complex_smol.set_ribosome_details(add_datails + complex_smol_in.get_prokaryote())
                            # AN ODDITY: proteins in 1.9 can have details given at two elements and only one element to write them in
                            # if there are details in p_in they should be added with a flag
                            set_oddity_details(component_in, complex_smol_in, complex_smol)

                        sup_mol_list.add_supramolecule(complex_smol)

            sample.set_supramolecule_list(sup_mol_list)

            if mol_list.hasContent_():
                sample.set_macromolecule_list(mol_list)

            xml_out.set_sample(sample)

        # element 4 - <xs:complexType name="entry_type">
        # XSD: <xs:element name="structure_determination_list"> has 1 element
        # element 1 - <xs:element name="structure_determination_list">
        # XSD: <xs:element name="structure_determination" type="structure_determination_type" maxOccurs="unbounded"/>
        sd_list = emdb30.structure_determination_listType()
        # XSD: <xs:complexType name="structure_determination_type"> has 6 elements and 1 attribute
        struct_det = emdb30.structure_determination_type()
        # attribute 1 - <xs:complexType name="structure_determination_type">
        # XSD: <xs:attribute name="id" type="xs:positiveInteger" use="required">
        # only ever 1 processing element - therefore assume 1
        struct_det.set_structure_determination_id(1)
        # element 1 - <xs:complexType name="structure_determination_type">
        # XSD: <xs:element name="method">
        process_in = xml_in.get_processing()
        if process_in is not None:
            em_method = process_in.get_method()
        else:
            # assume single particle
            em_method = const.EMM_SP
        if em_method in [const.EMM_SP, const.EMM_STOM, const.EMM_TOM]:
            struct_det.set_method(em_method)
        elif em_method == 'twoDCrystal':
            struct_det.set_method(const.EMM_EC)
        elif em_method == const.EMM_HEL:
            struct_det.set_method(const.EMM_HEL)

        exp_in = xml_in.get_experiment()
        spec_prep_in = exp_in.get_specimenPreparation()
        # element 2 - <xs:complexType name="structure_determination_type">
        # XSD: <xs:element name="aggregation_state">
        if spec_prep_in is not None:
            self.check_set(spec_prep_in.get_specimenState, struct_det.set_aggregation_state)
        # element 3 - <xs:complexType name="structure_determination_type">
        # XSD: <xs:element name="macromolecules_and_complexes" type="macromolecules_and_complexes_type" minOccurs="0">
        # element 4 - <xs:complexType name="structure_determination_type">
        # XSD: <xs:element name="specimen_preparation_list"> has 1 element
        spec_prep_list = emdb30.specimen_preparation_listType()
        # element 1 - <xs:element name="specimen_preparation_list">
        # XSD: <xs:element ref="specimen_preparation" maxOccurs="unbounded"/>
        # XSD: <xs:element name="specimen_preparation" type="base_preparation_type" abstract="true"> base for 5 substitution groups
        vitr_in = exp_in.get_vitrification()
        n_sp = max(1, len(vitr_in))
        j = 1
        for i in range(0, n_sp):
            prep = None
            if em_method == const.EMM_TOM:
                # substitution group 1 - <xs:element name="specimen_preparation" type="base_preparation_type" abstract="true">
                # XSD: <xs:element name="tomography_preparation" type="tomography_preparation_type" substitutionGroup="specimen_preparation">
                prep = emdb30.tomography_preparation_type()
                prep.original_tagname_ = 'tomography_preparation'
                # XSD: <xs:complexType name="tomography_preparation_type"> has a base and 5 elements
                # base - <xs:complexType name="tomography_preparation_type">
                set_base_preparation(prep, spec_prep_in, vitr_in)
                # XSD: <xs:extension base="base_preparation_type">
                # element 1 - <xs:complexType name="tomography_preparation_type">
                # XSD: <xs:element name="fiducial_markers_list" minOccurs="0">
                # element 2 - <xs:complexType name="tomography_preparation_type">
                # XSD: <xs:element name="high_pressure_freezing" minOccurs="0">
                # element 3 - <xs:complexType name="tomography_preparation_type">
                # XSD: <xs:element name="embedding_material" type="xs:token" minOccurs="0">
                # element 4 - <xs:complexType name="tomography_preparation_type">
                # XSD: <xs:element name="cryo_protectant" type="xs:token" minOccurs="0">
                # element 5 - <xs:complexType name="tomography_preparation_type">
                # XSD: <xs:element name="sectioning" minOccurs="0">
            elif em_method == const.EMM_SP:
                # substitution group 2 - <xs:element name="specimen_preparation" type="base_preparation_type" abstract="true">
                # XSD: <xs:element name="single_particle_preparation" type="single_particle_preparation_type" substitutionGroup="specimen_preparation"/>
                prep = emdb30.single_particle_preparation_type()
                prep.original_tagname_ = 'single_particle_preparation'
                # XSD: <xs:complexType name="single_particle_preparation_type"> has a base
                # base - <xs:complexType name="single_particle_preparation_type">
                set_base_preparation(prep, spec_prep_in, vitr_in)
            elif em_method == const.EMM_STOM:
                # substitution group 3 - <xs:element name="specimen_preparation" type="base_preparation_type" abstract="true">
                # XSD: <xs:element name="subtomogram_averaging_preparation" type="subtomogram_averaging_preparation_type" substitutionGroup="specimen_preparation"/>
                prep = emdb30.subtomogram_averaging_preparation_type()
                prep.original_tagname_ = 'subtomogram_averaging_preparation'
                # XSD: <xs:complexType name="subtomogram_averaging_preparation_type"> has a base
                # base - <xs:complexType name="subtomogram_averaging_preparation_type">
                set_base_preparation(prep, spec_prep_in, vitr_in)
            elif em_method == const.EMM_HEL:
                # substitution group 4 - <xs:element name="specimen_preparation" type="base_preparation_type" abstract="true">
                # XSD: <xs:element name="helical_preparation" type="helical_preparation_type" substitutionGroup="specimen_preparation"/>
                prep = emdb30.helical_preparation_type()
                prep.original_tagname_ = 'helical_preparation'
                # XSD: <xs:complexType name="helical_preparation_type"> has a base
                # base - <xs:complexType name="helical_preparation_type">
                set_base_preparation(prep, spec_prep_in, vitr_in)
            elif em_method == 'twoDCrystal':
                # substitution group 5 - <xs:element name="specimen_preparation" type="base_preparation_type" abstract="true">
                # XSD: <xs:element name="crystallography_preparation" type="crystallography_preparation_type" substitutionGroup="specimen_preparation">
                prep = emdb30.crystallography_preparation_type()
                prep.original_tagname_ = 'crystallography_preparation'
                # XSD: <xs:complexType name="crystallography_preparation_type"> has a base and 1 element
                # base - <xs:complexType name="crystallography_preparation_type">
                set_base_preparation(prep, spec_prep_in, vitr_in)
                # element 1 - <xs:complexType name="crystallography_preparation_type">
                # XSD: <xs:element name="crystal_formation"> has 7 elements
                x_form = emdb30.crystal_formationType()
                # element 1 - <xs:element name="crystal_formation">
                # XSD: <xs:element name="lipid_protein_ratio" type="xs:float" minOccurs="0"/>
                # element 2 - <xs:element name="crystal_formation">
                # XSD: <xs:element name="lipid_mixture" type="xs:token" minOccurs="0"/>
                # element 3 - <xs:element name="crystal_formation">
                # XSD: <xs:element name="instrument" minOccurs="0">
                # element 4 - <xs:element name="crystal_formation">
                # XSD: <xs:element name="atmosphere" type="xs:token" minOccurs="0"/>
                # element 5 - <xs:element name="crystal_formation">
                # XSD: <xs:element name="temperature" type="crystal_formation_temperature_type" minOccurs="0"/>
                # element 6 - <xs:element name="crystal_formation">
                # XSD: <xs:element name="time" type="crystal_formation_time_type" minOccurs="0"/>
                # element 7 - <xs:element name="crystal_formation">
                # XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                self.check_set(spec_prep_in.get_crystalGrowDetails, x_form.set_details)
                prep.set_crystal_formation(x_form)

            if prep is not None:
                spec_prep_list.add_specimen_preparation(prep)
                prep.set_preparation_id(j)
            j = j + 1
        if spec_prep_list.hasContent_():
            struct_det.set_specimen_preparation_list(spec_prep_list)
        # element 5 - <xs:complexType name="structure_determination_type">
        # XSD: <xs:element name="microscopy_list"> has 1 element
        # element 1 - <xs:element name="microscopy_list">
        microscopy_list = emdb30.microscopy_listType()
        # XSD: <xs:element ref="microscopy" maxOccurs="unbounded"/>
        # XSD: <xs:element name="microscopy" type="base_microscopy_type" abstract="true"> for 5 substitution groups
        im_ac_in = exp_in.get_imageAcquisition()
        imaging_list_in = exp_in.get_imaging()
        # forward reference that will be used in tomography processing
        i = 1
        for img in imaging_list_in:
            if em_method == const.EMM_SP:
                # substitution group 1 - <xs:element name="microscopy" type="base_microscopy_type" abstract="true">
                # XSD: <xs:element name="single_particle_microscopy" type="single_particle_microscopy_type" substitutionGroup="microscopy"/>
                mic = emdb30.single_particle_microscopy_type()
                mic.original_tagname_ = 'single_particle_microscopy'
                # base - <xs:element name="microscopy" type="base_microscopy_type" abstract="true">
                # XSD: <xs:complexType name="single_particle_microscopy_type">
                # XSD: <xs:extension base="base_microscopy_type">
                set_base_microscopy(mic, img, im_ac_in)
            elif em_method == const.EMM_HEL:
                # substitution group 2 - <xs:element name="microscopy" type="base_microscopy_type" abstract="true">
                # XSD: <xs:element name="helical_microscopy" type="helical_microscopy_type" substitutionGroup="microscopy"/> has a base
                mic = emdb30.helical_microscopy_type()
                mic.original_tagname_ = 'helical_microscopy'
                # base - <xs:element name="microscopy" type="base_microscopy_type" abstract="true">
                # XSD: <xs:complexType name="helical_microscopy_type">
                # XSD: <xs:extension base="base_microscopy_type">
                set_base_microscopy(mic, img, im_ac_in)
            elif em_method == const.EMM_TOM:
                # substitution group 3 - <xs:element name="microscopy" type="base_microscopy_type" abstract="true">
                # XSD: <xs:element name="tomography_microscopy" type="tomography_microscopy_type" substitutionGroup="microscopy">
                mic = emdb30.tomography_microscopy_type()
                mic.original_tagname_ = 'tomography_microscopy'
                # base - <xs:element name="microscopy" type="base_microscopy_type" abstract="true">
                # XSD: <xs:extension base="base_microscopy_type">
                set_base_microscopy(mic, img, im_ac_in)
                # element 1 -  <xs:element name="tomography_microscopy" type="tomography_microscopy_type" substitutionGroup="microscopy">
                # XSD: <xs:element name="tilt_series" type="tilt_series_type" maxOccurs="unbounded" minOccurs="0">
                process_in = xml_in.get_processing()
                if process_in is not None:
                    tom_proc = process_in.get_tomography()
                set_tilt_series(mic, img, tom_proc)
            elif em_method == const.EMM_STOM:
                # substitution group 4 - <xs:element name="microscopy" type="base_microscopy_type" abstract="true">
                # XSD: <xs:element name="subtomogram_averaging_microscopy" type="tomography_microscopy_type" substitutionGroup="microscopy"/>
                mic = emdb30.tomography_microscopy_type()
                mic.original_tagname_ = 'subtomogram_averaging_microscopy'
                # base - <xs:element name="microscopy" type="base_microscopy_type" abstract="true">
                # XSD: <xs:extension base="base_microscopy_type">
                set_base_microscopy(mic, img, im_ac_in)
                # element 1 -  <xs:element name="tomography_microscopy" type="tomography_microscopy_type" substitutionGroup="microscopy">
                # XSD: <xs:element name="tilt_series" type="tilt_series_type" maxOccurs="unbounded" minOccurs="0">
                set_tilt_series(mic, img)
            elif em_method == 'twoDCrystal':
                # substitution group 5 - <xs:element name="microscopy" type="base_microscopy_type" abstract="true">
                # XSD: <xs:element name="crystallography_microscopy" type="crystallography_microscopy_type" substitutionGroup="microscopy"> has base and 2 elements one of which is a choice of 2 elements
                mic = emdb30.crystallography_microscopy_type()
                mic.original_tagname_ = 'crystallography_microscopy'
                # base - <xs:element name="microscopy" type="base_microscopy_type" abstract="true">
                # XSD: <xs:extension base="base_microscopy_type">
                set_base_microscopy(mic, img, im_ac_in)
                # element 1 - <xs:element name="crystallography_microscopy" type="crystallography_microscopy_type" substitutionGroup="microscopy">
                # XSD: <xs:element name="camera_length">
                # element 2 - <xs:element name="crystallography_microscopy" type="crystallography_microscopy_type" substitutionGroup="microscopy">
                # 2 choices
                # element 2 - choice 1
                # XSD: <xs:element name="tilt_list" minOccurs="0">
                # element 2 - choice 2
                # XSD: <xs:element name="tilt_series" type="tilt_series_type" maxOccurs="unbounded" minOccurs="0">
                set_tilt_series(mic, img)

            microscopy_list.add_microscopy(mic)
            i += 1

        struct_det.set_microscopy_list(microscopy_list)
        # element 6 - <xs:complexType name="structure_determination_type">
        # XSD: <xs:element ref="image_processing" maxOccurs="unbounded">
        # XSD: <xs:element name="image_processing" type="base_image_processing_type" abstract="true"/> has 5 substitution groups
        # In 1.9 reconstruction is a list and independent of method
        # In 3.0 each reconstruction is mapped to an additional image_processing element
        if process_in is not None:

            def set_ctf_correction(reconstruction, im_proc):
                """
                Method that sets ctf correction - only details are set from v1.9 to v3.0
                """
                ctf_in = reconstruction.get_ctfCorrection()
                if ctf_in is not None:
                    ctf = emdb30.ctf_correction_type()
                    # XSD: <xs:complexType name="ctf_correction_type"> has 5 elements
                    # element 1 - <xs:complexType name="ctf_correction_type">
                    # XSD: <xs:element name="phase_reversal" minOccurs="0">
                    # element 2 - <xs:complexType name="ctf_correction_type">
                    # XSD: <xs:element name="amplitude_correction" minOccurs="0">
                    # element 3 - <xs:complexType name="ctf_correction_type">
                    # XSD: <xs:element name="correction_operation" minOccurs="0">
                    # element 4 - <xs:complexType name="ctf_correction_type">
                    # XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                    # element 5 - <xs:complexType name="ctf_correction_type">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    ctf.set_details(ctf_in)
                    im_proc.set_ctf_correction(ctf)

            def set_crystal_parameters(spec_prep_in, setter):
                """
                Set crystal parameters element in v3.0 based on 2D/3D crystal parameter info in v1.9

                Parameters:
                @param spec_prep_in: Object wrapping specimen preparation element in v1.9
                @param setter: Setter function for setting crystal parameters object in v_in 3.0
                """
                if spec_prep_in is not None:
                    cryst_par_in = spec_prep_in.get_twoDCrystalParameters()
                    two_dcryst = True
                    if cryst_par_in is None:
                        cryst_par_in = spec_prep_in.get_threeDCrystalParameters()
                        if cryst_par_in is None:
                            return
                        two_dcryst = False

                    cryst_par = emdb30.crystal_parameters_type()
                    if two_dcryst:
                        self.check_set(cryst_par_in.get_planeGroup, cryst_par.set_plane_group)
                    else:
                        self.check_set(cryst_par_in.get_spaceGroup, cryst_par.set_space_group)
                    # XSD: <xs:complexType name="unit_cell_type"> has 7 elements
                    unit_cell = emdb30.unit_cell_type()
                    # element 1 - <xs:complexType name="unit_cell_type">
                    # XSD: <xs:element name="a" type="cell_type"/>
                    self.set_value_and_units(cryst_par_in.get_aLength, unit_cell.set_a, emdb30.cell_type, units=const.U_ANG)
                    # element 2 - <xs:complexType name="unit_cell_type">
                    # XSD: <xs:element name="b" type="cell_type"/>
                    self.set_value_and_units(cryst_par_in.get_bLength, unit_cell.set_b, emdb30.cell_type, units=const.U_ANG)
                    # element 3 - <xs:complexType name="unit_cell_type">
                    # XSD: <xs:element name="c" type="cell_type" minOccurs="0"/>
                    self.set_value_and_units(cryst_par_in.get_cLength, unit_cell.set_c, emdb30.cell_type, units=const.U_ANG)
                    # element 4 - <xs:complexType name="unit_cell_type">
                    # XSD: <xs:element name="c_sampling_length" type="cell_type" minOccurs="0"/>
                    # element 5 - <xs:complexType name="unit_cell_type">
                    # XSD: <xs:element name="gamma" type="cell_angle_type"/>
                    self.set_value_and_units(cryst_par_in.get_gamma, unit_cell.set_gamma, emdb30.cell_angle_type, units=const.U_DEG)
                    # element 6 - <xs:complexType name="unit_cell_type">
                    # XSD: <xs:element name="alpha" type="cell_angle_type" minOccurs="0"/>
                    self.set_value_and_units(cryst_par_in.get_alpha, unit_cell.set_alpha, emdb30.cell_angle_type, units=const.U_DEG)
                    # element 7 - <xs:complexType name="unit_cell_type">
                    # XSD: <xs:element name="beta" type="cell_angle_type" minOccurs="0"/>
                    self.set_value_and_units(cryst_par_in.get_beta, unit_cell.set_beta, emdb30.cell_angle_type, units=const.U_DEG)

                    if unit_cell.hasContent_():
                        cryst_par.set_unit_cell(unit_cell)
                    if cryst_par.hasContent_():
                        setter(cryst_par)

            def set_helical_symmetry(spec_prep_in, rec, symm):
                """
                1.9 -> 3.0: Set helical symmetry parameters of reconstruction.

                Parameters:
                @param spec_prep_in: Object wrapping specimen preparation element in v1.9
                @param rec: Reconstruction object (v3.0) assumed to have [set/get]_applied_symmetry methods
                """
                if spec_prep_in is not None:
                    # XSD: <xs:complexType name="helical_parameters_type"> has 3 elements
                    hx_par = emdb30.helical_parameters_type()

                    hx_par_in = spec_prep_in.get_helicalParameters()
                    if hx_par_in is not None:
                        # element 1 - <xs:complexType name="helical_parameters_type">
                        # XSD: <xs:element name="delta_z" minOccurs="0">
                        self.set_value_and_units(hx_par_in.get_deltaZ, hx_par.set_delta_z, emdb30.delta_zType, units=const.U_ANG)
                        # element 2 - <xs:complexType name="helical_parameters_type">
                        # XSD: <xs:element name="delta_phi" minOccurs="0">
                        self.set_value_and_units(hx_par_in.get_deltaPhi, hx_par.set_delta_phi, emdb30.delta_phiType, units=const.U_DEG)
                        # element 3 - <xs:complexType name="helical_parameters_type">
                        # XSD: <xs:element name="axial_symmetry" minOccurs="0">
                        self.check_set(hx_par_in.get_axialSymmetry, hx_par.set_axial_symmetry)
                        # element 3 - <xs:complexType name="helical_parameters_type">
                        # REMOVED - XSD: <xs:element name="hand" minOccurs="0"> added as required by old v1.9s
                        #hnd = hx_par_in.get_hand()
                        #hx_par.set_hand(hnd)

                    if hx_par.hasContent_():
                        symm.set_helical_parameters(hx_par)
                    rec.set_applied_symmetry(symm)

            def set_base_image_processing(im_proc, i):
                """
                Method for setting elements common to all image processing

                Parameters:
                @param: im_proc - image processing object
                @param: i - reconstruction index
                """
                # XSD: <xs:complexType name="base_image_processing_type"> has 2 elements and 1 attribute
                # attribute 1 - <xs:complexType name="base_image_processing_type">
                im_proc.set_image_processing_id(i)
                # XSD: <xs:attribute name="id" type="xs:positiveInteger" use="required"/>
                # element 1 - <xs:complexType name="base_image_processing_type">
                # XSD: <xs:element name="image_recording_id" type="xs:positiveInteger"/
                # im_proc.set_image_recording_processing_id(i)
                # element 2 - <xs:complexType name="base_image_processing_type">
                # XSD: <xs:element name="details" type="xs:token" minOccurs="0"/>
                # details are set in separate add ons

            def set_final_reconstruction(final_rec, reconstruction, proc, spec_prep_in=None, no_apply_symm=False):
                """
                Method that sets final reconstruction elements

                Parameters:
                @params: final_rec - final reconstruction object
                @params: reconstruction - reconstruction object from v1.9
                @params: proc - processing object from v1.9
                """
                # XSD: <xs:complexType name="final_reconstruction_type"> has 8 elements
                # element 1 - <xs:complexType name="final_reconstruction_type">
                # XSD: <xs:element name="number_classes_used" type="xs:positiveInteger" minOccurs="0"/>
                # self.check_set(?, final_rec.set_number_classes_used)
                # element 2 - <xs:complexType name="final_reconstruction_type">
                # XSD: <xs:element name="applied_symmetry" type="applied_symmetry_type" minOccurs="0">
                # XSD: <xs:complexType name="applied_symmetry_type"> has 3 elements
                symm = emdb30.applied_symmetry_type()
                if self.roundtrip:
                    # applied_symmetry_type is a sequence
                    if proc is not None and not no_apply_symm:
                        symm_in = proc.get_appliedSymmetry()
                        if symm_in is not None:
                            # element 1 - <xs:complexType name="applied_symmetry_type">
                            # XSD: <xs:element name="space_group" type="xs:token"/>
                            # element 2 - <xs:complexType name="applied_symmetry_type">
                            # XSD: <xs:element name="point_group">
                            symm.set_point_group(symm_in)
                    if spec_prep_in is not None:
                        # element 3 - <xs:complexType name="applied_symmetry_type">
                        # XSD: <xs:element name="helical_parameters" type="helical_parameters_type">
                        set_helical_symmetry(spec_prep_in, final_rec, symm)
                else:
                    # applied_symmetry_type elements are choices
                    if spec_prep_in is None:
                        symm_in = None
                        if proc is not None:
                            symm_in = proc.get_appliedSymmetry()
                        if symm_in is not None:
                            # element 1 - <xs:complexType name="applied_symmetry_type">
                            # XSD: <xs:element name="space_group" type="xs:token"/>
                            # element 2 - <xs:complexType name="applied_symmetry_type">
                            # XSD: <xs:element name="point_group">
                            symm.set_point_group(symm_in)
                    else:
                        # element 3 - <xs:complexType name="applied_symmetry_type">
                        # XSD: <xs:element name="helical_parameters" type="helical_parameters_type">
                        set_helical_symmetry(spec_prep_in, final_rec, symm)
                if symm.hasContent_():
                    final_rec.set_applied_symmetry(symm)
                # element 3 - <xs:complexType name="final_reconstruction_type">
                # XSD: <xs:element name="algorithm" type="reconstruction_algorithm_type" minOccurs="0">
#                 add_to_reconstruction_details = ''
                allowed_rec_alg = ['ALGEBRAIC (ARTS)', 'BACK PROJECTION', 'EXACT BACK PROJECTION', 'FOURIER SPACE', 'SIMULTANEOUS ITERATIVE (SIRT)']
#                 add_details = None
#                 if helical:
#                     alg = reconstruction.get_algorithm()
#                     if proc is None:
#                         h_proc = process_in.get_singleParticle()
#                         self.check_set(h_proc.get_numProjections, final_rec.set_number_images_used)
#                         self.check_set(h_proc.get_numClassAverages, final_rec.set_number_classes_used)
#                         final_rec.set_algorithm(const.SP_TAG + (alg or ''))
#                         add_details = const.SP_TAG + ' ' + h_proc.get_details()
#                     else:
#                         final_rec.set_algorithm(const.HEL_TAG + (alg or ''))
#                 else:
                rec_alg = reconstruction.get_algorithm()
                if rec_alg is not None and rec_alg in allowed_rec_alg:
                    final_rec.set_algorithm(rec_alg)
                #else:
                #add_to_reconstruction_details = 'Algorithm given: %s. ' % rec_alg
                # element 4 - <xs:complexType name="final_reconstruction_type">
                # XSD: <xs:element name="resolution" minOccurs="0">
                resolution_in = reconstruction.get_resolutionByAuthor()
                if resolution_in is not None:
                    res = emdb30.resolutionType()
                    if self.roundtrip:
                        res.set_valueOf_(resolution_in)
                    else:
                        res.set_valueOf_(float(resolution_in))
                    res.set_units(const.U_ANG)
                    res.set_res_type("BY AUTHOR")
                    final_rec.set_resolution(res)
                # element 5 - <xs:complexType name="final_reconstruction_type">
                # XSD: <xs:element name="resolution_method" minOccurs="0">
                allowed_res_methods = ['DIFFRACTION PATTERN/LAYERLINES', 'FSC 0.143 CUT-OFF', 'FSC 0.33 CUT-OFF', 'FSC 0.5 CUT-OFF', 'FSC 1/2 BIT CUT-OFF', 'FSC 3 SIGMA CUT-OFF', 'OTHER']
                res_method = reconstruction.get_resolutionMethod()
                # add_to_reconstruction_details = ''
                if reconstruction is not None:
                    if res_method in allowed_res_methods:
                        self.check_set(reconstruction.get_resolutionMethod, final_rec.set_resolution_method)
                    else:
                        # add_to_reconstruction_details = 'resolutionMethod: ' + reconstruction.get_resolutionMethod()
                        final_rec.set_resolution_method('OTHER')
                # element 6 - <xs:complexType name="final_reconstruction_type">
                # XSD: <xs:element name="reconstruction_filtering" type="reconstruction_filtering_type" minOccurs="0">
                # element 7 - <xs:complexType name="final_reconstruction_type">
                # XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                soft_list = make_software_list(reconstruction.get_software())
                if soft_list is not None:
                    final_rec.set_software_list(soft_list)
                # element 8 - <xs:complexType name="final_reconstruction_type">
                # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
#                if add_to_reconstruction_details == '':
#                     self.check_set(reconstruction.get_details, final_rec.set_details)
#                else:
                rec_details = reconstruction.get_details()
                final_rec.set_details(rec_details)
#                 all_details = ''
#                 if rec_details is not None:
#                     all_details += rec_details
#                 if add_details is not None:
#                     all_details += add_details
#                 if all_details is not None:
#                    final_rec.set_details(all_details)
#                     all_details = add_to_reconstruction_details
#                     if rec_details is not None:
#                         all_details = all_details + 'Other details given: ' + rec_details
#                     final_rec.set_details(all_details)

            def set_final_angle_assignment(im_proc, reconstruction):
                """
                Method that sets angle assignment type - only details are set from v1.9 to v3.0
                """
                ang_in = reconstruction.get_eulerAnglesDetails()
                if ang_in is not None:
                    ang = emdb30.angle_assignment_type()
                    # XSD: <xs:complexType name="angle_assignment_type"> has 4 elements
                    # element 1 - <xs:complexType name="angle_assignment_type">
                    # XSD: <xs:element name="type">
                    # element 2 - <xs:complexType name="angle_assignment_type">
                    # XSD: <xs:element name="projection_matching_processing" minOccurs="0">
                    # element 3 - <xs:complexType name="angle_assignment_type">
                    # XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                    # element 4 - <xs:complexType name="angle_assignment_type">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    ang.set_details(ang_in)
                    im_proc.set_final_angle_assignment(ang)

            def adjust_for_Euler_angles_details(reconstruction, im_proc):
                """
                Helper method
                """
                ang_in = reconstruction.get_eulerAnglesDetails()
                if ang_in is not None and ang_in != '':
                    add_details = '{eulerAnglesDetails}: ' + ang_in
                    all_details = add_details
                    current_details = im_proc.get_details()
                    if current_details is not None:
                        all_details = current_details + all_details
                    im_proc.set_details(all_details)

            def set_single_particle_add_on(im_proc, process_in, reconstruction, spec_prep_in):
                """
                Method that sets all single particle image processing elements not common with other image processing

                Parameters:
                @param: im_proc - image processing object
                @param: process_in - processing from xml 1.9
                """
                sp_proc = process_in.get_singleParticle()
                # XSD: <xs:group name="single_particle_proc_add_group"> has 9 elements
                # element 1 - <xs:group name="single_particle_proc_add_group">
                # XSD: <xs:element name="particle_selection" type="particle_selection_type" maxOccurs="unbounded" minOccurs="0"/>
                # XSD: <xs:complexType name="particle_selection_type"> has 5 elements
                part_select = emdb30.particle_selection_type()
                # element 1 - <xs:complexType name="particle_selection_type">
                # XSD: <xs:element name="number_particles_selected" type="xs:positiveInteger" minOccurs="0"/>
                # element 2 - <xs:complexType name="particle_selection_type">
                # XSD: <xs:element name="reference_model" type="xs:token" minOccurs="0">
                # element 3 - <xs:complexType name="particle_selection_type">
                # XSD: <xs:element name="method" type="xs:string" minOccurs="0">
                # element 4 - <xs:complexType name="particle_selection_type">
                # XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                # element 5 - <xs:complexType name="particle_selection_type">
                # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                if part_select.hasContent_():
                    im_proc.set_particle_selection(part_select)
                # element 2 - <xs:group name="single_particle_proc_add_group">
                # XSD: <xs:element name="ctf_correction" type="ctf_correction_type" minOccurs="0"/>
                set_ctf_correction(reconstruction, im_proc)
                # element 3 - <xs:group name="single_particle_proc_add_group">
                # XSD: <xs:element name="startup_model" type="starting_map_type" minOccurs="0" maxOccurs="unbounded">
                # element 4 - <xs:group name="single_particle_proc_add_group">
                # XSD: <xs:element name="final_reconstruction" type="non_subtom_final_reconstruction_type" minOccurs="0"/>
                # XSD: <xs:complexType name="non_subtom_final_reconstruction_type"> has base and 1 element
                non_subtom_rec = emdb30.non_subtom_final_reconstruction_type()
                # base - <xs:complexType name="non_subtom_final_reconstruction_type">
                # XSD: <xs:complexType name="final_reconstruction_type"> has 8 elements
                # final_rec = emdb30.final_reconstruction_type()
                set_final_reconstruction(non_subtom_rec, reconstruction, sp_proc, spec_prep_in)
                # element 1 - <xs:complexType name="non_subtom_final_reconstruction_type">
                # XSD: <xs:element name="number_images_used" type="xs:positiveInteger" minOccurs="0">
                self.check_set(sp_proc.get_numProjections, non_subtom_rec.set_number_images_used)
                im_proc.set_final_reconstruction(non_subtom_rec)
                # element 5 - <xs:group name="single_particle_proc_add_group">
                # XSD: <xs:element name="initial_angle_assignment" type="angle_assignment_type" minOccurs="0"/>
                # element 6 - <xs:group name="single_particle_proc_add_group">
                # XSD: <xs:element name="final_angle_assignment" type="angle_assignment_type" minOccurs="0"/>
                set_final_angle_assignment(im_proc, reconstruction)
                # element 7 - <xs:group name="single_particle_proc_add_group">
                # XSD: <xs:element name="final_multi_reference_alignment" minOccurs="0">
                # element 8 - <xs:group name="single_particle_proc_add_group">
                # XSD: <xs:element name="final_two_d_classification" type="classification_type" minOccurs="0"/>
                num_cls_in = sp_proc.get_numClassAverages()
                if num_cls_in is not None:
                    final_cls = emdb30.classification_type()
                    final_cls.set_number_classes(num_cls_in)
                    if final_cls.hasContent_():
                        im_proc.set_final_two_d_classification(final_cls)
                # element 9 - <xs:group name="single_particle_proc_add_group">
                # XSD: <xs:element name="final_three_d_classification" type="classification_type" minOccurs="0"/>

                # details for base
                self.check_set(sp_proc.get_details, im_proc.set_details)

            def set_helical_add_on(im_proc, process_in, reconstruction, spec_prep_in):
                """
                1.9 -> 3.0: Method that sets all helical image processing elements not common with other image processing

                Parameters:
                @param: im_proc - image processing object
                @param: process_in - processing from xml 1.9
                @param: reconstruction - reconstruction from v1.9
                """
                h_proc = process_in.get_helical()
                # XSD: <xs:group name="helical_processing_add_group"> has 9 elements
                # element 1 - <xs:group name="helical_processing_add_group">
                # XSD: <xs:element name="final_reconstruction" type="non_subtom_final_reconstruction_type" minOccurs="0"/>
                # XSD: <xs:complexType name="non_subtom_final_reconstruction_type"> has base and 1 element
                non_subtom_rec = emdb30.non_subtom_final_reconstruction_type()
                # base - <xs:complexType name="non_subtom_final_reconstruction_type">
                # XSD: <xs:complexType name="final_reconstruction_type"> has 8 elements
                # final_rec = emdb30.final_reconstruction_type()
                set_final_reconstruction(non_subtom_rec, reconstruction, h_proc, spec_prep_in, no_apply_symm=True, helical=True)
                # element 1 - <xs:complexType name="non_subtom_final_reconstruction_type">
                # XSD: <xs:element name="number_images_used" type="xs:positiveInteger" minOccurs="0">
                im_proc.set_final_reconstruction(non_subtom_rec)
                # element 2 - <xs:group name="helical_processing_add_group">
                # XSD: <xs:element name="ctf_correction" type="ctf_correction_type" minOccurs="0"/>
                set_ctf_correction(reconstruction, im_proc)
                # element 3 - <xs:group name="helical_processing_add_group">
                # XSD: <xs:element name="segment_selection" type="segment_selection_type" maxOccurs="unbounded" minOccurs="0">
                # element 4 - <xs:group name="helical_processing_add_group">
                # XSD: <xs:element name="refinement" type="refinement_type" minOccurs="0"/>
                # element 5 - <xs:group name="helical_processing_add_group">
                # XSD: <xs:element name="startup_model" type="starting_map_type" maxOccurs="unbounded" minOccurs="0">
                # element 6 - <xs:group name="helical_processing_add_group">
                # XSD: <xs:element name="helical_layer_lines" type="layer_lines_type" minOccurs="0"/>
                # element 7 - <xs:group name="helical_processing_add_group">
                # XSD: <xs:element name="initial_angle_assignment" type="angle_assignment_type" minOccurs="0"/>
                # element 8 - <xs:group name="helical_processing_add_group">
                # XSD: <xs:element name="final_angle_assignment" type="angle_assignment_type" minOccurs="0"/>
                set_final_angle_assignment(im_proc, reconstruction)
                # element 9 - <xs:group name="helical_processing_add_group">
                # XSD: <xs:element name="crystal_parameters" type="crystal_parameters_type" minOccurs="0"/>
                set_crystal_parameters(spec_prep_in, im_proc.set_crystal_parameters)
                # details for base
                if h_proc is not None:
                    self.check_set(h_proc.get_details, im_proc.set_details)

            def set_tomography_add_on(im_proc, process_in, reconstruction, spec_prep_in):
                """
                Method that sets all tomography image processing elements not common with other image processing

                Parameters:
                @param: im_proc - image processing object
                @param: process_in - processing from xml 1.9
                @param: reconstruction - reconstruction from v1.9
                """
                tom_proc = process_in.get_tomography()
                # XSD: <xs:group name="tomography_proc_add_group"> has 4 elements
                # element 1 - <xs:group name="tomography_proc_add_group">
                # XSD: <xs:element name="final_reconstruction" type="non_subtom_final_reconstruction_type" minOccurs="0"/>
                # XSD: <xs:complexType name="non_subtom_final_reconstruction_type"> has base and 1 element
                non_subtom_rec = emdb30.non_subtom_final_reconstruction_type()
                # base - <xs:complexType name="non_subtom_final_reconstruction_type">
                # XSD: <xs:complexType name="final_reconstruction_type"> has 8 elements
                set_final_reconstruction(non_subtom_rec, reconstruction, tom_proc, spec_prep_in)
                # element 1 - <xs:complexType name="non_subtom_final_reconstruction_type">
                # XSD: <xs:element name="number_images_used" type="xs:positiveInteger" minOccurs="0">
                self.check_set(tom_proc.get_numSections, non_subtom_rec.set_number_images_used)
                im_proc.set_final_reconstruction(non_subtom_rec)
                # element 2 - <xs:group name="tomography_proc_add_group">
                # XSD: <xs:element name="series_aligment_software_list" type="software_list_type" minOccurs="0"/>
                # element 3 - <xs:group name="tomography_proc_add_group">
                # XSD: <xs:element name="ctf_correction" type="ctf_correction_type" minOccurs="0"/>
                set_ctf_correction(reconstruction, im_proc)
                # element 4 - <xs:group name="tomography_proc_add_group">
                # XSD: <xs:element name="crystal_parameters" type="crystal_parameters_type" minOccurs="0"/>
                set_crystal_parameters(spec_prep_in, im_proc.set_crystal_parameters)
                # details for base
                self.check_set(tom_proc.get_details, im_proc.set_details)
                # in case <eulerAnglesDetails> exists
                adjust_for_Euler_angles_details(reconstruction, im_proc)

            def set_subtomography_add_on(im_proc, process_in, reconstruction, spec_prep_in):
                """
                Method that sets all subtomogram image processing elements not common with other image processing

                Parameters:
                @param: im_proc - image processing object
                @param: process_in - processing from xml 1.9
                @param: reconstruction - reconstruction from v1.9
                """
                st_proc = process_in.get_subtomogramAveraging()
                # XSD: <xs:group name="subtomogram_averaging_proc_add_group"> has 7 elements
                # element 1 - <xs:group name="subtomogram_averaging_proc_add_group">
                # XSD: <xs:element name="final_reconstruction" type="subtomogram_final_reconstruction_type" minOccurs="0"/>
                # XSD: <xs:complexType name="subtomogram_final_reconstruction_type"> has base and 1 element
                subtom_rec = emdb30.subtomogram_final_reconstruction_type()
                # base - <xs:complexType name="subtomogram_final_reconstruction_type">
                # XSD: <xs:extension base="final_reconstruction_type">
                set_final_reconstruction(subtom_rec, reconstruction, st_proc, spec_prep_in)
                # element 1 - <xs:complexType name="subtomogram_final_reconstruction_type">
                # XSD: <xs:element name="number_subtomograms_used" type="xs:positiveInteger" minOccurs="0">
                self.check_set(st_proc.get_numSubtomograms, subtom_rec.set_number_subtomograms_used)
                im_proc.set_final_reconstruction(subtom_rec)
                # element 2 - <xs:group name="subtomogram_averaging_proc_add_group">
                # XSD: <xs:element name="extraction"> has 6 elements
                subtom_extr = emdb30.extractionType()
                # element 1 - <xs:element name="extraction">
                # XSD: <xs:element name="number_tomograms" type="xs:positiveInteger"/>
                # element 2 - <xs:element name="extraction">
                # XSD: <xs:element name="number_images_used" type="xs:positiveInteger"/>
                # element 3 - <xs:element name="extraction">
                # XSD: <xs:element name="reference_model" type="xs:token" minOccurs="0"/>
                # element 4 - <xs:element name="extraction">
                # XSD: <xs:element name="method" type="xs:string" minOccurs="0"/>
                # element 5 - <xs:element name="extraction">
                # XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                # element 6 - <xs:element name="extraction">
                # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                if subtom_extr.hasContent_():
                    im_proc.set_extraction(subtom_extr)
                # element 3 - <xs:group name="subtomogram_averaging_proc_add_group">
                # XSD: <xs:element name="ctf_correction" type="ctf_correction_type" minOccurs="0"/>
                set_ctf_correction(reconstruction, im_proc)
                # element 4 - <xs:group name="subtomogram_averaging_proc_add_group">
                # XSD: <xs:element name="final_multi_reference_alignment" minOccurs="0">
                # element 5 - <xs:group name="subtomogram_averaging_proc_add_group">
                # XSD: <xs:element name="final_three_d_classification" type="classification_type" minOccurs="0"/>
                num_cls_in = st_proc.get_numClassAverages()
                if num_cls_in is not None:
                    final_cls = emdb30.classification_type()
                    final_cls.set_number_classes(num_cls_in)
                    if final_cls.hasContent_():
                        im_proc.set_final_three_d_classification(final_cls)
                # element 6 - <xs:group name="subtomogram_averaging_proc_add_group">
                # XSD: <xs:element name="final_angle_assignment" type="angle_assignment_type" minOccurs="0"/>
                set_final_angle_assignment(im_proc, reconstruction)
                # element 7 - <xs:group name="subtomogram_averaging_proc_add_group">
                # XSD: <xs:element name="crystal_parameters" type="crystal_parameters_type" minOccurs="0"/>
                set_crystal_parameters(spec_prep_in, im_proc.set_crystal_parameters)

                # details for base
                self.check_set(st_proc.get_details, im_proc.set_details)

            def set_crystallography_add_on(im_proc, process_in, reconstruction, spec_prep_in):
                """
                Method that sets all crystallography image processing elements not common with other image processing

                Parameters:
                @param: im_proc - image processing object
                @param: process_in - processing from xml 1.9
                @param: reconstruction - reconstruction from v1.9
                """
                x_proc = process_in.get_twoDCrystal()
                # XSD: <xs:group name="crystallography_proc_add_group"> has 9 elements
                # element 1 - <xs:group name="crystallography_proc_add_group">
                # XSD: <xs:element name="final_reconstruction" type="non_subtom_final_reconstruction_type" minOccurs="0"/>
                # XSD: <xs:complexType name="non_subtom_final_reconstruction_type"> has base and 1 element
                non_subtom_rec = emdb30.non_subtom_final_reconstruction_type()
                # base - <xs:complexType name="non_subtom_final_reconstruction_type">
                # XSD: <xs:complexType name="final_reconstruction_type"> has 8 elements
                set_final_reconstruction(non_subtom_rec, reconstruction, x_proc, spec_prep_in, no_apply_symm=True)
                # element 1 - <xs:complexType name="non_subtom_final_reconstruction_type">
                # XSD: <xs:element name="number_images_used" type="xs:positiveInteger" minOccurs="0">
                # self.check_set(x_proc.get_numSections, non_subtom_rec.set_number_images_used)
                im_proc.set_final_reconstruction(non_subtom_rec)
                # element 2 - <xs:group name="crystallography_proc_add_group">
                # XSD: <xs:element name="crystal_parameters" type="crystal_parameters_type" minOccurs="0"/>
                set_crystal_parameters(spec_prep_in, im_proc.set_crystal_parameters)
                # element 3 - <xs:group name="crystallography_proc_add_group">
                # XSD: <xs:element name="startup_model" type="starting_map_type" maxOccurs="unbounded" minOccurs="0">
                # element 4 - <xs:group name="crystallography_proc_add_group">
                # XSD: <xs:element name="ctf_correction" type="ctf_correction_type" minOccurs="0"/>
                set_ctf_correction(reconstruction, im_proc)
                # element 5 - <xs:group name="crystallography_proc_add_group">
                # XSD: <xs:element name="molecular_replacement" type="molecular_replacement_type" minOccurs="0" >
                # element 6 - <xs:group name="crystallography_proc_add_group">
                # XSD: <xs:element name="lattice_distortion_correction_software_list" type="software_list_type" minOccurs="0"/>
                # element 7 - <xs:group name="crystallography_proc_add_group">
                # XSD: <xs:element name="symmetry_determination_software_list" type="software_list_type" minOccurs="0"/>
                # element 8 - <xs:group name="crystallography_proc_add_group">
                # XSD: <xs:element name="merging_software_list" type="software_list_type" minOccurs="0"/>
                # element 9 - <xs:group name="crystallography_proc_add_group">
                # XSD: <xs:element name="crystallography_statistics" type="crystallography_statistics_type" minOccurs="0"/>

                # details for base
                self.check_set(x_proc.get_details, im_proc.set_details)

            recon_in = process_in.get_reconstruction()
            reconstruction_index = 1
            spec_prep_in = exp_in.get_specimenPreparation()
            for reconstruction in recon_in:
                if em_method == const.EMM_SP:
                    # substitution group 1 - <xs:element ref="image_processing" maxOccurs="unbounded">
                    # XSD: <xs:element name="singleparticle_processing" substitutionGroup="image_processing" type="singleparticle_processing_type"/> has a base and add on
                    im_proc = emdb30.singleparticle_processing_type()
                    im_proc.original_tagname_ = 'singleparticle_processing'
                    # base - <xs:extension base="base_image_processing_type">
                    set_base_image_processing(im_proc, reconstruction_index)
                    # add on - <xs:group ref="single_particle_proc_add_group"/>
                    set_single_particle_add_on(im_proc, process_in, reconstruction, spec_prep_in)
                    struct_det.add_image_processing(im_proc)

                elif em_method == const.EMM_HEL:
                    # substitution group 2 - <xs:element ref="image_processing" maxOccurs="unbounded">
                    # XSD: <xs:element name="helical_processing" substitutionGroup="image_processing" type="helical_processing_type"/> has a base and add on
                    im_proc = emdb30.helical_processing_type()
                    im_proc.original_tagname_ = 'helical_processing'
                    # base - <xs:extension base="base_image_processing_type">
                    set_base_image_processing(im_proc, reconstruction_index)
                    # add on - <xs:group ref="helical_processing_add_group"/>
                    set_helical_add_on(im_proc, process_in, reconstruction, spec_prep_in)
                    struct_det.add_image_processing(im_proc)

                elif em_method == const.EMM_TOM:
                    # substitution group 3 - <xs:element ref="image_processing" maxOccurs="unbounded">
                    # XSD: <xs:element name="tomography_processing" substitutionGroup="image_processing" type="tomography_processing_type"/>
                    im_proc = emdb30.tomography_processing_type()
                    im_proc.original_tagname_ = 'tomography_processing'
                    # base - <xs:extension base="base_image_processing_type">
                    set_base_image_processing(im_proc, reconstruction_index)
                    # add on - <xs:group ref="tomography_processing_add_group"/>
                    set_tomography_add_on(im_proc, process_in, reconstruction, spec_prep_in)
                    struct_det.add_image_processing(im_proc)

                elif em_method == const.EMM_STOM:
                    # substitution group 4 - <xs:element ref="image_processing" maxOccurs="unbounded">
                    # XSD: <xs:element name="subtomogram_averaging_processing" substitutionGroup="image_processing" type="subtomogram_averaging_processing_type"/>
                    im_proc = emdb30.subtomogram_averaging_processing_type()
                    im_proc.original_tagname_ = 'subtomogram_averaging_processing'
                    # base - <xs:extension base="base_image_processing_type">
                    set_base_image_processing(im_proc, reconstruction_index)
                    # add on - <xs:group ref="subtomogram_averaging_processing_add_group"/>
                    set_subtomography_add_on(im_proc, process_in, reconstruction, spec_prep_in)
                    struct_det.add_image_processing(im_proc)

                elif em_method == 'twoDCrystal':
                    # substitution group 5 - <xs:element ref="image_processing" maxOccurs="unbounded">
                    # XSD: <xs:element name="crystallography_processing" substitutionGroup="image_processing" type="crystallography_processing_type"/>
                    im_proc = emdb30.crystallography_processing_type()
                    im_proc.original_tagname_ = 'crystallography_processing'
                    # base - <xs:extension base="base_image_processing_type">
                    set_base_image_processing(im_proc, reconstruction_index)
                    # add on - <xs:group ref="crystallography_proc_add_group"/>
                    set_crystallography_add_on(im_proc, process_in, reconstruction, spec_prep_in)
                    struct_det.add_image_processing(im_proc)

                reconstruction_index = reconstruction_index + 1

        sd_list.add_structure_determination(struct_det)

        xml_out.set_structure_determination_list(sd_list)

        # element 5 - <xs:complexType name="entry_type">
        # XSD: <xs:element name="map" type="map_type">
        map_out = emdb30.map_type()
        map_in = xml_in.get_map()
        spec_prep_in = exp_in.get_specimenPreparation()
        if map_in is not None:
            copy_map_19_to_30(map_in, map_out, is_map=True, spec_prep_in=spec_prep_in)
        xml_out.set_map(map_out)
        # element 6 - <xs:complexType name="entry_type">
        # XSD: <xs:element name="interpretation" type="interpretation_type" minOccurs="0"/>
        # XSD: <xs:complexType name="interpretation_type"> has 6 elements
        intrp = None
        # element 1 - <xs:complexType name="interpretation_type">
        # XSD: <xs:element name="modelling_list" minOccurs="0"> has 1 element
        fitting_list_in = exp_in.get_fitting()
        if fitting_list_in is not None and len(fitting_list_in) > 0:
            intrp = emdb30.interpretation_type()
            modelling_list = emdb30.modelling_listType()
            for fit in fitting_list_in:
                # element 1 - <xs:element name="modelling_list" minOccurs="0">
                # XSD: <xs:element name="modelling" type="modelling_type" maxOccurs="unbounded">
                # XSD: <xs:complexType name="modelling_type"> has 8 elements
                modelling = emdb30.modelling_type()
                PDB_entry_id_given = False
                pdb_list_in = fit.get_pdbEntryIdList()
                if pdb_list_in is not None:
                    pdb_in = pdb_list_in.get_pdbEntryId()
                    chains_in = pdb_list_in.get_pdbChainId()
                    if len(pdb_in) > 0:
                        for p_in in pdb_in:
                            # element 1 - <xs:complexType name="modelling_type">
                            # XSD: <xs:element name="initial_model" maxOccurs="unbounded" minOccurs="0"> has 3 elements
                            pdb_model = emdb30.initial_modelType()
                            # element 1 - <xs:element name="initial_model" maxOccurs="unbounded" minOccurs="0">
                            # XSD: <xs:element name="access_code">
                            pdb_model.set_access_code(p_in)
                            # element 2 - <xs:element name="initial_model" maxOccurs="unbounded" minOccurs="0">
                            # XSD: <xs:element name="chain" maxOccurs="unbounded" minOccurs="0"> extension of base="chain_type"
                            # Map all chains on a best effort basis - if it matches the pattern PDBID_CHAIN - check if PDBID matches
                            for ch_in in chains_in:
                                chain = emdb30.chainType()
                                mtch = re.match(const.PDB_CHAIN_PAT, ch_in)
                                if mtch is not None:
                                    PDB_entry_id_given = True
                                    match_groups = mtch.groups()
                                    pdb_code = match_groups[0]
                                    ch = match_groups[2]
                                    if ch != '':
                                        if pdb_code == p_in:
                                            chain.set_chain_id(ch)
                                        else:
                                            # this is an annotation problem
                                            pass
                                else:
                                    ch_minus_commas = ch_in.replace(',', '')
                                    chain.set_chain_id(ch_minus_commas)
                                if chain.hasContent_():
                                    pdb_model.add_chain(chain)
                            # element 3 - <xs:element name="initial_model" maxOccurs="unbounded" minOccurs="0">
                            # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>

                            modelling.add_initial_model(pdb_model)
                    elif len(chains_in) > 0:
                        # Pathological case when chains are specified but no PDB entry
                        # in this case use the first element as the PDB entry - try and parse the first chain element to see if the PDB ID is embedded
                        # element 1 - <xs:complexType name="modelling_type">
                        # XSD: <xs:element name="initial_model" maxOccurs="unbounded" minOccurs="0"> has 3 elements
                        pdb_model = emdb30.initial_modelType()
                        self.warn(1, "Chain IDs specified but no PDB ID! Will try and parse PDB ID from first chain ID!")
                        chain_in = chains_in[0]
                        mtch = re.match(const.PDB_CHAIN_PAT, chain_in)
                        if mtch is not None:
                            match_groups = mtch.groups()
                            pdb_code = match_groups[0]
                            chain_in = match_groups[2]
                        else:
                            pdb_code = chain_in
                        # element 1 - <xs:element name="initial_model" maxOccurs="unbounded" minOccurs="0">
                        # XSD: <xs:element name="access_code">
                        pdb_model.set_access_code(pdb_code)
                        # element 2 - <xs:element name="initial_model" maxOccurs="unbounded" minOccurs="0">
                        # XSD: <xs:element name="chain" maxOccurs="unbounded" minOccurs="0"> extension of base="chain_type"
                        chain = emdb30.chainType()
                        chain.set_id(chain_in)
                        pdb_model.add_chain(chain)
                        # Rest of chains
                        for const in chains_in[1:]:
                            chain = emdb30.chainType()
                            chain.set_id(const)
                            pdb_model.add_chain(chain)
                        # element 3 - <xs:element name="initial_model" maxOccurs="unbounded" minOccurs="0">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        modelling.add_initial_model(pdb_model)

                # element 2 - <xs:complexType name="modelling_type">
                # XSD: <xs:element name="final_model" minOccurs="0">
                # element 3 - <xs:complexType name="modelling_type">
                # XSD: <xs:element name="refinement_protocol" minOccurs="0">
                allowed_rec_protocols = ['AB INITIO MODEL', 'BACKBONE TRACE', 'FLEXIBLE FIT', 'OTHER', 'RIGID BODY FIT']
                ref_prot = fit.get_refProtocol()
                if ref_prot is not None:
                    if ref_prot in const.FITTING_19_to_30:
                        ref_prot = const.FITTING_19_to_30.get(ref_prot)
                    if self.roundtrip:
                        modelling.set_refinement_protocol(ref_prot)
                    else:
                        if ref_prot in allowed_rec_protocols:
                            self.check_set(fit.get_refProtocol, modelling.set_refinement_protocol)
                        else:
                            modelling.set_refinement_protocol('OTHER')
                # element 4 - <xs:complexType name="modelling_type">
                # XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                soft_list = make_software_list(fit.get_software())
                if soft_list is not None:
                    modelling.set_software_list(soft_list)
                # element 5 - <xs:complexType name="modelling_type">
                # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                if PDB_entry_id_given:
                    add_details = 'PDBEntryID_givenInChain. '
                    all_details = add_details
                    fit_details = fit.get_details()
                    if fit_details is not None:
                        all_details = add_details + fit_details
                    modelling.set_details(all_details)
                else:
                    self.check_set(fit.get_details, modelling.set_details)
                # element 6 - <xs:complexType name="modelling_type">
                # XSD: <xs:element name="target_criteria" type="xs:token" minOccurs="0">
                self.check_set(fit.get_targetCriteria, modelling.set_target_criteria)
                # element 7 - <xs:complexType name="modelling_type">
                # XSD: <xs:element name="refinement_space" type="xs:token" minOccurs="0">
                self.check_set(fit.get_refSpace, modelling.set_refinement_space)
                # element 8 - <xs:complexType name="modelling_type">
                # XSD: <xs:element name="overall_bvalue" type="xs:float" minOccurs="0">
                self.check_set(fit.get_overallBValue, modelling.set_overall_bvalue)

                modelling_list.add_modelling(modelling)
            intrp.set_modelling_list(modelling_list)

        # element 2 - <xs:complexType name="interpretation_type">
        # XSD: <xs:element name="figure_list" minOccurs="0"> has 1 element
        supp_in = xml_in.get_supplement()
        if supp_in is not None:
            if intrp is None:
                intrp = emdb30.interpretation_type()
            # element 1 - <xs:element name="figure_list" minOccurs="0">
            # XSD: <xs:element name="figure" type="figure_type" maxOccurs="unbounded"/>
            fig_list_in = supp_in.get_figureSet()
            if fig_list_in is not None:
                figs_in = fig_list_in.get_figure()
                fig_list = emdb30.figure_listType()
                for fig_in in figs_in:
                    # XSD: <xs:complexType name="figure_type"> has2 elements
                    fig = emdb30.figure_type()
                    # element 1 - <xs:complexType name="figure_type">
                    # XSD: <xs:element name="file">
                    fig.set_file(fig_in.get_file())
                    # element 2 - <xs:complexType name="figure_type">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    fig.set_details(fig_in.get_details())
                    fig_list.add_figure(fig)
                if fig_list.hasContent_():
                    intrp.set_figure_list(fig_list)
            # element 3 - <xs:complexType name="interpretation_type">
            # XSD: <xs:element name="segmentation_list" minOccurs="0"> has 1 element
            mask_set_in = supp_in.get_maskSet()
            if mask_set_in is not None:
                masks_in = mask_set_in.get_mask()
                seg_list = emdb30.segmentation_listType()
                for mask_in in masks_in:
                    # element 1 - <xs:element name="segmentation_list" minOccurs="0">
                    # XSD: <xs:element name="segmentation" maxOccurs="unbounded"> has 3 elements
                    seg = emdb30.segmentationType()
                    # element 1 - <xs:element name="segmentation" maxOccurs="unbounded">
                    # XSD: <xs:element name="file">
                    mask_in_filename = mask_in.get_file().valueOf_
                    seg.set_file(mask_in_filename)
                    # element 2 - <xs:element name="segmentation" maxOccurs="unbounded">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    # element 3 - <xs:element name="segmentation" maxOccurs="unbounded">
                    # XSD: <xs:element name="mask_details" type="map_type" minOccurs="0">
                    # deprecated
                    seg_map = emdb30.map_type()
                    copy_map_19_to_30(mask_in, seg_map)
                    seg.set_mask_details(seg_map)
                    seg_list.add_segmentation(seg)
                if seg_list.hasContent_():
                    intrp.set_segmentation_list(seg_list)
            # element 4 - <xs:complexType name="interpretation_type">
            # XSD: <xs:element name="slices_list" minOccurs="0"> has 1 element
            slice_set_in = supp_in.get_sliceSet()
            if slice_set_in is not None:
                slices_in = slice_set_in.get_slice()
                slc_list = emdb30.slices_listType()
                # element 1 - <xs:element name="slices_list" minOccurs="0">
                # XSD: <xs:element name="slice" type="map_type" maxOccurs="unbounded"/>
                for slice_in in slices_in:
                    slice_map = emdb30.map_type()
                    copy_map_19_to_30(slice_in, slice_map)
                    if slice_map.hasContent_():
                        slc_list.add_slice(slice_map)
                if slc_list.hasContent_():
                    intrp.set_slices_list(slc_list)
            # element 5 - <xs:complexType name="interpretation_type">
            # XSD: <xs:element name="additional_map_list" minOccurs="0">
            # element 6 - <xs:complexType name="interpretation_type">
            # XSD: <xs:element name="half_map_list" minOccurs="0">
        if intrp is not None and intrp.hasContent_():
            xml_out.set_interpretation(intrp)
        # element 7 - <xs:complexType name="entry_type">
        # XSD: <xs:element name="validation" minOccurs="0">
        if supp_in is not None:
            fsc_set_in = supp_in.get_fscSet()
            if fsc_set_in is not None:
                valid_list = emdb30.validationType()
                fsc_list_in = fsc_set_in.get_fsc()
                if fsc_list_in is not None:
                    for fsc_in in fsc_list_in:
                        # XSD: <xs:complexType name="fsc_curve_validation_type"> extends <xs:complexType name="validation_type">
                        # XSD: <xs:complexType name="validation_type"> has 2 elements
                        fsc = emdb30.fsc_curve_validation_type()
                        fsc.original_tagname_ = 'fsc_curve'
                        # element 1 - <xs:complexType name="validation_type">
                        # XSD: <xs:element name="file">
                        fsc.set_file(fsc_in.get_file())
                        # element 2 - <xs:complexType name="validation_type">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        self.check_set(fsc_in.get_details, fsc.set_details)
                        valid_list.add_validation_method(fsc)
                    if valid_list.hasContent_():
                        xml_out.set_validation(valid_list)

        # Write XML to file
        file_out = open(output_file, 'w') if output_file else sys.stdout
        file_out.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        xml_out.export(file_out, 0, name_='emd')

        if file_out is not sys.stdout:
            file_out.close()

        # validate the output v3.0 file
        if self.validate_xml_out:
            print '--validation v3.0' + '-'*100
            if self.relaxed:
                self.validate(output_file, EMDBSettings.schema30relaxed)
            else:
                self.validate(output_file, EMDBSettings.schema30)
            print '-'*117

    def translate_1_9_to_1_9(self, input_file, output_file):
        """
        Convert input file from 1.9 to 1.9 schema. This gets tags into the same ordering as other >1.9 converters.

        Parameters:
        @param input_file: Name of input file
        @param output_file: Name of output file
        """
        xml_out = emdb_19.parse(input_file, silence=True)
        # Write XML to file
        file_hdl = open(output_file, 'w') if output_file else sys.stdout
        file_hdl.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        xml_out.export(file_hdl, 0, name_='emdEntry')

        if file_hdl is not sys.stdout:
            file_hdl.close()

    def translate_3_0_to_1_9(self, input_file, output_file):
        """
        Convert input file from 3.0 to 1.9 schema

        Parameters:
        @param input_file: Name of input file
        @param output_file: Name of output file
        """
        const = self.Constants
        emdb30 = emdb_30
        if self.relaxed:
            emdb30 = emdb_30relaxed

        # def set_el_single_particle(imp_in, final_reconstruct_in, proc_spec):
        #     """
        #     Helper function - groups calls
        #     """
        #     # element 1 - <xs:complexType name="singPartType">
        #     # XSD: <xs:element name="appliedSymmetry" type="pointGroupSymmetryType" minOccurs="0"/>
        #     symm_in = final_reconstruct_in.get_applied_symmetry()
        #     if symm_in is not None:
        #         self.check_set(symm_in.get_point_group, proc_spec.set_appliedSymmetry)
        #     # element 2 - <xs:complexType name="singPartType">
        #     # XSD: <xs:element name="numProjections" type="xs:positiveInteger" minOccurs="0"/
        #     self.check_set(final_reconstruct_in.get_number_images_used, proc_spec.set_numProjections)
        #     # element 3 - <xs:complexType name="singPartType">
        #     # XSD: <xs:element name="numClassAverages" type="xs:positiveInteger" minOccurs="0"/>
        #     sp_cls_in = imp_in.get_final_two_d_classification()
        #     if sp_cls_in is not None:
        #         self.check_set(sp_cls_in.get_number_classes, proc_spec.set_numClassAverages)
        #     # element 4 - <xs:complexType name="singPartType">
        #     # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
        #     self.check_set(imp_in.get_details, proc_spec.set_details)

        def add_external_references(ref_in, ref_out):
            """
            Copy over reference list for journals or non-journals

            Parameters:
            @param ref_in: Input citation with reference list
            @param ref_out: Output citation to which reference list is added.
            """
            ext_refs_in = ref_in.get_external_references()
            for ext_ref_in in ext_refs_in:
                ext_ref_out = emdb_19.externalRefType()
                ref_type = ext_ref_in.get_type()
                if ref_type is not None:
                    ref_type_low = ref_type.lower()
                ext_ref_out.set_type(ref_type_low)
                ext_ref_out.set_valueOf_(ext_ref_in.get_valueOf_())
                ref_out.add_externalReference(ext_ref_out)

        def get_authors(auth_list_in, simple=False):
            """
            Get authors from 3.0 -> 1.9 while reformatting them and creating a string

            Parameters
            @param auth_list_in: list object of 3.0 author objects
            @param simple: boolean - True means that the authors in 3.0 are simple strings, otherwise they are journal authors
            @return:
            """
            auth_list = []
            for auth_in in auth_list_in:
                if simple:
                    ang_in_details = auth_in
                else:
                    ang_in_details = auth_in.get_valueOf_()
                auth_list.append(ang_in_details)
#                 authCompIn = ang_in_details.split(', ')
#                 lenAuthCompIn = len(authCompIn)
#                 if lenAuthCompIn < 2:
#                     self.warn(1, "Author name has more (or less) than two comma separated strings (%d) - will be ignored!" % lenAuthCompIn)
#                 else:
#                     auth_list.append('%slc_in %slc_in' % (authCompIn[0], authCompIn[1].strip('.')))
            if len(auth_list) > 0:
                auth_str = ', '.join(auth_list)
            else:
                auth_str = ''
            return auth_str

        def copy_citation(cite_in, cite_out):
            """
            Copy over citation from 3.0 to 1.9

            Parameters:
            @param cite_in: Input citation in 3.0 schema: prRefType
            @param cite_out: Output citation in 1.9 schema
            """
            # XSD: <xs:complexType name="prRefType"> is ..
            # an extension of <xs:extension base="pubType"> and has a
            # ... 1 attribute
            ref_in = cite_in.get_citation_type()
            if ref_in is not None and ref_in.hasContent_():
                if ref_in.original_tagname_ == 'journal_citation':
                    # XSD: <xs:complexType name="jrnlArtType"> has
                    # ... 8 elements
                    jrnl = emdb_19.jrnlArtType()
                    # element 1 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="authors" type="xs:string"/>
                    jrnl.set_authors(get_authors(ref_in.get_author()))
                    # element 2 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="articleTitle" type="xs:string"/>
                    jrnl.set_articleTitle(ref_in.get_title())
                    # element 3 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="journal" type="xs:string"/>
                    if self.roundtrip:
                        jrnl_name = ref_in.get_journal()
                    else:
                        jrnl_name = ref_in.get_journal() or ref_in.get_journal_abbreviation() or 'n/a'
                    if jrnl_name is not None:
                        jrnl.set_journal(jrnl_name)
                    # element 4 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="volume" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_volume, jrnl.set_volume)
                    # element 5 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="firstPage" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_first_page, jrnl.set_firstPage)
                    # element 6 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="lastPage" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_last_page, jrnl.set_lastPage)
                    # element 7 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="year" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_year, jrnl.set_year)
                    # element 8 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="externalReference" type="externalRefType" minOccurs="0" maxOccurs="unbounded"/>
                    add_external_references(ref_in, jrnl)
                    # attribute 1 - <xs:complexType name="prRefType">
                    # XSD: <xs:attribute name="published" type="xs:boolean" use="required"/>
                    cite_out.set_published(ref_in.get_published())

                    cite_out.set_journalArticle(jrnl)
                else:
                    # XSD: <xs:complexType name="jrnlArtType"> has
                    # .. 12 elements
                    non_jrnl = emdb_19.nonJrnlArtType()
                    # element 1 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="authors" type="xs:string"/>
                    non_jrnl.set_authors(get_authors(ref_in.get_author()))
                    # element 2 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="chapterTitle" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_chapter_title, non_jrnl.set_chapterTitle)
                    # element 3 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="book" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_title, non_jrnl.set_book)
                    # element 4 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="thesisTitle" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_thesis_title, non_jrnl.set_thesisTitle)
                    # element 5 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="editor" type="xs:string" minOccurs="0"/>
                    non_jrnl.set_editor(get_authors(ref_in.get_editor()))
                    # element 6 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="publisher" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_publisher, non_jrnl.set_publisher)
                    # element 7 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="publisherLocation" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_publisher_location, non_jrnl.set_publisherLocation)
                    # element 8 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="volume" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_volume, non_jrnl.set_volume)
                    # element 9 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="firstPage" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_first_page, non_jrnl.set_firstPage)
                    # element 10 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="lastPage" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_last_page, non_jrnl.set_lastPage)
                    # element 11 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="year" type="xs:string" minOccurs="0"/>
                    self.check_set(ref_in.get_year, non_jrnl.set_year)
                    # element 12 - <xs:complexType name="jrnlArtType">
                    # XSD: <xs:element name="externalReference" type="externalRefType" minOccurs="0" maxOccurs="unbounded"/>
                    add_external_references(ref_in, non_jrnl)
                    # attribute 1 - <xs:complexType name="prRefType">
                    # XSD: <xs:attribute name="published" type="xs:boolean" use="required"/>
                    cite_out.set_published(ref_in.get_published())

                    cite_out.set_nonJournalArticle(non_jrnl)

        def make_software_from_list(soft_list_in):
            """
            Take a software list (3.0 construct) and create a single string with software names (1.9 construct)

            Parameters:
            @param soft_list_in: software list as software_list_type (3.0)
            @return: Comma (', ') separated string of software
            """
            if soft_list_in is not None and len(soft_list_in) > 0:
                soft_str = ', '.join([soft.get_name() for soft in soft_list_in])
                return soft_str
            else:
                return None

        def copy_natural_source(src_in, src_out, cell=True, organelle=True, tissue=True, cellular_location=True, organ=True, only_common=False):
            """
            Copy natural source from 3.0 to 1.9

            Parameters:
            @param src_in: Instance of molecule or supramolecule
            @param src_out: Instance of molecule type (e.g. protein) in 1.9
            @param cell: Whether to generate cell field
            @param organelle: Whether to generate organelle field
            @param tissue: Whether to generate tissue field
            @param cellular_location: Whether to generate cellular_location field
            """
            ns_in = src_in.get_natural_source()
            # Natural source can have any number of elements
            ns_1_in = None
            if isinstance(ns_in, list) and len(ns_in) > 0:
                ns_1_in = ns_in[0]
            else:
                ns_1_in = ns_in
            if ns_1_in is not None and ns_1_in != []:
                org_in = ns_1_in.get_organism()
                if org_in is not None:
                    # element 1 - <xs:complexType name="proteinType/virusType/cellCompType/nuclAcidType/ligandType/riboTypeEu/riboTypePro">
                    # XSD: <xs:element name="sciSpeciesName" type="sciSpeciesType" minOccurs="0"/>
                    orn_val = org_in.get_valueOf_()
                    ncbi_val = org_in.get_ncbi()
                    src_out.set_sciSpeciesName(emdb_19.sciSpeciesType(valueOf_=orn_val, ncbiTaxId=ncbi_val))
                strain_in = ns_1_in.get_strain()
                if strain_in is not None:
                    # element 2 - <xs:complexType name="proteinType/virusType/cellCompType/nuclAcidType/ligandType/riboTypeEu/riboTypePro">
                    # XSD: <xs:element name="sciSpeciesStrain" type="xs:string" minOccurs="0"/>
                    src_out.set_sciSpeciesStrain(strain_in)

                # element 3 - <xs:complexType name="proteinType/virusType/cellCompType/nuclAcidType/ligandType/riboTypeEu/riboTypePro">
                # XSD: <xs:element name="synSpeciesName" type="xs:string" minOccurs="0"/>
                self.check_set(ns_1_in.get_synonym_organism, src_out.set_synSpeciesName)

                if not only_common:
                    nat_src = emdb_19.natSrcType()
                    if cell or organelle or tissue or organ or cellular_location:
                        # XSD: <xs:complexType name="natSrcType"> has 4 elements
                        # element 1 - <xs:complexType name="natSrcType">
                        # XSD: <xs:element name="cell" type="xs:string" minOccurs="0"/>
                        if cell:
                            self.check_set(ns_1_in.get_cell, nat_src.set_cell)
                        # element 2 - <xs:complexType name="natSrcType">
                        # XSD: <xs:element name="organelle" type="xs:string" minOccurs="0"/>
                        if organelle:
                            self.check_set(ns_1_in.get_organelle, nat_src.set_organelle)
                        # element 3 - <xs:complexType name="natSrcType">
                        # XSD: <xs:element name="organOrTissue" type="xs:string" minOccurs="0"/>
                        if tissue:
                            self.check_set(ns_1_in.get_tissue, nat_src.set_organOrTissue)
                        if organ:
                            self.check_set(ns_1_in.get_organ, nat_src.set_organOrTissue)
                        # element 4 - <xs:complexType name="natSrcType">
                        # XSD: <xs:element name="cellLocation" type="xs:string" minOccurs="0"/>
                        if cellular_location:
                            self.check_set(ns_1_in.get_cellular_location, nat_src.set_cellLocation)

                    src_out.set_natSource(nat_src)

        def copy_ctf_and_euler_angles(im_proc_in, rec_obj, im_proc_out):
            """
            Copy CTF and Euler angle info from 3.0 to 1.9 elements

            Parameters:
            @param im_proc_in: image_processing object from 3.0
            @param rec_obj: reconstruction object from 1.9
            @param im_proc_out: image_processing object from 1.9
            """
            ctf_in = im_proc_in.get_ctf_correction()
            if ctf_in is not None:
                ctf_in_details = ctf_in.get_details()
                if ctf_in_details is not None:
                    rec_obj.set_ctfCorrection(ctf_in_details)

            # Not all elements have a euler angle element
            try:
                ang_in = im_proc_in.get_final_angle_assignment()
                if ang_in is not None:
                    ang_in_details = ang_in.get_details()
                    if ang_in_details is not None:
                        rec_obj.set_eulerAnglesDetails(ang_in_details)

            except:
                # Check if info has been stored in details section
                details = im_proc_in.get_details()
                if details is not None:
                    if details.find('{eulerAnglesDetails}: ') != -1:
                        split_details = details.split('{eulerAnglesDetails}: ')
                        if len(split_details) == 2:
                            rec_obj.set_eulerAnglesDetails(split_details[1])
                            im_proc_out.set_details(split_details[0])

        def set_mol_weight(comp, wt_in, meth=False):
            """
            Set molecular weight if provided

            Parameters:
            @param comp: v1.9 component with setter functions for setting mW
            @param wt_in: Molecular weight object from 3.0
            @param meth: Whether set_molWtMethod should be called (true only for sample)
            """
            # XSD: <xs:element name="molWtTheo" type="mwType" minOccurs="0"/>
            if wt_in is not None:
                exp_weight = wt_in.get_experimental()
                if exp_weight is not None:
                    comp.set_molWtExp(emdb_19.mwType(valueOf_=exp_weight.get_valueOf_(), units=const.U_MDA))
                theo_weight = wt_in.get_theoretical()
                if theo_weight is not None:
                    comp.set_molWtTheo(emdb_19.mwType(valueOf_=theo_weight.get_valueOf_(), units=const.U_MDA))
                if meth:
                    self.check_set(wt_in.get_method, comp.set_molWtMethod)

        def copy_external_references(getter, setter):
            """
            Copy external references from 3.0 to 1.9

            Parameters:
            @param getter: function to get external references from 3.0
            @param setter: function to set external references from 1.9
            """
            x_refs_in = getter()
            if x_refs_in is not None and len(x_refs_in) > 0:
                x_refs = emdb_19.externalReferencesType()
                for x_ref_in in x_refs_in:
                    ref_type = x_ref_in.get_type()
                    if ref_type == 'UNIPROTKB':
                        x_refs.add_refUniProt(x_ref_in.get_valueOf_().upper())
                    elif ref_type == 'GO':
                        x_refs.add_refGo(x_ref_in.get_valueOf_())
                    elif ref_type == 'INTERPRO':
                        x_refs.add_refInterpro(x_ref_in.get_valueOf_())
                setter(x_refs)

        def set_crystal_parameters(imp_in, spec_prep):
            """
            Get v3.0 crystal parameters element and set 2D/3D crystal parameters element in v1.9

            Parameters
            @param imp_in: v3.0 object of image processing extension class which is assumed to have crystal parameters
            @param spec_prep: v1.9 specimen preparation object which will have the appropriate 2D/3D crystal parameters set
            """
            cryst_par_in = imp_in.get_crystal_parameters()
            if cryst_par_in is not None:
                two_dcryst = True
                cryst_par = None
                plane_group_in = cryst_par_in.get_plane_group()
                if plane_group_in is None:
                    space_group_in = cryst_par_in.get_space_group()
                    if space_group_in is None:
                        return
                    two_dcryst = False
                if two_dcryst:
                    # XSD: <xs:complexType name="twoDxtalParamType"> has 7 elements
                    cryst_par = emdb_19.twoDxtalParamType()
                else:
                    # XSD: <xs:complexType name="threeDxtalParamType"> has 7 elements
                    cryst_par = emdb_19.threeDxtalParamType()

                unit_cell_in = cryst_par_in.get_unit_cell()
                if unit_cell_in is not None and cryst_par is not None:
                    # element 1 - <xs:complexType name="twoDxtalParamType"> and <xs:complexType name="threeDxtalParamType">
                    # XSD: <xs:element name="aLength" type="lengthType" minOccurs="0"/>
                    self.set_value_and_units(unit_cell_in.get_a, cryst_par.set_aLength, emdb_19.lengthType, units='A')
                    # element 2 - <xs:complexType name="twoDxtalParamType"> and <xs:complexType name="threeDxtalParamType">
                    # XSD: <xs:element name="bLength" type="lengthType" minOccurs="0"/>
                    self.set_value_and_units(unit_cell_in.get_b, cryst_par.set_bLength, emdb_19.lengthType, units='A')
                    # element 3 - <xs:complexType name="twoDxtalParamType"> and <xs:complexType name="threeDxtalParamType">
                    # XSD: <xs:element name="cLength" type="lengthType" minOccurs="0"/>
                    self.set_value_and_units(unit_cell_in.get_c, cryst_par.set_cLength, emdb_19.lengthType, units="A")
                    # element 4 - <xs:complexType name="twoDxtalParamType"> and <xs:complexType name="threeDxtalParamType">
                    # XSD: <xs:element name="alpha" type="lengthType" minOccurs="0"/>
                    self.set_value_and_units(unit_cell_in.get_alpha, cryst_par.set_alpha, emdb_19.anglType, units=const.U_DEGF)
                    # element 4 - <xs:complexType name="twoDxtalParamType"> and <xs:complexType name="threeDxtalParamType">
                    # XSD: <xs:element name="beta" type="lengthType" minOccurs="0"/>
                    self.set_value_and_units(unit_cell_in.get_beta, cryst_par.set_beta, emdb_19.anglType, units=const.U_DEGF)
                    # element 4 - <xs:complexType name="twoDxtalParamType"> and <xs:complexType name="threeDxtalParamType">
                    # XSD: <xs:element name="gamma" type="lengthType" minOccurs="0"/>
                    self.set_value_and_units(unit_cell_in.get_gamma, cryst_par.set_gamma, emdb_19.anglType, units=const.U_DEGF)
                if two_dcryst:
                    # element 7 - <xs:complexType name="twoDxtalParamType">
                    # XSD: <xs:element name="planeGroup" type="plGrpType"/>
                    cryst_par.set_planeGroup(plane_group_in)

                    spec_prep.set_twoDCrystalParameters(cryst_par)
                else:
                    cryst_par.set_spaceGroup(space_group_in)

                    spec_prep.set_threeDCrystalParameters(cryst_par)

        def set_helical_parameters(final_reconstruct_in, spec_prep):
            """
            3.0 -> 1.9: Set v1.9 helical symmetry parameters of a specimen preparation object

            Parameters:
            @param final_reconstruct_in: v3.0 reconstruction object
            @param spec_prep: v1.9 specimen preparation object
            """
            if final_reconstruct_in is not None:
                symm_in = final_reconstruct_in.get_applied_symmetry()
                if symm_in is not None:
                    # XSD: <xs:complexType name="helixParamType"> has 4 elements
                    hx_par_in = symm_in.get_helical_parameters()
                    if hx_par_in is not None:
                        hx_par = emdb_19.helixParamType()
                        # element 1 - <xs:complexType name="helixParamType">
                        # XSD: <xs:element name="deltaPhi" type="anglType" minOccurs="0"/>
                        self.set_value_and_units(hx_par_in.get_delta_phi, hx_par.set_deltaPhi, emdb_19.anglType, units=const.U_DEGF)
                        # element 2 - <xs:complexType name="helixParamType">
                        # XSD: <xs:element name="deltaZ" type="lengthType" minOccurs="0"/>
                        self.set_value_and_units(hx_par_in.get_delta_z, hx_par.set_deltaZ, emdb_19.lengthType, units='A')
                        # element 3 - <xs:complexType name="helixParamType">
                        # XSD: <xs:element name="hand" type="handType" minOccurs="0"/>
                        d_phi = hx_par_in.get_delta_phi()
                        if d_phi is not None:
                            hand = None
                            if d_phi < 0:
                                hand = 'LEFT HANDED'
                            else:
                                hand = 'RIGHT HANDED'
                            if hand is not None:
                                hx_par.set_hand(hand)

                        # element 4 - <xs:complexType name="helixParamType">
                        # XSD: <xs:element name="axialSymmetry" type="xs:string" minOccurs="0"/>
                        self.check_set(hx_par_in.get_axial_symmetry, hx_par.set_axialSymmetry)
                        if self.roundtrip:
                            spec_prep.set_helicalParameters(hx_par)
                        else:
                            if hx_par.hasContent_():
                                spec_prep.set_helicalParameters(hx_par)

        def copy_map_30_to_19(map_in, map_out):
            """
            Copy map from 3.0 to 1.9

            Parameters:
            @param map_in: input 3.0 map
            @param map_out: output 1.9 map
            """
            # XSD: <xs:complexType name="mapType"> has 14 elements
            # element 1 - <xs:complexType name="mapType">
            # XSD: <xs:element name="file" type="mapFileType"/>
            map_file = emdb_19.mapFileType()
            # XSD: <xs:complexType name="mapFileType"> is an extension of base="mapNamePattern" and has 3 attributes
            # extension
            # XSD: <xs:simpleType name="mapNamePattern">, <xs:restriction base="xs:token">, <xs:pattern value="emd_\d+\.map\.gz"/>
            map_in_file = map_in.get_file()
            if map_in_file is not None:
                map_in_file_lower = map_in_file.lower()
                if map_in_file_lower != '':
                    map_file.set_valueOf_(map_in_file_lower)
                else:
                    self.warn(1, 'No file name given for a map')
                    if not self.roundtrip:
                        map_file.set_valueOf_('emd_0000.map.gz')
                    else:
                        map_file.set_valueOf_('')
            # attribute 1 - <xs:complexType name="mapFileType">
            # XSD: <xs:attribute name="type" type="xs:string" use="required" fixed="map"/>
            map_file.set_type("map")
            # attribute 2 - <xs:complexType name="mapFileType">
            # XSD: <xs:attribute name="format" type="xs:string" use="required" fixed="CCP4"/>
            map_file.set_format(map_in.get_format())
            # attribute 3 - <xs:complexType name="mapFileType">
            # XSD: <xs:attribute name="sizeKb" type="xs:positiveInteger" use="required"/>
            map_file.set_sizeKb(map_in.get_size_kbytes())
            map_out.set_file(map_file)

            # element 2 - <xs:complexType name="mapType">
            # XSD: <xs:element name="dataType" type="mapDataType"/>
            map_in_data_type = map_in.get_data_type()
            #data_type_dict_inv = {v: k for k, v in const.DATA_TYPE_DICT_19_TO_30.iteritems()}
            data_type_dict_inv = {}
            for k, v in const.DATA_TYPE_DICT_19_TO_30.iteritems():
                data_type_dict_inv[v] = k
            map_out_data_type = data_type_dict_inv.get(map_in_data_type)
            map_out.set_dataType(map_out_data_type)
            # element 3 - <xs:complexType name="mapType">
            # XSD: <xs:element name="dimensions" type="dimensionType"/>
            dim_in = map_in.get_dimensions()
            if dim_in is not None:
                num_rows = dim_in.get_row()
                num_columns = dim_in.get_col()
                num_sections = dim_in.get_sec()
                dim = emdb_19.dimensionType(numRows=num_rows, numColumns=num_columns, numSections=num_sections)
                map_out.set_dimensions(dim)
            # element 4 - <xs:complexType name="mapType">
            # XSD: <xs:element name="origin" type="originType"/>
            orig_in = map_in.get_origin()
            if orig_in is not None:
                origin_row = -1
                orig_row = orig_in.get_row()
                if orig_row is not None:
                    origin_row = int(orig_row)
                origin_col = -1
                orig_col = orig_in.get_col()
                if orig_col is not None:
                    origin_col = int(orig_col)
                orig_sec = orig_in.get_sec()
                origin_sec = -1
                if orig_sec is not None:
                    origin_sec = int(orig_sec)
                orig = emdb_19.originType(originRow=origin_row, originCol=origin_col, originSec=origin_sec)
                map_out.set_origin(orig)
            # element 5 - <xs:complexType name="mapType">
            # XSD: <xs:element name="limit" type="limitType"/>
            if dim_in is not None and orig_in is not None:
                limit_row = origin_row + num_rows - 1
                limit_col = origin_col + num_columns - 1
                limit_sec = origin_sec + num_sections - 1
                lim = emdb_19.limitType(limitRow=int(limit_row), limitCol=int(limit_col), limitSec=int(limit_sec))
                map_out.set_limit(lim)
            # element 6 - <xs:complexType name="mapType">
            # XSD: <xs:element name="spacing" type="spacingType"/>
            spc_in = map_in.get_spacing()
            if spc_in is not None:
                spc = emdb_19.spacingType(spacingRow=spc_in.get_x(), spacingCol=spc_in.get_y(), spacingSec=spc_in.get_z())
                map_out.set_spacing(spc)
            # element 7 - <xs:complexType name="mapType">
            # XSD: <xs:element name="cell" type="cellType"/>
            cell_in = map_in.get_cell()
            if cell_in is not None and cell_in.hasContent_():
                cell = emdb_19.cellType(cellA=emdb_19.cType(valueOf_=cell_in.get_a().get_valueOf_(), units='A'),
                                        cellB=emdb_19.cType(valueOf_=cell_in.get_b().get_valueOf_(), units='A'),
                                        cellC=emdb_19.cType(valueOf_=cell_in.get_c().get_valueOf_(), units='A'),
                                        cellAlpha=emdb_19.cAngleType(valueOf_=cell_in.get_alpha().get_valueOf_(), units=const.U_DEGF),
                                        cellBeta=emdb_19.cAngleType(valueOf_=cell_in.get_beta().get_valueOf_(), units=const.U_DEGF),
                                        cellGamma=emdb_19.cAngleType(valueOf_=cell_in.get_gamma().get_valueOf_(), units=const.U_DEGF))
                map_out.set_cell(cell)
            # element 8 - <xs:complexType name="mapType">
            # XSD: <xs:element name="axisOrder" type="axisOrderType"/>
            ax_in = map_in.get_axis_order()
            if ax_in is not None and ax_in.hasContent_():
                axis_order = emdb_19.axisOrderType(axisOrderFast=ax_in.get_fast().upper(), axisOrderMedium=ax_in.get_medium().upper(), axisOrderSlow=ax_in.get_slow().upper())
                map_out.set_axisOrder(axis_order)
            # element 9 - <xs:complexType name="mapType">
            # XSD: <xs:element name="statistics" type="statisticsType"/>
            map_out.set_statistics(map_in.get_statistics())
            # element 10 - <xs:complexType name="mapType">
            # XSD: <xs:element name="spaceGroupNumber" type="xs:string"/>
            symm_in = map_in.get_symmetry()
            if symm_in is not None:
                self.check_set(symm_in.get_space_group, map_out.set_spaceGroupNumber)
            # element 11 - <xs:complexType name="mapType">
            # XSD: <xs:element name="details" type="xs:string"/>
            # self.check_set(map_in.get_details, map_out.set_details)
            # element 12 - <xs:complexType name="mapType">
            # XSD: <xs:element name="pixelSpacing" type="pixelSpacingType"/>
            pix_in = map_in.get_pixel_spacing()
            if pix_in is not None and pix_in.hasContent_():
                pix = emdb_19.pixelSpacingType(emdb_19.pixType(valueOf_=pix_in.get_x().get_valueOf_(), units='A'),
                                               emdb_19.pixType(valueOf_=pix_in.get_y().get_valueOf_(), units='A'),
                                               emdb_19.pixType(valueOf_=pix_in.get_z().get_valueOf_(), units='A'))
                map_out.set_pixelSpacing(pix)
            # element 13 - <xs:complexType name="mapType">
            # XSD: <xs:element name="contourLevel" minOccurs="0">
            # In 1.9 contour level is only defined for primary map, not for masks etc
            map_details = map_in.get_details()
            if hasattr(map_out, 'set_contourLevel'):
                cntr_list_in = map_in.get_contour_list()
                if cntr_list_in is not None:
                    for cntr_in in cntr_list_in.get_contour():
                        if cntr_in.get_primary():
                            cntr = emdb_19.contourLevelType()
                            cnt_level = cntr_in.get_level()
                            if cnt_level is not None:
                                if map_details is not None:
                                    if map_details.find('{level is a whole number}') != -1:
                                        cntr.set_valueOf_(int(cnt_level))
                                    else:
                                        cntr.set_valueOf_(float(cnt_level))
                                else:
                                    cntr.set_valueOf_(float(cnt_level))
                            self.check_set(cntr_in.get_source, cntr.set_source, string.lower)
                            if cntr.hasContent_():
                                map_out.set_contourLevel(cntr)
            # element 11 - <xs:complexType name="mapType">
            # XSD: <xs:element name="details" type="xs:string"/>
            if map_details is not None:
                if map_details.find('{level is a whole number}') != -1:
                    map_details = map_details.replace('{level is a whole number}', '')
                if map_details != '':
                    map_out.set_details(map_details)
            elif not self.roundtrip:
                map_out.set_details('')
            # element 14 - <xs:complexType name="mapType">
            # XSD: <xs:element name="annotationDetails" type="xs:string" minOccurs="0"/>
            annot_details = map_in.get_annotation_details()
            if annot_details is not None:
                map_out.set_annotationDetails(annot_details)

        def create_eng_source(mol_or_smol_in):
            """
            Method that converts v3.0 recombinant expression into v1.9 engineered source endSource
            """
            eng_src = None
            rec_exp = None
            rec_exps = mol_or_smol_in.get_recombinant_expression()
            if rec_exps is not None:
                eng_src = emdb_19.engSrcType()
                if isinstance(rec_exps, list) and len(rec_exps) > 0:
                    rec_exp = rec_exps[0]
                elif rec_exps != []:
                    rec_exp = rec_exps
                # XSD: <xs:complexType name="engSrcType"> has 4 elements
                # element 1 - <xs:complexType name="engSrcType">
                # XSD: <xs:element name="expSystem" type="sciSpeciesType" minOccurs="0"/>
                if rec_exp is not None:
                    rec_org = rec_exp.get_recombinant_organism()
                    if rec_org is not None:
                        eng_src.set_expSystem(emdb_19.sciSpeciesType(valueOf_=rec_org.valueOf_, ncbiTaxId=rec_org.get_ncbi()))
                        # element 2 - <xs:complexType name="engSrcType">
                        # XSD: <xs:element name="expSystemStrain" type="xs:string" minOccurs="0"/>
                        eng_src.set_expSystemStrain(rec_exp.get_recombinant_strain())
                        # element 3 - <xs:complexType name="engSrcType">
                        # XSD: <xs:element name="expSystemCell" type="xs:string" minOccurs="0"/>
                        eng_src.set_expSystemCell(rec_exp.get_recombinant_cell())
                        # element 4 - <xs:complexType name="engSrcType">
                        # XSD: <xs:element name="vector" type="xs:string" minOccurs="0"/>
                        eng_src.set_vector(rec_exp.get_recombinant_plasmid())

                return eng_src

        def unpack_odd_details(mol_or_smol_in, comp, spec_comp):
            """
            """
            details_in = mol_or_smol_in.get_details()
            if details_in is not None:
                if details_in.find(' {comp details}: ') != -1:
                    # split on flag
                    all_details = details_in.split(' {comp details}: ')
                    if all_details[0] != '':
                        # this goes to component detials
                        comp.set_details(all_details[0])
                    if all_details[1] != '':
                        spec_comp.set_details(all_details[1])
                else:
                    # this comes from base/component only
                    comp.set_details(details_in)

        xml_in = None
        try:
            xml_in = emdb30.parse(input_file, silence=True)
        except:
            self.warn(1, 'Cannot open file: %s' % input_file)
        # XSD: <xs:complexType name="entryType"> has
        # .. 7 elements and 2 attributes
        xml_out = emdb_19.entryType()
        # attribute 1 - <xs:complexType name="entryType">
        # XSD: <xs:attribute name="accessCode" type="xs:string" use="required"/>
        if xml_in is not None:
            xml_out.set_accessCode(self.format_emdb_code(xml_in.get_emdb_id(), True))
        # attribute 2 - <xs:complexType name="entryType">
        # XSD: <xs:attribute name="version" type="xs:string" fixed="1.9.6"/>
        xml_out.set_version('1.9.6')
        adm_in = None
        if xml_in is not None:
            adm_in = xml_in.get_admin()
        # element 1 - <xs:complexType name="entryType">
        # XSD: <xs:element name="admin" type="adminType" minOccurs="1" maxOccurs="1"/>
        # XSD: <xs:complexType name="adminType"> has 1 element
        adm_out = emdb_19.adminType()
        if adm_in is not None:
            # element 1 - <xs:complexType name="adminType">
            # XSD: <xs:element name="lastUpdate" type="xs:date"/>
            dates_in = adm_in.get_key_dates()
            adm_out.set_lastUpdate(dates_in.get_update())
            xml_out.set_admin(adm_out)

            # element 2 - <xs:complexType name="entryType">
            # XSD: <xs:element name="deposition" type="depType" minOccurs="1" maxOccurs="1"/>
            # XSD: <xs:complexType name="depType"> has 18 elements
            dep = emdb_19.depType()
            # element 1 - <xs:complexType name="depType">
            # XSD: <xs:element name="status" minOccurs="1" maxOccurs="1
            status_in = adm_in.get_current_status()
            status_code_in = status_in.get_code().get_valueOf_()
            # Remap HOLD to HOLD1
            if status_code_in == const.STS_HOLD:
                status_code_in = const.STS_HOLD1
            prior_status_in = None
            if status_code_in == const.STS_OBS:
                # check status_history_list for prior status
                status_history_list_in = adm_in.get_status_history_list()
                if status_history_list_in is not None:
                    status_history_in = status_history_list_in.get_status()
                    if len(status_history_in) > 0:
                        prior_in = status_history_in[0]
                        prior_code_in = prior_in.get_code()
                        if prior_code_in is not None:
                            prior_status_in = prior_code_in.get_valueOf_()
            if prior_status_in is None:
                dep.set_status(emdb_19.statusType(valueOf_=status_code_in))
            else:
                dep.set_status(emdb_19.statusType(valueOf_=status_code_in, prior=prior_status_in))

            xml_out.set_deposition(dep)
            # element 2 - <xs:complexType name="depType">
            # XSD: <xs:element name="depositionDate" type="xs:date" minOccurs="1" maxOccurs="1"/>
            dep.set_depositionDate(dates_in.get_deposition())
            # element 3 - <xs:complexType name="depType">
            # XSD: <xs:element name="depositionSite" minOccurs="1" maxOccurs="1">
            sites_in = adm_in.get_sites()
            dep.set_depositionSite(const.PROC_SITE_30_TO_19[sites_in.get_deposition().lower()])
            # element 4 - <xs:complexType name="depType">
            # XSD: <xs:element name="processingSite" minOccurs="1" maxOccurs="1">
            proc_site_in = sites_in.get_last_processing()
            if proc_site_in is not None:
                dep.set_processingSite(const.PROC_SITE_30_TO_19[proc_site_in.lower()])
            # element 5 - <xs:complexType name="depType">
            # XSD: <xs:element name="headerReleaseDate" type="xs:date" minOccurs="1" maxOccurs="1"/>
            hdr_rel_date = dates_in.get_header_release()
            if hdr_rel_date is not None:
                dep.set_headerReleaseDate(hdr_rel_date)
            # element 6 - <xs:complexType name="depType">
            # XSD: <xs:element name="mapReleaseDate" type="xs:date" minOccurs="0" maxOccurs="1"/>
            self.check_set(dates_in.get_map_release, dep.set_mapReleaseDate)
            # element 7 - <xs:complexType name="depType">
            # XSD: <xs:element name="obsoletedDate" type="xs:date" minOccurs="0" maxOccurs="1"/>
            self.check_set(dates_in.get_obsolete, dep.set_obsoletedDate)
            # element 8 - <xs:complexType name="depType">
            # XSD: <xs:element name="supersededByList" type="emdbListType" minOccurs="0" maxOccurs="1"/>
            supersede_list_in = adm_in.get_superseded_by_list()
            if supersede_list_in is not None:
                supersede_list = emdb_19.emdbListType()
                for supersede_in in supersede_list_in.get_entry():
                    supersede_list.add_entry(supersede_in.get_entry())
                if supersede_list.hasContent_():
                    dep.set_supersededByList(supersede_list)
            # element 9 - <xs:complexType name="depType">
            # XSD: <xs:element name="replaceExistingEntry" type="xs:boolean" minOccurs="0" maxOccurs="1"/>
            self.check_set(adm_in.get_replace_existing_entry, dep.set_replaceExistingEntry)
            # element 10 - <xs:complexType name="depType">
            # XSD: <xs:element name="obsoleteList" type="emdbListType" minOccurs="0" maxOccurs="1"/>
            obs_list_in = adm_in.get_obsolete_list()
            if obs_list_in is not None:
                obs_list = emdb_19.emdbListType()
                for obs_in in obs_list_in.get_entry():
                    obs_list.add_entry(obs_in.get_entry())
                if obs_list.hasContent_():
                    # element 9 - <xs:complexType name="depType">
                    # XSD: <xs:element name="replaceExistingEntry" type="xs:boolean" minOccurs="0" maxOccurs="1"/>
                    dep.set_replaceExistingEntry(True)
                    dep.set_obsoleteList(obs_list)
            # element 11 - <xs:complexType name="depType">
            # XSD: <xs:element name="details" type="xs:string" minOccurs="0" maxOccurs="1"/>
            # <xs:element name="key_dates"> in emdb30.xsd doesn't have details
            # element 12 - <xs:complexType name="depType">
            # XSD: <xs:element name="inFrameEMDBId" type="emdbEntryIdType" minOccurs="0"/>
            x_ref_in = xml_in.get_crossreferences()
            emdb_list_in = x_ref_in.get_emdb_list()
            if emdb_list_in is not None:
                refs_in = emdb_list_in.get_emdb_reference()
                infr_list = []
                for ref_in in refs_in:
                    rel_in = ref_in.get_relationship()
                    emdb_id_in = ref_in.get_emdb_id()
                    if emdb_id_in is not None:
                        if rel_in is None:
                            # Assume full overlap
                            infr_list.append(emdb_id_in)
                        elif rel_in.get_in_frame() == 'FULLOVERLAP':
                            infr_list.append(emdb_id_in)
                if infr_list != []:
                    if len(infr_list) > 0:
                        infr_text = ', '.join(infr_list)
                        if infr_text != '' or infr_text.isspace():
                            dep.set_inFrameEMDBId(infr_text)
            # element 13 - <xs:complexType name="depType">
            # XSD: <xs:element name="title" type="xs:string" minOccurs="1" maxOccurs="1"/>
            title = adm_in.get_title()
            if title is not None:
                dep.set_title(title)
            else:
                dep.set_title('')
            # element 14 - <xs:complexType name="depType">
            # XSD: <xs:element name="authors" type="xs:string" minOccurs="1" maxOccurs="1"/>
            auth_list_in = adm_in.get_authors_list()
            dep.set_authors(get_authors(auth_list_in.get_author(), simple=True))
            # element 15 - <xs:complexType name="depType">
            # XSD: <xs:element name="keywords" type="xs:string" minOccurs="0" maxOccurs="1"/>
            self.check_set(adm_in.get_keywords, dep.set_keywords)
            # element 16 - <xs:complexType name="depType">
            # XSD: <xs:element name="fittedPDBEntryIdList" type="pdbidListType" minOccurs="0" maxOccurs="1"/>
            pdb_list_in = x_ref_in.get_pdb_list()
            if pdb_list_in is not None:
                fit_list = emdb_19.pdbidListType()
                fit_in = pdb_list_in.get_pdb_reference()
                for map_file in fit_in:
                    fit_list.add_fittedPDBEntryId(map_file.get_pdb_id())
                dep.set_fittedPDBEntryIdList(fit_list)
            # element 17 - <xs:complexType name="depType">
            # XSD: <xs:element name="primaryReference" type="prRefType" minOccurs="1" maxOccurs="1"/>
            cite_list_in = x_ref_in.get_citation_list()
            cite_in = cite_list_in.get_primary_citation()
            cite_pr = emdb_19.prRefType()
            dep.set_primaryReference(cite_pr)
            copy_citation(cite_in, cite_pr)
            # element 18 - <xs:complexType name="depType">
            # XSD: <xs:element name="secondaryReference" type="prRefType" minOccurs="0" maxOccurs="unbounded"/>
            sec_cites_in = cite_list_in.get_secondary_citation()
            if sec_cites_in is not None:
                for sec_cite_in in sec_cites_in:
                    sec_cite = emdb_19.prRefType()
                    copy_citation(sec_cite_in, sec_cite)
                    if sec_cite.hasContent_():
                        dep.add_secondaryReference(sec_cite)

        # element 3 - <xs:complexType name="entryType">
        # XSD: <xs:element name="map" type="mapType" maxOccurs="1"/>
        if xml_in is not None:
            map_in = xml_in.get_map()
            if map_in is not None:
                map_out = emdb_19.mapType()
                copy_map_30_to_19(map_in, map_out)
                xml_out.set_map(map_out)

        # element 4 - <xs:complexType name="entryType">
        # XSD: <xs:element name="supplement" type="supplType" minOccurs="0" maxOccurs="1"/>
        supp = None
        # XSD: <xs:complexType name="supplType"> has 4 elements
        intrp_in = None
        if xml_in is not None:
            intrp_in = xml_in.get_interpretation()
        if intrp_in is not None:
            supp = emdb_19.supplType()
            # element 1 - <xs:complexType name="supplType">
            # XSD: <xs:element name="maskSet" type="mskSetType" minOccurs="0"/>
            seg_list_in = intrp_in.get_segmentation_list()
            if seg_list_in is not None:
                segs_in = seg_list_in.get_segmentation()
                mask_set = emdb_19.mskSetType()
                for slc_in in segs_in:
                    m_seg_in = slc_in.get_mask_details()
                    if m_seg_in is not None:
                        mask = emdb_19.mskType()
                        copy_map_30_to_19(m_seg_in, mask)
                        mask_set.add_mask(mask)
                if mask_set.hasContent_():
                    supp.set_maskSet(mask_set)
            # element 2 - <xs:complexType name="supplType">
            # XSD: <xs:element name="sliceSet" type="slcSetType" minOccurs="0"/>
            slc_list_in = intrp_in.get_slices_list()
            if slc_list_in is not None:
                slcs_in = slc_list_in.get_slice()
                slc_set = emdb_19.slcSetType()
                for slc_in in slcs_in:
                    slc = emdb_19.slcType()
                    copy_map_30_to_19(slc_in, slc)
                    if slc.hasContent_():
                        slc_set.add_slice(slc)
                if slc_set.hasContent_():
                    supp.set_sliceSet(slc_set)
            # element 3 - <xs:complexType name="supplType">
            # XSD: <xs:element name="figureSet" type="figSetType" minOccurs="0"/>
            fig_list_in = intrp_in.get_figure_list()
            # write this out even if it is empty - this is to minimise unnecessary elements showing in the diff during round-trip conversion
            fig_set = emdb_19.figSetType()
            if fig_list_in is not None:
                figs_in = fig_list_in.get_figure()
                for map_file in figs_in:
                    fig = emdb_19.figType(map_file.get_file(), map_file.get_details())
                    if fig.hasContent_():
                        fig_set.add_figure(fig)
                if fig_set.hasContent_():
                    supp.set_figureSet(fig_set)
        # element 4 - <xs:complexType name="supplType">
        # XSD: <xs:element name="fscSet" type="fscSetType" minOccurs="0"/>
        fsc_set = emdb_19.fscSetType()
        valid_list_in = None
        if xml_in is not None:
            valid_list_in = xml_in.get_validation()
        if valid_list_in is not None:
            if supp is None:
                supp = emdb_19.supplType()
            vals_in = valid_list_in.get_validation_method()
            for val_in in vals_in:
                if val_in.original_tagname_ == 'fsc_curve':
                    fsc = emdb_19.fscType()
                    fsc_filename = val_in.get_file()
                    fsc.set_file(fsc_filename)
                    self.check_set(val_in.get_details, fsc.set_details)
                    fsc_set.add_fsc(fsc)
            if fsc_set.hasContent_():
                supp.set_fscSet(fsc_set)
        if supp is not None and supp.hasContent_():
            xml_out.set_supplement(supp)

        # element 5 - <xs:complexType name="entryType">
        # XSD: <xs:element name="sample" type="samplType" maxOccurs="1"/>
        sample_in = None
        if xml_in is not None:
            sample_in = xml_in.get_sample()
        if sample_in is not None:
            # XSD: <xs:complexType name="samplType"> has 8 elements
            sample = emdb_19.samplType()
            xml_out.set_sample(sample)
            # element 1 - <xs:complexType name="samplType">
            # XSD: <xs:element name="numComponents" type="xs:positiveInteger"/>
            sup_mols_in = []
            sup_mol_list_in = sample_in.get_supramolecule_list()
            if sup_mol_list_in is not None:
                sup_mols_in = sup_mol_list_in.get_supramolecule()
                n_sup_mols_in = len(sup_mols_in)
            else:
                n_sup_mols_in = 0
            mols_in = []
            mol_list_in = sample_in.get_macromolecule_list()
            if mol_list_in is not None:
                mols_in = mol_list_in.get_macromolecule()
                n_mols_in = len(mols_in)
            else:
                n_mols_in = 0
            num_comp_in = n_sup_mols_in + n_mols_in
            num_comp_set = False

            if num_comp_set is False:
                sample.set_numComponents(num_comp_in)
                for smol_in in sup_mols_in:
                    smol_type_in = smol_in.original_tagname_
                    if smol_type_in == 'sample_supramolecule':
                        self.check_set(smol_in.get_number_unique_components, sample.set_numComponents)
                        num_comp_set = True
                        # num_comp_in -= 1
            if num_comp_set is False:
                sample.set_numComponents(num_comp_in)
            # element 2 - <xs:complexType name="samplType">
            # XSD: <xs:element name="name" type="xs:string"/>
            smpl_name = sample_in.get_name()
            if smpl_name is not None:
                sample.set_name(smpl_name.get_valueOf_())
                # Override previously set sample name if empty
                if not sample.get_name():
                    if num_comp_in > 0:
                        for smol_in in sup_mols_in:
                            smol_parent = smol_in.get_parent()
                            if smol_parent is None or smol_parent == 0:
                                sample_name_in = smol_in.get_name()
                                if sample_name_in is not None:
                                    sample.set_name(sample_name_in.get_valueOf_())
            elif not self.roundtrip:
                sample.set_name('')

            # element 3 - <xs:complexType name="samplType">
            # XSD: <xs:element name="compDegree" type="xs:string" minOccurs="0"/>
            if num_comp_in > 0:
                for smol_in in sup_mols_in:
                    smol_type_in = smol_in.original_tagname_
                    if smol_type_in == 'sample_supramolecule':
                        self.check_set(smol_in.get_oligomeric_state, sample.set_compDegree)
                    else:
                        smol_parent = smol_in.get_parent()
                        if smol_parent is not None or smol_parent == 0:
                            self.check_set(smol_in.get_oligomeric_state, sample.set_compDegree)

            # element 4 - <xs:complexType name="samplType">
            # XSD: <xs:element name="molWtTheo" type="mwType" minOccurs="0"/>
            for smol_in in sup_mols_in:
                smol_type_in = smol_in.original_tagname_
                if smol_type_in == 'sample_supramolecule':
                    weight = smol_in.get_molecular_weight()
                    set_mol_weight(sample, weight, meth=True)

            # element 5 - <xs:complexType name="samplType">
            # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
            if num_comp_in > 0:
                for smol_in in sup_mols_in:
                    smol_type_in = smol_in.original_tagname_
                    if smol_type_in == 'sample_supramolecule':
                        self.check_set(smol_in.get_details, sample.set_details)
                    else:
                        smol_parent = smol_in.get_parent()
                        if smol_parent == 0:
                            self.check_set(smol_in.get_details, sample.set_details)

            # element 6 - <xs:complexType name="samplType">
            # XSD: <xs:element name="molWtMethod" type="xs:string" minOccurs="0"/>

            # element 7 - <xs:complexType name="samplType">
            # XSD: <xs:element name="molWtExp" type="mwType" minOccurs="0"/>

            # element 8 - <xs:complexType name="samplType">
            # XSD: <xs:element name="sampleComponentList" type="smplCompListType"/>
            if num_comp_in > 0:
                # XSD: <xs:complexType name="smplCompListType"> has 1 element
                comp_list = emdb_19.smplCompListType()
                comp_id = 1
                for smol_in in sup_mols_in:
                    smol_type_in = smol_in.original_tagname_
                    if smol_type_in != 'sample_supramolecule':
                        # element 1 - <xs:complexType name="smplCompListType">
                        # XSD: <xs:element name="sampleComponent" type="smplCompType" maxOccurs="unbounded"/>
                        # XSD: <xs:complexType name="smplCompType"> has 6 elements, 1 attribute and 1 choice of 8 elements
                        comp = emdb_19.smplCompType()
                        # attribute 1 - <xs:complexType name="smplCompType">
                        # XSD: <xs:attribute name="componentID" type="xs:positiveInteger" use="required"/>
                        comp.set_componentID(comp_id)
                        comp_id += 1
                        # element 1 - <xs:complexType name="smplCompType">
                        # XSD: <xs:element name="entry" type="cmpntClassType"/>
                        if smol_type_in == 'virus_supramolecule':
                            comp.set_entry('virus')
                        elif smol_type_in in ['organelle_or_cellular_component_supramolecule', 'cell_supramolecule', 'tissue_supramolecule']:
                            comp.set_entry('cellular-component')
                        elif smol_type_in == 'complex_supramolecule':
                            rib_detail = smol_in.get_ribosome_details()
                            if rib_detail is not None:
                                if rib_detail.find('eukaryo') != -1:
                                    comp.set_entry('ribosome-eukaryote')
                                elif rib_detail.find('prokaryo') != -1:
                                    comp.set_entry('ribosome-prokaryote')
                                else:
                                    comp.set_entry('protein')
                            else:
                                comp.set_entry('protein')

                        sci_name = smol_in.get_name()
                        if sci_name is not None:
                            # element 2 - <xs:complexType name="smplCompType">
                            # XSD: <xs:element name="sciName" type="xs:string"/>
                            name = sci_name.get_valueOf_()
                            comp.set_sciName(name)
                            # element 3 - <xs:complexType name="smplCompType">
                            # XSD: <xs:element name="synName" type="xs:string" minOccurs="0"/>
                            syn = sci_name.get_synonym()
                            if syn is not None:
                                comp.set_synName(syn)
                        # element 4 - <xs:complexType name="smplCompType">
                        # XSD: <xs:element name="molWtTheo" type="mwType" minOccurs="0"/>
                        # element 5 - <xs:complexType name="smplCompType">
                        # XSD: <xs:element name="molWtExp" type="mwType" minOccurs="0"/>
                        if smol_parent is None or smol_parent == 0:
                            # IL 1/Mar/2015 - tissue and cell do not have get_molecular_weight method()
                            if smol_type_in not in ['tissue_supramolecule', 'cell_supramolecule']:
                                weight = smol_in.get_molecular_weight()
                                #set_mol_weight(sample, weight, meth=False)
                        if smol_type_in in ['virus_supramolecule', 'organelle_or_cellular_component_supramolecule',
                                            'complex_supramolecule']:
                            weight = smol_in.get_molecular_weight()
                            set_mol_weight(comp, weight, meth=False)
                        # element 6 - <xs:complexType name="smplCompType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        if smol_parent == 0:
                            self.check_set(smol_in.get_details, sample.set_details)
                        if smol_type_in in ['sample_supramolecule', 'virus_supramolecule']:
                            self.check_set(smol_in.get_details, comp.set_details)
#                         if smol_type_in in ['organelle_or_cellular_component_supramolecule',
#                                             'cell_supramolecule', 'complex_supramolecule']:
#                             unpack_odd_details(smol_in, comp, protein)

                    # choice 1 - <xs:complexType name="smplCompType"> of 8 elements
                    # element 1 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="protein" type="proteinType"/>
                    if smol_type_in == 'complex_supramolecule':
                        rib_detail = smol_in.get_ribosome_details()
                        not_euk = True
                        not_pro = True
                        if rib_detail is not None:
                            not_euk = rib_detail.find('eukaryo') == -1
                            not_pro = rib_detail.find('prokaryo') == -1

                        if not_euk and not_pro:
                            # XSD: <xs:complexType name="proteinType"> has 10 elements
                            protein = emdb_19.proteinType()
                            # element 1 - <xs:complexType name="proteinType">
                            # XSD: <xs:element name="sciSpeciesName" type="sciSpeciesType" minOccurs="0"/>
                            # element 2 - <xs:complexType name="proteinType">
                            # XSD: <xs:element name="sciSpeciesStrain" type="xs:string" minOccurs="0"/>
                            # element 3 - <xs:complexType name="proteinType">
                            # XSD: <xs:element name="synSpeciesName" type="xs:string" minOccurs="0"/>
                            # element 4 - <xs:complexType name="proteinType">
                            # XSD: <xs:element name="oligomericDetails" type="xs:string" minOccurs="0"/>
                            self.check_set(smol_in.get_oligomeric_state, protein.set_oligomericDetails)
                            # element 5 - <xs:complexType name="proteinType">
                            # XSD: <xs:element name="numCopies" type="xs:string" minOccurs="0"/>
                            self.check_set(smol_in.get_number_of_copies, protein.set_numCopies)
                            # element 6 - <xs:complexType name="proteinType">
                            # XSD: <xs:element name="recombinantExpFlag" type="xs:boolean"/>
                            smol_rec_flag = smol_in.get_recombinant_exp_flag() or False
                            protein.set_recombinantExpFlag(smol_rec_flag)
                            # copy_recombinant_source(smol_in, protein)
                            # element 7 - <xs:complexType name="proteinType">
                            # XSD: <xs:element name="natSource" type="natSrcType" minOccurs="0"/>
                            if self.relaxed:
                                copy_natural_source(smol_in, protein)
                            else:
                                copy_natural_source(smol_in, protein, cell=False, organelle=False, tissue=False, cellular_location=False, organ=False)
                            # element 8 - <xs:complexType name="proteinType">
                            # XSD: <xs:element name="engSource" type="engSrcType" minOccurs="0"/>
                            eng_src = create_eng_source(smol_in)
                            if eng_src is not None and eng_src.hasContent_:
                                protein.set_engSource(eng_src)
                            # element 9 - <xs:complexType name="proteinType">
                            # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            unpack_odd_details(smol_in, comp, protein)
                            # element 10 - <xs:complexType name="proteinType">
                            # XSD: <xs:element name="externalReferences" type="externalReferencesType" minOccurs="0"/>
                            copy_external_references(smol_in.get_external_references, protein.set_externalReferences)

                            comp.set_protein(protein)
                    # element 2 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="cellular-component" type="cellCompType"/>
                    if smol_type_in in ['organelle_or_cellular_component_supramolecule', 'cell_supramolecule', 'tissue_supramolecule']:
                        # Treat this as a cellular component as there is no better mapping
                        # XSD: <xs:complexType name="cellCompType"> has 10 elements
                        cell = emdb_19.cellCompType()
                        # element 1 - <xs:complexType name="cellCompType">
                        # XSD: <xs:element name="sciSpeciesName" type="sciSpeciesType" minOccurs="0"/>
                        # element 2 - <xs:complexType name="cellCompType">
                        # XSD: <xs:element name="sciSpeciesStrain" type="xs:string" minOccurs="0"/>
                        # element 3 - <xs:complexType name="cellCompType">
                        # XSD: <xs:element name="synSpeciesName" type="xs:string" minOccurs="0"/>
                        # element 4 - <xs:complexType name="cellCompType">
                        # XSD: <xs:element name="oligomericDetails" type="xs:string" minOccurs="0"/>
                        self.check_set(smol_in.get_oligomeric_state, cell.set_oligomericDetails)
                        # element 5 - <xs:complexType name="cellCompType">
                        # XSD: <xs:element name="numCopies" type="xs:string" minOccurs="0"/>
                        self.check_set(smol_in.get_number_of_copies, cell.set_numCopies)
                        # element 6 - <xs:complexType name="cellCompType">
                        # XSD: <xs:element name="recombinantExpFlag" type="xs:boolean"/>
                        smol_rec_flag = smol_in.get_recombinant_exp_flag() or False
                        cell.set_recombinantExpFlag(smol_rec_flag)
                        # element 7 - <xs:complexType name="cellCompType">
                        # XSD: <xs:element name="natSource" type="natSrcType" minOccurs="0"/>
                        if smol_type_in == 'organelle_or_cellular_component_supramolecule':
                            copy_natural_source(smol_in, cell)
                        if smol_type_in == 'cell_supramolecule':
                            copy_natural_source(smol_in, cell, organelle=False, cellular_location=False)
                        # element 8 - <xs:complexType name="cellCompType">
                        # XSD: <xs:element name="engSource" type="engSrcType" minOccurs="0"/>
                        # XSD: <xs:complexType name="engSrcType"> has 4 elements
                        if smol_type_in == 'organelle_or_cellular_component_supramolecule':
                            eng_src = create_eng_source(smol_in)
                            if eng_src is not None and eng_src.hasContent_:
                                cell.set_engSource(eng_src)
                        # element 9 - <xs:complexType name="cellCompType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        self.check_set(smol_in.get_details, comp.set_details)
                        # element 10 - <xs:complexType name="cellCompType">
                        # XSD: <xs:element name="externalReferences" type="externalReferencesType" minOccurs="0"/>
                        copy_external_references(smol_in.get_external_references, cell.set_externalReferences)

                        comp.set_cellular_component(cell)

                    # element 3 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="virus" type="virusType"/>
                    if smol_type_in == 'virus_supramolecule':
                        # XSD: <xs:complexType name="virusType"> has unbound number of 14 choices ????!!!
                        vir = emdb_19.virusType()
                        unpack_odd_details(smol_in, comp, vir)
                        # element 1 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="sciSpeciesName" type="sciSpeciesType" minOccurs="1" maxOccurs="1"/>
                        ang_in_details = smol_in.get_sci_species_name()
                        if ang_in_details is not None:
                            vir.set_sciSpeciesName(emdb_19.sciSpeciesType(valueOf_=ang_in_details.get_valueOf_(), ncbiTaxId=ang_in_details.get_ncbi()))
                        # element 2 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="synSpeciesName" type="xs:string" minOccurs="0" maxOccurs="1"/>
                        self.check_set(smol_in.get_syn_species_name, vir.set_synSpeciesName)
                        # element 3 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="sciSpeciesSerotype" type="xs:string" minOccurs="0" maxOccurs="1"/>
                        self.check_set(smol_in.get_sci_species_serotype, vir.set_sciSpeciesSerotype)
                        # element 4 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="sciSpeciesSerocomplex" type="xs:string" minOccurs="0" maxOccurs="1"/>
                        self.check_set(smol_in.get_sci_species_serocomplex, vir.set_sciSpeciesSerocomplex)
                        # element 5 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="sciSpeciesSubspecies" type="xs:string" minOccurs="0" maxOccurs="1"/>
                        self.check_set(smol_in.get_sci_species_subspecies, vir.set_sciSpeciesSubspecies)
                        # element 6 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="sciSpeciesStrain" type="xs:string" minOccurs="0" maxOccurs="1"/>
                        self.check_set(smol_in.get_sci_species_strain, vir.set_sciSpeciesStrain)
                        # element 7 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="empty" type="xs:boolean" minOccurs="1" maxOccurs="1"/>
                        vir.set_empty(smol_in.get_virus_empty())
                        # element 8 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="enveloped" type="xs:boolean" minOccurs="1" maxOccurs="1"/>
                        vir.set_enveloped(smol_in.get_virus_enveloped())
                        # element 9 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="isolate" type="virusIsolType" minOccurs="1" maxOccurs="1"/>
                        vir.set_isolate(smol_in.get_virus_isolate())
                        # element 10 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="class" type="virusClassType" minOccurs="1" maxOccurs="1"/>
                        vir.set_class(smol_in.get_virus_type())
                        # element 11 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="externalReferences" type="externalReferencesType" minOccurs="0"/>
                        copy_external_references(smol_in.get_external_references, vir.set_externalReferences)
                        # element 12 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="natSource" type="natSrcVirusType" minOccurs="0"/>
                        # XSD: <xs:complexType name="natSrcVirusType"> has 3 elements
                        ns_in = smol_in.get_natural_host()
                        if ns_in is not None and len(ns_in) > 0:
                            n_in = ns_in[0]
                            nat_src_virus = emdb_19.natSrcVirusType()
                            # element 1 - <xs:complexType name="natSrcVirusType">
                            # XSD: <xs:element name="hostCategory" type="hostCategoryType" minOccurs="0" maxOccurs="1"/>
                            self.check_set(n_in.get_synonym_organism, nat_src_virus.set_hostCategory)
                            # element 2 - <xs:complexType name="natSrcVirusType">
                            # XSD: <xs:element name="hostSpecies" type="sciSpeciesType" minOccurs="0" maxOccurs="1"/>
                            org_in = n_in.get_organism()
                            if org_in is not None:
                                # XSD: <xs:complexType name="sciSpeciesType"> extension of string and has 1 attribute
                                sci_species = emdb_19.sciSpeciesType()
                                # XSD: <xs:extension base="xs:string">
                                sci_species.set_valueOf_(org_in.get_valueOf_())
                                # attribute 1 - <xs:complexType name="sciSpeciesType">
                                # XSD: <xs:attribute name="ncbiTaxId" type="xs:integer"/>
                                sci_species.set_ncbiTaxId(org_in.get_ncbi())
                                if sci_species.hasContent_():
                                    nat_src_virus.set_hostSpecies(sci_species)
                            # element 3 - <xs:complexType name="natSrcVirusType">
                            # XSD: <xs:element name="hostSpeciesStrain" type="xs:string" minOccurs="0" maxOccurs="1"/>
                            ang_in_details = n_in.get_strain()
                            if ang_in_details is not None:
                                nat_src_virus.set_hostSpeciesStrain(ang_in_details)

                            vir.add_natSource(nat_src_virus)

                        # element 13 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="engSource" type="engSrcType" minOccurs="0"/>
                        hs_in = smol_in.get_host_system()
                        if hs_in is not None:
                            # XSD: <xs:complexType name="engSrcType"> has 4 elements
                            esrc = emdb_19.engSrcType()
                            # element 1 - <xs:complexType name="engSrcType">
                            # XSD: <xs:element name="expSystem" type="sciSpeciesType" minOccurs="0"/>
                            org_in = hs_in.get_recombinant_organism()
                            if org_in is not None:
                                # XSD: <xs:complexType name="sciSpeciesType"> base of string and has 1 attribute
                                species = emdb_19.sciSpeciesType()
                                # XSD: <xs:extension base="xs:string">
                                self.check_set(org_in.get_valueOf_, species.set_valueOf_)
                                # attribute 1 - <xs:complexType name="sciSpeciesType">
                                # XSD: <xs:attribute name="ncbiTaxId" type="xs:integer"/>
                                self.check_set(org_in.get_ncbi, species.set_ncbiTaxId)

                                esrc.set_expSystem(species)
                            # element 2 - <xs:complexType name="engSrcType">
                            # XSD: <xs:element name="expSystemStrain" type="xs:string" minOccurs="0"/>
                            self.check_set(hs_in.get_recombinant_strain, esrc.set_expSystemStrain)
                            # element 3 - <xs:complexType name="engSrcType">
                            # XSD: <xs:element name="expSystemCell" type="xs:string" minOccurs="0"/>
                            self.check_set(hs_in.get_recombinant_cell, esrc.set_expSystemCell)
                            # element 4 - <xs:complexType name="engSrcType">
                            # XSD: <xs:element name="vector" type="xs:string" minOccurs="0"/>
                            self.check_set(hs_in.get_recombinant_plasmid, esrc.set_vector)

                            vir.add_engSource(esrc)

                        # element 14 - <xs:complexType name="virusType">
                        # XSD: <xs:element name="shell" type="shellType" minOccurs="0"/>
                        shell_list_in = smol_in.get_virus_shell()
                        for shell_in in shell_list_in:
                            # XSD <xs:complexType name="shellType"> has 3 element and 1 attribute
                            shell = emdb_19.shellType()
                            # element 1 - <xs:complexType name="shellType">
                            # XSD: <xs:element name="nameElement" type="xs:string" minOccurs="0" maxOccurs="1"/>
                            self.check_set(shell_in.get_name, shell.set_nameElement)
                            # element 2 - <xs:complexType name="shellType">
                            # XSD: <xs:element name="diameter" type="diamType" minOccurs="0" maxOccurs="1"/>
                            self.set_value_and_units(shell_in.get_diameter, shell.set_diameter, emdb_19.diamType, units='A')
                            # element 3 - <xs:complexType name="shellType">
                            # XSD: <xs:element name="tNumber" type="floatOrNAType" minOccurs="0" maxOccurs="1"/>
                            self.check_set(shell_in.get_triangulation, shell.set_tNumber)
                            # attribute 1 - <xs:complexType name="shellType">
                            # XSD: <xs:attribute name="id" type="xs:positiveInteger" use="required"/>
                            self.check_set(shell_in.get_shell_id, shell.set_id)

                            vir.add_shell(shell)

                        comp.set_virus(vir)

                    # element 4 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="nucleic-acid" type="nuclAcidType"/>
                    # element 5 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="ligand" type="ligandType"/>
                    # element 6 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="label" type="labelType"/>

                    # elements 7 and 8
                    if smol_type_in == 'complex_supramolecule':
                        rib_detail = smol_in.get_ribosome_details()
                        if rib_detail is not None:
                            rib = None
                            if rib_detail.find('eukaryo') != -1:
                                # element 7 in choice 1 - <xs:complexType name="smplCompType">
                                # XSD: <xs:element name="ribosome-eukaryote" type="riboTypeEu"/>
                                # XSD: <xs:complexType name="riboTypeEu"> has 11 elements
                                rib = emdb_19.riboTypeEu()
                                # element 1 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="eukaryote" type="xs:string" minOccurs="1"/>
                                rib_details = smol_in.get_ribosome_details()
                                if rib_details.find('ribosome-eukaryote:') != -1:
                                    rib_details = rib_details.replace('ribosome-eukaryote: ', '')
                                rib.set_eukaryote(rib_details)
                                # element 2 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="sciSpeciesName" type="sciSpeciesType" minOccurs="0" maxOccurs="1"/>
                                self.check_set(smol_in.get_name, rib.set_sciSpeciesName)
                                # element 3 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="sciSpeciesStrain" type="xs:string" minOccurs="0"/>
                                # element 4 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="synSpeciesName" type="xs:string" minOccurs="0" maxOccurs="1"/>
                                # element 5 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="oligomericDetails" type="xs:string" minOccurs="0" maxOccurs="1"/>
                                # element 6 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="numCopies" type="xs:string" minOccurs="0" maxOccurs="1"/>
                                # element 7 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="recombinantExpFlag" type="xs:boolean" minOccurs="0" maxOccurs="1"/>
                                self.check_set(smol_in.get_recombinant_exp_flag, rib.set_recombinantExpFlag)
                                # copy_recombinant_source(smol_in, rib)
                                # element 8 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="natSource" type="natSrcType" minOccurs="0" maxOccurs="1"/>
                                copy_natural_source(smol_in, rib, cell=False, organelle=False, tissue=False, cellular_location=False, organ=False)
                                # element 9 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="engSource" type="engSrcType" minOccurs="0" maxOccurs="1"/>
                                eng_src = create_eng_source(smol_in)
                                if eng_src is not None and eng_src.hasContent_:
                                    rib.set_engSource(eng_src)
                                # element 10 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="details" type="xs:string" minOccurs="0" maxOccurs="1"/>
                                unpack_odd_details(smol_in, comp, rib)
                                # element 11 - <xs:complexType name="riboTypeEu">
                                # XSD: <xs:element name="externalReferences" type="externalReferencesType" minOccurs="0" maxOccurs="1"/>
                                copy_external_references(smol_in.get_external_references, rib.set_externalReferences)

                                comp.set_ribosome_eukaryote(rib)

                            if rib_detail.find('prokaryo') != -1:
                                # element 8 in choice 1 - <xs:complexType name="smplCompType">
                                # XSD: <xs:element name="ribosome-prokaryote" type="riboTypePro"/>
                                # XSD: <xs:complexType name="riboTypePro"> has 11 elements
                                rib = emdb_19.riboTypePro()
                                # element 1 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="prokaryote" type="xs:string" minOccurs="1"/>
                                rib_details = smol_in.get_ribosome_details()
                                if rib_details.find('ribosome-prokaryote:') != -1:
                                    rib_details = rib_details.replace('ribosome-prokaryote: ', '')
                                rib.set_prokaryote(rib_details)
                                # element 2 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="sciSpeciesName" type="sciSpeciesType" minOccurs="0" maxOccurs="1"/>
                                self.check_set(smol_in.get_name, rib.set_sciSpeciesName)
                                # element 3 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="sciSpeciesStrain" type="xs:string" minOccurs="0"/>
                                # element 4 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="synSpeciesName" type="xs:string" minOccurs="0" maxOccurs="1"/>
                                # element 5 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="oligomericDetails" type="xs:string" minOccurs="0" maxOccurs="1"/>
                                # element 6 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="numCopies" type="xs:string" minOccurs="0" maxOccurs="1"/>
                                # element 7 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="recombinantExpFlag" type="xs:boolean" minOccurs="0" maxOccurs="1"/>
                                self.check_set(smol_in.get_recombinant_exp_flag, rib.set_recombinantExpFlag)
                                # copy_recombinant_source(smol_in, rib)
                                # element 8 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="natSource" type="natSrcType" minOccurs="0" maxOccurs="1"/>
                                copy_natural_source(smol_in, rib)#, cell=False, organelle=False, tissue=False, cellular_location=False, organ=False)
                                # element 9 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="engSource" type="engSrcType" minOccurs="0" maxOccurs="1"/>
                                eng_src = create_eng_source(smol_in)
                                if eng_src is not None and eng_src.hasContent_:
                                    rib.set_engSource(eng_src)
                                # element 10 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="details" type="xs:string" minOccurs="0" maxOccurs="1"/>
                                unpack_odd_details(smol_in, comp, rib)
                                # element 11 - <xs:complexType name="riboTypePro">
                                # XSD: <xs:element name="externalReferences" type="externalReferencesType" minOccurs="0" maxOccurs="1"/>
                                copy_external_references(smol_in.get_external_references, rib.set_externalReferences)

                                comp.set_ribosome_prokaryote(rib)
                    if smol_type_in != 'sample_supramolecule':
                        comp_list.add_sampleComponent(comp)
                for mol_in in mols_in:
                    mol_type_in = mol_in.original_tagname_
                    # element 1 - <xs:complexType name="smplCompListType">
                    # XSD: <xs:element name="sampleComponent" type="smplCompType" maxOccurs="unbounded"/>
                    # XSD: <xs:complexType name="smplCompType"> has 6 elements, 1 attribute and 1 choice of 8 elements
                    comp = emdb_19.smplCompType()
                    # attribute 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:attribute name="componentID" type="xs:positiveInteger" use="required"/>
                    comp.set_componentID(comp_id)
                    comp_id += 1
                    other_mol_nuc_acid = False
                    # element 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="entry" type="cmpntClassType"/>
                    if mol_type_in == 'protein_or_peptide':
                        comp.set_entry('protein')
                    elif mol_type_in == 'ligand':
                        comp.set_entry('ligand')
                    elif mol_type_in == 'em_label':
                        comp.set_entry('label')
                    elif mol_type_in in ['dna', 'rna']:
                        comp.set_entry('nucleic-acid')
                    elif mol_type_in == 'other_macromolecule':
                        mol_class_in = mol_in.get_classification()
                        if mol_class_in is not None:
                            if mol_class_in in ['DNA/RNA', 'OTHER_NA', 'other', 'polydeoxyribonucleotide/polyribonucleotide hybrid']:
                                comp.set_entry('nucleic-acid')
                                other_mol_nuc_acid = True

                    sci_name = mol_in.get_name()
                    if sci_name is not None:
                        # element 2 - <xs:complexType name="smplCompType">
                        # XSD: <xs:element name="sciName" type="xs:string"/>
                        name = sci_name.get_valueOf_()
                        comp.set_sciName(name)
                        # element 3 - <xs:complexType name="smplCompType">
                        # XSD: <xs:element name="synName" type="xs:string" minOccurs="0"/>
                        syn = sci_name.get_synonym()
                        if syn is not None:
                            comp.set_synName(syn)
                    # element 4 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="molWtTheo" type="mwType" minOccurs="0"/>
                    # element 5 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="molWtExp" type="mwType" minOccurs="0"/>
                    set_mol_weight(comp, mol_in.get_molecular_weight(), meth=False)
                    # element 6 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>

                    # choice 1 - <xs:complexType name="smplCompType"> of 8 elements
                    # element 1 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="protein" type="proteinType"/>
                    if mol_type_in == 'protein_or_peptide':
                        # XSD: <xs:complexType name="proteinType"> has 10 elements
                        protein = emdb_19.proteinType()
                        unpack_odd_details(mol_in, comp, protein)
                        # element 1 - <xs:complexType name="proteinType">
                        # XSD: <xs:element name="sciSpeciesName" type="sciSpeciesType" minOccurs="0"/>
                        # element 2 - <xs:complexType name="proteinType">
                        # XSD: <xs:element name="sciSpeciesStrain" type="xs:string" minOccurs="0"/>
                        # element 3 - <xs:complexType name="proteinType">
                        # XSD: <xs:element name="synSpeciesName" type="xs:string" minOccurs="0"/>
                        # element 4 - <xs:complexType name="proteinType">
                        # XSD: <xs:element name="oligomericDetails" type="xs:string" minOccurs="0"/>
                        self.check_set(mol_in.get_oligomeric_state, protein.set_oligomericDetails)
                        # element 5 - <xs:complexType name="proteinType">
                        # XSD: <xs:element name="numCopies" type="xs:string" minOccurs="0"/>
                        self.check_set(mol_in.get_number_of_copies, protein.set_numCopies)
                        # element 6 - <xs:complexType name="proteinType">
                        # XSD: <xs:element name="recombinantExpFlag" type="xs:boolean"/>
                        mol_rec_flag = mol_in.get_recombinant_exp_flag() or False
                        protein.set_recombinantExpFlag(mol_rec_flag)
                        # copy_recombinant_source(smol_in, protein)
                        # element 7 - <xs:complexType name="proteinType">
                        # XSD: <xs:element name="natSource" type="natSrcType" minOccurs="0"/>
                        copy_natural_source(mol_in, protein)#, cell=False, organelle=False, tissue=False, cellular_location=False, organ=False)
                        # element 8 - <xs:complexType name="proteinType">
                        # XSD: <xs:element name="engSource" type="engSrcType" minOccurs="0"/>
                        eng_src = create_eng_source(mol_in)
                        if eng_src is not None and eng_src.hasContent_:
                            protein.set_engSource(eng_src)
                        # element 9 - <xs:complexType name="proteinType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        unpack_odd_details(mol_in, comp, protein)
                        # element 10 - <xs:complexType name="proteinType">
                        # XSD: <xs:element name="externalReferences" type="externalReferencesType" minOccurs="0"/>
                        seq = mol_in.get_sequence()
                        copy_external_references(seq.get_external_references, protein.set_externalReferences)

                        comp.set_protein(protein)
                    # element 2 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="cellular-component" type="cellCompType"/>
                    # element 3 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="virus" type="virusType"/>
                    # element 4 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="nucleic-acid" type="nuclAcidType"/>

                    if mol_type_in in ['dna', 'rna'] or other_mol_nuc_acid:
                        # XSD: <xs:complexType name="nuclAcidType"> has 7 elements
                        nuc_acid = emdb_19.nuclAcidType()
                        unpack_odd_details(mol_in, comp, nuc_acid)
                        copy_natural_source(mol_in, nuc_acid, cell=False, organelle=False, tissue=False, cellular_location=False, organ=False, only_common=True)
                        # element 1 - <xs:complexType name="nuclAcidType">
                        # XSD: <xs:element name="sciSpeciesName" type="sciSpeciesType" minOccurs="0"/>
                        # element 2 - <xs:complexType name="nuclAcidType">
                        # XSD: <xs:element name="sciSpeciesStrain" type="xs:string" minOccurs="0"/>
                        # element 3 - <xs:complexType name="nuclAcidType">
                        # XSD: <xs:element name="synSpeciesName" type="xs:string" minOccurs="0"/>
                        # element 4 - <xs:complexType name="nuclAcidType">
                        # XSD: <xs:element name="syntheticFlag" type="xs:boolean"/>
                        na_synth_in = mol_in.get_synthetic_flag() or False
                        nuc_acid.set_syntheticFlag(na_synth_in)
                        # element 5 - <xs:complexType name="nuclAcidType">
                        # XSD: <xs:element name="sequence" type="xs:string" minOccurs="0"/>
                        seq_in = mol_in.get_sequence()
                        if seq_in is not None:
                            nuc_acid.set_sequence(seq_in.get_string())
                        # element 6 - <xs:complexType name="nuclAcidType">
                        # XSD: <xs:element name="class" type="naClassType"/>
                        if mol_type_in == 'rna':
                            na_class_in = mol_in.get_classification()
                            if na_class_in == 'TRANSFER':
                                nuc_acid.set_class('T-RNA')
                            else:
                                nuc_acid.set_class('RNA')
                        elif mol_type_in == 'dna':
                            nuc_acid.set_class('DNA')
                        else:
                            if other_mol_nuc_acid:
                                mol_class_in = mol_in.get_classification()
                                if mol_class_in is not None:
                                    if mol_class_in in ['DNA/RNA', 'polydeoxyribonucleotide/polyribonucleotide hybrid']:
                                        nuc_acid.set_class('DNA/RNA')
                                    if mol_class_in in ['OTHER_NA', 'other']:
                                        nuc_acid.set_class('OTHER')
                        # element 7 - <xs:complexType name="nuclAcidType">
                        # XSD: <xs:element name="structure" type="naStructType"/>
                        na_struct_in = mol_in.get_structure() or 'OTHER'
                        nuc_acid.set_structure(na_struct_in)

                        comp.set_nucleic_acid(nuc_acid)

                    # element 5 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="ligand" type="ligandType"/>
                    if mol_type_in == 'ligand':
                        # XSD: <xs:complexType name="ligandType"> has 10 elements
                        lig = emdb_19.ligandType()
                        unpack_odd_details(mol_in, comp, lig)
                        # element 1 - <xs:complexType name="ligandType">
                        # XSD: <xs:element name="sciSpeciesName" type="sciSpeciesType" minOccurs="0"/>
                        # element 2 - <xs:complexType name="ligandType">
                        # XSD: <xs:element name="sciSpeciesStrain" type="xs:string" minOccurs="0"/>
                        # element 3 - <xs:complexType name="ligandType">
                        # XSD: <xs:element name="synSpeciesName" type="xs:string" minOccurs="0"/>
                        # element 4 - <xs:complexType name="ligandType">
                        # XSD: <xs:element name="oligomericDetails" type="xs:string" minOccurs="0"/>
                        self.check_set(mol_in.get_oligomeric_state, lig.set_oligomericDetails)
                        # element 5 - <xs:complexType name="ligandType">
                        # XSD: <xs:element name="numCopies" type="xs:string" minOccurs="0"/>
                        self.check_set(mol_in.get_number_of_copies, lig.set_numCopies)
                        # element 6 - <xs:complexType name="ligandType">
                        # XSD: <xs:element name="recombinantExpFlag" type="xs:boolean"/>
                        ligant_rec_flag = mol_in.get_recombinant_exp_flag() or False
                        lig.set_recombinantExpFlag(ligant_rec_flag)
                        # copy_recombinant_source(mol_in, lig)
                        # element 7 - <xs:complexType name="ligandType">
                        # XSD: <xs:element name="natSource" type="natSrcType" minOccurs="0"/>
                        copy_natural_source(mol_in, lig)
                        # element 8 - <xs:complexType name="ligandType">
                        # XSD: <xs:element name="engSource" type="engSrcType" minOccurs="0"/>
                        eng_src = create_eng_source(mol_in)
                        if eng_src is not None and eng_src.hasContent_:
                            lig.set_engSource(eng_src)
                        # element 9 - <xs:complexType name="ligandType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        # element 10 - <xs:complexType name="ligandType">
                        # XSD: <xs:element name="externalReferences" type="externalReferencesType" minOccurs="0"/>
                        copy_external_references(mol_in.get_external_references, lig.set_externalReferences)

                        comp.set_ligand(lig)

                    # element 6 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="label" type="labelType"/>
                    if mol_type_in == 'em_label':
                        # XSD: <xs:complexType name="labelType"> has 3 elements
                        lab = emdb_19.labelType()
                        unpack_odd_details(mol_in, comp, lab)
                        # element 1 - <xs:complexType name="labelType">
                        # XSD: <xs:element name="formula" type="xs:string" minOccurs="0"/>
                        self.check_set(mol_in.get_formula, lab.set_formula)
                        # element 2 - <xs:complexType name="labelType">
                        # XSD: <xs:element name="oligomericDetails" type="xs:string" minOccurs="0"/>
                        self.check_set(mol_in.get_oligomeric_state, lab.set_oligomericDetails)
                        # element 3 - <xs:complexType name="labelType">
                        # XSD: <xs:element name="numCopies" type="xs:string" minOccurs="0"/>
                        self.check_set(mol_in.get_number_of_copies, lab.set_numCopies)

                        comp.set_label(lab)

                    # element 7 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="ribosome-eukaryote" type="riboTypeEu"/>
                    # element 8 in choice 1 - <xs:complexType name="smplCompType">
                    # XSD: <xs:element name="ribosome-prokaryote" type="riboTypePro"/>

                    comp_list.add_sampleComponent(comp)

                if comp_list.hasContent_():
                    sample.set_sampleComponentList(comp_list)

        # xml_out.set_sample(sample)

        # element 6 - <xs:complexType name="entryType">
        # XSD: <xs:element name="experiment" type="expType" maxOccurs="1"/>
        # Assume that this element exists!
        sd_in = None
        if xml_in is not None:
            sd_in = xml_in.get_structure_determination_list().get_structure_determination()[0]
        if sd_in is not None:
            # XSD <xs:complexType name="expType"> has 5 elements
            exp = emdb_19.expType()
            # element 1 - <xs:complexType name="expType">
            # XSD: <xs:element name="vitrification" type="vitrifType" maxOccurs="unbounded"/>
            # in 3.0 vitrification and specimen preparation are combined - therefore deal with them at the same time
            sd_in_specimen_preparation_list = sd_in.get_specimen_preparation_list()
            if sd_in_specimen_preparation_list is not None:
                spec_prep_list_in = sd_in_specimen_preparation_list.get_specimen_preparation()
                for sp_in in spec_prep_list_in:
                    # XSD: <xs:complexType name="vitrifType"> has 7 elements
                    vitr = emdb_19.vitrifType()
                    vitr_in = sp_in.get_vitrification()
                    if vitr_in is not None:
                        # element 1 - <xs:complexType name="vitrifType">
                        # XSD: <xs:element name="cryogenName" type="cryogenType"/>
                        cryogen_out = vitr_in.get_cryogen_name()
                        if cryogen_out is not None:
                            if self.roundtrip:
                                vitr.set_cryogenName(cryogen_out)
                            else:
                                allowed_cryo_names = ['ETHANE', 'ETHANE-PROPANE MIXTURE', 'METHANE', 'NITROGEN', 'HELIUM', 'PROPANE', 'FREON 12', 'FREON 22', 'NONE', 'OTHER']
                                if cryogen_out in allowed_cryo_names:
                                    vitr.set_cryogenName(cryogen_out)
                                else:
                                    vitr.set_cryogenName('OTHER')
                        # element 2 - <xs:complexType name="vitrifType">
                        # XSD: <xs:element name="humidity" type="xs:string" minOccurs="0"/>
                        chamber_humidity = vitr_in.get_chamber_humidity()
                        if chamber_humidity is not None:
                            hum = chamber_humidity.get_valueOf_()
                            vitr.set_humidity(hum)
                        # element 3 - <xs:complexType name="vitrifType">
                        # XSD: <xs:element name="temperature" type="tempType" minOccurs="0"/>
                        ang_in_details = vitr_in.get_chamber_temperature()
                        if ang_in_details is not None:
                            vitr.set_temperature(emdb_19.tempType(valueOf_=ang_in_details.get_valueOf_(), units=const.U_KELF))
                        # element 4 - <xs:complexType name="vitrifType">
                        # XSD: <xs:element name="instrument" type="vitrInstrType" minOccurs="0"/>
                        vitrification_instrument = vitr_in.get_instrument()
                        # allowed_vitrification_instruments = ['FEI VITROBOT MARK I', 'FEI VITROBOT MARK II', 'FEI VITROBOT MARK III', 'FEI VITROBOT MARK IV', 'GATAN CRYOPLUNGE 3', 'HOMEMADE PLUNGER', 'LEICA EM CPC', 'LEICA EM GP', 'LEICA KF80', 'LEICA PLUNGER', 'REICHERT-JUNG PLUNGER', 'OTHER']
                        # if vitrification_instrument in allowed_vitrification_instruments:
                        vitr.set_instrument(vitrification_instrument)
                        #else:
                        # vitr.set_instrument('OTHER')
                        # element 5 - <xs:complexType name="vitrifType">
                        # XSD: <xs:element name="method" type="xs:string" minOccurs="0"/>
                        self.check_set(vitr_in.get_method, vitr.set_method)
                        # element 6 - <xs:complexType name="vitrifType">
                        # XSD: <xs:element name="timeResolvedState" type="xs:string" minOccurs="0"/>
                        self.check_set(vitr_in.get_timed_resolved_state, vitr.set_timeResolvedState)
                        # element 7 - <xs:complexType name="vitrifType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        self.check_set(vitr_in.get_details, vitr.set_details)
                    else:
                        # negative staining
                        vitr.set_cryogenName('NONE')

                    exp.add_vitrification(vitr)

            # element 2 - <xs:complexType name="expType">
            # XSD: <xs:element name="imaging" type="imgType" maxOccurs="unbounded"/>
            mic_list_in = sd_in.get_microscopy_list().get_microscopy()
            for mic_in in mic_list_in:
                # XSD: <xs:complexType name="imgType"> has 26 elements
                img = emdb_19.imgType()
                # element 1 - <xs:complexType name="imgType">
                # XSD: <xs:element name="astigmatism" type="xs:string" minOccurs="0"/>
                ali_in = mic_in.get_alignment_procedure()
                if ali_in is not None:
                    leg_ali_in = ali_in.get_legacy()
                    if leg_ali_in is not None:
                        self.check_set(leg_ali_in.get_astigmatism, img.set_astigmatism)
                # element 2 - <xs:complexType name="imgType">
                # XSD: <xs:element name="electronSource" type="eSourceType"/>
                self.check_set(mic_in.get_electron_source, img.set_electronSource)
                # element 3 - <xs:complexType name="imgType">
                # XSD: <xs:element name="electronDose" type="eDoseType" minOccurs="0"/>
                im_rec_list_in = mic_in.get_image_recording_list()
                if im_rec_list_in is not None and im_rec_list_in != []:
                    im_recs = im_rec_list_in.get_image_recording()
                    for im_rec_in in im_recs:
                        self.set_value_and_units(im_rec_in.get_average_electron_dose_per_image, img.set_electronDose, constructor=emdb_19.eDoseType, units=const.U_EL_A2)
                # element 4 - <xs:complexType name="imgType">
                # XSD: <xs:element name="energyFilter" type="xs:string" minOccurs="0"/>
                sp_op_in = mic_in.get_specialist_optics()
                if sp_op_in is not None:
                    egf = sp_op_in.get_energy_filter()
                    if egf is not None:
                        img.set_energyFilter(egf.get_name())
                # element 5 - <xs:complexType name="imgType">
                # XSD: <xs:element name="imagingMode" type="imgModeType"/>
                self.check_set(mic_in.get_imaging_mode, img.set_imagingMode)
                # element 6 - <xs:complexType name="imgType">
                # XSD: <xs:element name="nominalDefocusMin" type="defocusType" minOccurs="0"/>
                nom_defocus_min = mic_in.get_nominal_defocus_min()
                mic_details = mic_in.get_details()
                if nom_defocus_min is not None:
                    nom_defocus_min_val = float(nom_defocus_min.valueOf_) * 1000
                    if self.roundtrip:
                        if mic_details is not None:
                            if mic_details.find('{nominal defocus min is int}') != -1:
                                nom_defocus_min_val = int(nom_defocus_min_val)
                    img.set_nominalDefocusMin(emdb_19.defocusType(valueOf_=nom_defocus_min_val, units='nm'))
                # element 7 - <xs:complexType name="imgType">
                # XSD: <xs:element name="nominalDefocusMax" type="defocusType" minOccurs="0"/>
                nom_defocus_max = mic_in.get_nominal_defocus_max()
                if nom_defocus_max is not None:
                    nom_defocus_max_val = float(nom_defocus_max.valueOf_) * 1000
                    if self.roundtrip:
                        if mic_details is not None:
                            if mic_details.find('{nominal defocus max is int}') != -1:
                                nom_defocus_max_val = int(nom_defocus_max_val)
                    img.set_nominalDefocusMax(emdb_19.defocusType(valueOf_=nom_defocus_max_val, units='nm'))
                # element 8 - <xs:complexType name="imgType">
                # XSD: <xs:element name="illuminationMode" type="illumType"/>
                self.check_set(mic_in.get_illumination_mode, img.set_illuminationMode)
                # element 9 - <xs:complexType name="imgType">
                # XSD: <xs:element name="specimenHolder" type="xs:string" minOccurs="0"/>
                self.check_set(mic_in.get_specimen_holder, img.set_specimenHolder)
                # element 10 - <xs:complexType name="imgType">
                # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                if mic_details is not None:
                    if mic_details.find('{nominal defocus min is int}') != -1:
                        mic_details = mic_details.replace('{nominal defocus min is int}', '')
                    if mic_details.find('{nominal defocus max is int}') != -1:
                        mic_details = mic_details.replace('{nominal defocus max is int}', '')
                    if mic_details != '':
                        img.set_details(mic_details)
                # element 11 - <xs:complexType name="imgType">
                # XSD: <xs:element name="detector" type="detectorType" minOccurs="0"/>
                if im_rec_list_in is not None and im_rec_list_in != []:
                    im_recs = im_rec_list_in.get_image_recording()
                    for im_rec_in in im_recs:
                        fod = im_rec_in.get_film_or_detector_model()
                        if fod is not None:
                            det_model_in = fod.get_valueOf_()
                            allowed_detectors = ['AGFA SCIENTA FILM', 'KODAK 4489 FILM', 'KODAK SO-163 FILM', 'GENERIC FILM',
                                                 'GENERIC IMAGE PLATES', 'DIRECT ELECTRON DE-10 (5k x 4k)', 'DIRECT ELECTRON DE-12 (4k x 3k)',
                                                 'DIRECT ELECTRON DE-16 (4k x 4k)', 'DIRECT ELECTRON DE-20 (5k x 3k)',
                                                 'DIRECT ELECTRON DE-64 (8k x 8k)', 'FALCON II', 'FEI CETA (4k x 4k)', 'FEI EAGLE (2k x 2k)',
                                                 'FEI EAGLE (4k x 4k)', 'FEI FALCON I (4k x 4k)', 'FEI FALCON II (4k x 4k)',
                                                 'FEI FALCON III (4k x 4k)', 'GATAN MULTISCAN', 'GATAN ORIUS SC200 (2k x 2k)',
                                                 'GATAN ORIUS SC600 (2.7k x 2.7k)', 'GATAN ORIUS SC1000 (4k x 2.7k)', 'GATAN ULTRASCAN 1000 (2k x 2k)',
                                                 'GATAN ULTRASCAN 4000 (4k x 4k)', 'GATAN ULTRASCAN 10000 (10k x 10k)', 'GATAN K2 (4k x 4k)',
                                                 'GATAN K2 BASE (4k x 4k)', 'GATAN K2 SUMMIT (4k x 4k)', 'GATAN K2 IS (4k x 4k)',
                                                 'GATAN K2 QUANTUM (4k x 4k)', 'GENERIC GATAN (2k x 2k)', 'GENERIC GATAN (4k x 4k)',
                                                 'GENERIC GATAN', 'PROSCAN TEM-PIV (2k x 2k)', 'SIA 15C (3k x 3k)', 'TVIPS TEMCAM-F816 (8k x 8k)',
                                                 'TVIPS TEMCAM-F415 (4k x 4k)', 'TVIPS TEMCAM-F416 (4k x 4k)', 'TVIPS TEMCAM-F216 (2k x 2k)',
                                                 'TVIPS TEMCAM-F224 (2k x 2k)', 'GENERIC TVIPS (2k x 2k)', 'GENERIC TVIPS (4k x 4k)',
                                                 'GENERIC TVIPS', 'GENERIC CCD (2k x 2k)', 'GENERIC CCD (4k x 4k)', 'GENERIC CCD', 'OTHER']
                            if det_model_in in allowed_detectors:
                                img.set_detector(det_model_in)
                            else:
                                img.set_detector('OTHER')
                # element 12 - <xs:complexType name="imgType">
                # XSD: <xs:element name="nominalCs" type="csType" minOccurs="0"/>
                ang_in_details = mic_in.get_nominal_cs()
                if ang_in_details is not None:
                    img.set_nominalCs(emdb_19.csType(valueOf_=ang_in_details.get_valueOf_(), units='mm'))
                # element 13 - <xs:complexType name="imgType">
                # XSD: <xs:element name="tiltAngleMin" type="tiltType" minOccurs="0"/>
                mic_type = mic_in.original_tagname_
                if mic_type == 'tomography_microscopy':
                    tilt_series_list_in = mic_in.get_tilt_series()
                    if len(tilt_series_list_in) > 0:
                        ts_in = tilt_series_list_in[0]
                        axis_in = ts_in.get_axis1()
                        self.set_value_and_units(axis_in.get_min_angle, img.set_tiltAngleMin, emdb_19.tiltType, const.U_DEGF)
                else:
                    tilt_min_in = mic_in.get_tilt_angle_min()
                    if tilt_min_in is not None:
                        img.set_tiltAngleMin(emdb_19.tiltType(valueOf_=tilt_min_in, units=const.U_DEGF))
                # element 14 - <xs:complexType name="imgType">
                # XSD: <xs:element name="calibratedMagnification" type="xs:float" minOccurs="0"/>
                self.check_set(mic_in.get_calibrated_magnification, img.set_calibratedMagnification)
                # element 15 - <xs:complexType name="imgType">
                # XSD: <xs:element name="tiltAngleMax" type="tiltType" minOccurs="0"/>
                if mic_type in ['subtomogram_averaging_microscopy', 'tomography_microscopy']:
                    tilt_series_list_in = mic_in.get_tilt_series()
                    if len(tilt_series_list_in) > 0:
                        ts_in = tilt_series_list_in[0]
                        axis_in = ts_in.get_axis1()
                        self.set_value_and_units(axis_in.get_max_angle, img.set_tiltAngleMax, emdb_19.tiltType, const.U_DEGF)
                        self.set_value_and_units(axis_in.get_min_angle, img.set_tiltAngleMin, emdb_19.tiltType, const.U_DEGF)
                else:
                    tilt_max_in = mic_in.get_tilt_angle_max()
                    if tilt_max_in is not None:
                        img.set_tiltAngleMax(emdb_19.tiltType(valueOf_=tilt_max_in, units=const.U_DEGF))
                    tilt_min_in = mic_in.get_tilt_angle_min()
                    if tilt_min_in is not None:
                        img.set_tiltAngleMin(emdb_19.tiltType(valueOf_=tilt_min_in, units=const.U_DEGF))
                # elements 16 - 18
                temp_in = mic_in.get_temperature()
                if temp_in is not None:
                    temp_av_in = temp_in.get_temperature_average()
                    temp_max_in = temp_in.get_temperature_max()
                    temp_min_in = temp_in.get_temperature_min()
                    if temp_av_in is not None:
                        # element 16 - <xs:complexType name="imgType">
                        # XSD: <xs:element name="temperature" type="tempType" minOccurs="0"/>
                        img.set_temperature(emdb_19.tempType(valueOf_=temp_av_in.get_valueOf_(), units=const.U_KELF))
                    if temp_max_in is not None:
                        # element 17 - <xs:complexType name="imgType">
                        # XSD: <xs:element name="temperatureMin" type="tempType" minOccurs="0"/>
                        img.set_temperatureMax(emdb_19.tempType(valueOf_=temp_max_in.get_valueOf_(), units=const.U_KELF))
                    if temp_min_in is not None:
                        # element 18 - <xs:complexType name="imgType">
                        # XSD: <xs:element name="temperatureMax" type="tempType" minOccurs="0"/>
                        img.set_temperatureMin(emdb_19.tempType(valueOf_=temp_min_in.get_valueOf_(), units=const.U_KELF))

                # element 19 - <xs:complexType name="imgType">
                # XSD: <xs:element name="microscope" type="microscopeType"/>
                img.set_microscope(mic_in.get_microscope())
                # element 20 - <xs:complexType name="imgType">
                # XSD: <xs:element name="date" type="xs:string" minOccurs="0"/>
                self.check_set(mic_in.get_date, img.set_date)
                # if ang_in_details is not None:
                # img.set_date(ang_in_details.strftime(const.EM_DATE_FORMAT).upper())
                # element 21 - <xs:complexType name="imgType">
                # XSD: <xs:element name="specimenHolderModel" type="specimenHolderType"/>
                specimen_holder_model_in = mic_in.get_specimen_holder_model() or 'OTHER'
                specimen_holder_model_out = const.SPECIMEN_HOLDER_30_to_19[specimen_holder_model_in] if specimen_holder_model_in in const.SPECIMEN_HOLDER_30_to_19 else specimen_holder_model_in
                img.set_specimenHolderModel(specimen_holder_model_out)
                # element 22 - <xs:complexType name="imgType">
                # XSD: <xs:element name="acceleratingVoltage" type="accVoltType" minOccurs="0"/>
                acc_vol_in = mic_in.get_acceleration_voltage()
                if acc_vol_in is not None:
                    img.set_acceleratingVoltage(emdb_19.accVoltType(valueOf_=acc_vol_in.get_valueOf_(), units='kV'))
                # element 23 - <xs:complexType name="imgType">
                # XSD: <xs:element name="nominalMagnification" type="xs:float" minOccurs="0"/>
                self.check_set(mic_in.get_nominal_magnification, img.set_nominalMagnification)
                # element 24 - <xs:complexType name="imgType">
                # XSD: <xs:element name="energyWindow" type="eWindowType" minOccurs="0"/>
                if sp_op_in is not None:
                    egf = sp_op_in.get_energy_filter()
                    if egf is not None:
                        eng_low_in = egf.get_lower_energy_threshold()
                        eng_high_in = egf.get_upper_energy_threshold()
                        e_text = None
                        if eng_low_in is not None:
                            if eng_low_in < 0:
                                e_text = 'none'
                            else:
                                # assume that both low and high are defined in this case
                                e_min_in = eng_low_in.get_valueOf_()
                                if e_min_in == '':
                                    e_min_in = 0
                                e_text = '%g-%g' % (float(e_min_in), float(eng_high_in.get_valueOf_()))
                            e_win = emdb_19.eWindowType(valueOf_=e_text, units='eV')
                            img.set_energyWindow(e_win)
                # element 25 - <xs:complexType name="imgType">
                # XSD: <xs:element name="detectorDistance" type="xs:string" minOccurs="0"/>
                if im_rec_list_in is not None and im_rec_list_in != []:
                    im_recs = im_rec_list_in.get_image_recording()
                    for im_rec_in in im_recs:
                        self.check_set(im_rec_in.get_detector_distance, img.set_detectorDistance)
                # element 26 - <xs:complexType name="imgType">
                # XSD: <xs:element name="electronBeamTiltParams" type="xs:string" minOccurs="0"/>
                if ali_in is not None:
                    leg_ali_in = ali_in.get_legacy()
                    if leg_ali_in is not None:
                        self.check_set(leg_ali_in.get_electron_beam_tilt_params, img.set_electronBeamTiltParams)

                exp.add_imaging(img)

            # element 3 - <xs:complexType name="expType">
            # XSD: <xs:element name="imageAcquisition" type="imgScanType" maxOccurs="unbounded"/>
            image_acquasitions = {}
            for mic_in in mic_list_in:
                im_rec_list_in = mic_in.get_image_recording_list()
                if im_rec_list_in is not None and im_rec_list_in != []:
                    im_recs = im_rec_list_in.get_image_recording()
                    for im_rec_in in im_recs:
                        image_acquasition = {}
                        # XSD: <xs:complexType name="imgScanType"> has 7 elements
                        im_ac = emdb_19.imgScanType()
                        # element 1 -<xs:complexType name="imgScanType">
                        # XSD: <xs:element name="numDigitalImages" type="xs:positiveInteger" minOccurs="0"/>
                        self.check_set(im_rec_in.get_number_real_images, im_ac.set_numDigitalImages)
                        image_acquasition['numDigitalImages'] = im_rec_in.get_number_real_images()
                        # elements 2 and 3
                        dig_in = im_rec_in.get_digitization_details()
                        if dig_in is not None:
                            # element 2 -<xs:complexType name="imgScanType">
                            # XSD: <xs:element name="scanner" type="scannerType" minOccurs="0"/>
                            self.check_set(dig_in.get_scanner, im_ac.set_scanner)
                            image_acquasition['scanner'] = dig_in.get_scanner()
                            # element 3 -<xs:complexType name="imgScanType">
                            # XSD: <xs:element name="samplingSize" type="samplSizeType" minOccurs="0"/>
                            self.set_value_and_units(dig_in.get_sampling_interval, im_ac.set_samplingSize, emdb_19.samplSizeType, const.U_MCRN)
                            sample_interval = dig_in.get_sampling_interval()
                            if sample_interval is not None:
                                image_acquasition['samplingSize'] = sample_interval.valueOf_
                        # element 4 -<xs:complexType name="imgScanType">
                        # XSD: <xs:element name="odRange" type="xs:float" minOccurs="0"/>
                        self.check_set(im_rec_in.get_od_range, im_ac.set_odRange)
                        image_acquasition['odRange'] = im_rec_in.get_od_range()
                        # element 5 -<xs:complexType name="imgScanType">
                        # XSD: <xs:element name="URLRawData" type="xs:string" minOccurs="0"/>
                        aux_list_in = x_ref_in.get_auxiliary_link_list()
                        if aux_list_in is not None:
                            ang_in_details = aux_list_in.get_auxiliary_link()
                            if len(ang_in_details) > 0:
                                im_ac.set_URLRawData(ang_in_details[0].get_link())
                                image_acquasition['URLRawData'] = ang_in_details[0].get_link()
                        # element 6 -<xs:complexType name="imgScanType">
                        # XSD: <xs:element name="quantBitNumber" type="xs:positiveInteger" minOccurs="0"/>
                        self.check_set(im_rec_in.get_bits_per_pixel, im_ac.set_quantBitNumber)
                        image_acquasition['quantBitNumber'] = im_rec_in.get_bits_per_pixel()
                        # element 7 -<xs:complexType name="imgScanType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        self.check_set(im_rec_in.get_details, im_ac.set_details)
                        image_acquasition['details'] = im_rec_in.get_details()
                        if image_acquasitions != {}:
                            if image_acquasition != {}:
                                # there are other acquisitions; check if any is the same with the current one
                                is_same = False
                                for im_acq in image_acquasitions:
                                    if cmp(image_acquasitions[im_acq], image_acquasition) == 0:
                                        is_same = True
                                        break
                                if is_same is False:
                                    # no same image acquisitions - add this one to the list
                                    image_acquasitions[len(image_acquasitions) + 1] = image_acquasition
                                    if self.roundtrip:
                                        if im_ac.hasContent_():
                                            exp.add_imageAcquisition(im_ac)
                                    else:
                                        exp.add_imageAcquisition(im_ac)
                        else:
                            # there are no other image acquisitions
                            if image_acquasition != {}:
                                # add this acquisition to the list
                                image_acquasitions[len(image_acquasitions) + 1] = image_acquasition
                            if self.roundtrip:
                                if im_ac.hasContent_():
                                    exp.add_imageAcquisition(im_ac)
                            else:
                                exp.add_imageAcquisition(im_ac)

            # element 4 - <xs:complexType name="expType">
            # XSD: <xs:element name="fitting" type="fittingType" minOccurs="0" maxOccurs="unbounded"/>
            if intrp_in is not None:
                fit_list_in = intrp_in.get_modelling_list()
                if fit_list_in is not None:
                    fits_in = fit_list_in.get_modelling()
                    for fit_in in fits_in:
                        # XSD: <xs:complexType name="fittingType"> has 7 elements
                        fit = emdb_19.fittingType()
                        # element 1 - <xs:complexType name="fittingType">
                        # XSD: <xs:element name="pdbEntryIdList" type="pdbidList2Type"/>
                        pdb_list = emdb_19.pdbidList2Type()
                        # XSD: <xs:complexType name="pdbidList2Type"> has 2 elements
                        mods_in = fit_in.get_initial_model()
                        if len(mods_in) > 0:
                            for mod_in in mods_in:
                                # element 1 - <xs:complexType name="pdbidList2Type">
                                # XSD: <xs:element name="pdbEntryId" type="pdbidType" maxOccurs="unbounded"/>
                                acc_code = mod_in.get_access_code()
                                pdb_list.add_pdbEntryId(acc_code)
                                # element 1 - <xs:complexType name="pdbidList2Type">
                                # XSD: <xs:element name="pdbChainId" type="xs:string" minOccurs="0" maxOccurs="unbounded"/>
                                chain_list = mod_in.get_chain()
                                for chain in chain_list:
                                    ch_ids = chain.get_chain_id()
                                    for ch_id in ch_ids:
                                        fit_details = fit_in.get_details()
                                        if fit_details is not None:
                                            if fit_details.find("PDBEntryID_givenInChain. ") == -1:
                                                pdb_list.add_pdbChainId(ch_id)
                                            else:
                                                pdb_list.add_pdbChainId('%s_%s' % (acc_code, ch_id))
                                        else:
                                            pdb_list.add_pdbChainId(ch_id)

                        if self.roundtrip:
                            if pdb_list.hasContent_():
                                fit.set_pdbEntryIdList(pdb_list)
                        else:
                            fit.set_pdbEntryIdList(pdb_list)
                        # element 2 - <xs:complexType name="fittingType">
                        # XSD: <xs:element name="software" type="xs:string" minOccurs="0"/>
                        soft_list_in = fit_in.get_software_list()
                        if soft_list_in is not None:
                            soft_str = make_software_from_list(soft_list_in.get_software())
                            if soft_str is not None:
                                fit.set_software(soft_str)
                        # element 3 - <xs:complexType name="fittingType">
                        # XSD: <xs:element name="refProtocol" type="refProtocolType" minOccurs="0"/>
                        ref_prot = fit_in.get_refinement_protocol()
                        if ref_prot is not None:
                            ref_prot_low = ref_prot.lower()
                            allowed_prots = ['rigid body', 'flexible']
                            known_issues = {'rigid body fit': 'rigid body'}
                            if ref_prot_low in allowed_prots:
                                fit.set_refProtocol(ref_prot_low)
                            elif ref_prot_low in known_issues:
                                fit.set_refProtocol(known_issues.get(ref_prot_low))
                        # element 4 - <xs:complexType name="fittingType">
                        # XSD: <xs:element name="targetCriteria" type="xs:string" minOccurs="0"/>
                        self.check_set(fit_in.get_target_criteria, fit.set_targetCriteria)
                        # element 5 - <xs:complexType name="fittingType">
                        # XSD: <xs:element name="overallBValue" type="xs:float" minOccurs="0"/>
                        self.check_set(fit_in.get_overall_bvalue, fit.set_overallBValue)
                        # element 6 - <xs:complexType name="fittingType">
                        # XSD: <xs:element name="refSpace" type="refSpaceType" minOccurs="0"/>
                        self.check_set(fit_in.get_refinement_space, fit.set_refSpace)
                        # element 7 - <xs:complexType name="fittingType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        fit_details = fit_in.get_details()
                        if fit_details is not None:
                            if fit_details.find('PDBEntryID_givenInChain. ') != -1:
                                fit_details = fit_details.replace('PDBEntryID_givenInChain. ', '')
                                if fit_details != '':
                                    fit.set_details(fit_details)
                            else:
                                self.check_set(fit_in.get_details, fit.set_details)

                        if fit.hasContent_():
                            exp.add_fitting(fit)

            # element 5 - <xs:complexType name="expType">
            # XSD: <xs:element name="specimenPreparation" type="smplPrepType"/>
            # Only the first element sets the specimen preparation element in 1.9 (as there is only one element and not an array as in 3.0)
            # in some cases for 2D arrays, a second element has been defined in the transfer from 1.9 -> 3.0.
            # ID 2 will in these cases contain crystal grow details
            sd_in_specimen_preparation_list = sd_in.get_specimen_preparation_list()
            if sd_in_specimen_preparation_list is not None:
                spec_prep_list_in = sd_in_specimen_preparation_list.get_specimen_preparation()
                for sp_in in spec_prep_list_in:
                    sp_prep_type = sp_in.original_tagname_
                    sp_in_id = sp_in.get_preparation_id()
                    if sp_in_id == 1:
                        # XSD: <xs:complexType name="smplPrepType"> has 9 elements
                        smpl_prep = emdb_19.smplPrepType()
                        # Forward reference for x-tal image processing
                        spec_prep_1 = smpl_prep
                        # element 1 - <xs:complexType name="smplPrepType">
                        # XSD: <xs:element name="specimenState" type="specimenType" minOccurs="0"/>
                        agg_state_in = sd_in.get_aggregation_state()
                        if agg_state_in is not None:
                            agg_state_in = agg_state_in.lower()
                            agg_state_out = const.SPECIMEN_STATE_30_to_19[agg_state_in] if agg_state_in in const.SPECIMEN_STATE_30_to_19 else const.SPECIMEN_STATE_30_to_19.itervalues().next()
                            smpl_prep.set_specimenState(agg_state_out)
                        # element 2 - <xs:complexType name="smplPrepType">
                        # XSD: <xs:element name="specimenConc" type="samplConcType" minOccurs="0"/>
                        conc_in = sp_in.get_concentration()
                        if conc_in is not None:
                            conc = emdb_19.samplConcType()
                            conc.set_units('mg/ml')
                            conc.set_valueOf_(conc_in.get_valueOf_())
                            smpl_prep.set_specimenConc(conc)
                        # element 3 - <xs:complexType name="smplPrepType">
                        # XSD: <xs:element name="buffer" type="bufferType" minOccurs="0"/>
                        buf_in = sp_in.get_buffer()
                        if buf_in is not None:
                            buf = emdb_19.bufferType()
                            self.check_set(buf_in.get_ph, buf.set_ph)
                            self.check_set(buf_in.get_details, buf.set_details)
                            smpl_prep.set_buffer(buf)
                        # element 4 - <xs:complexType name="smplPrepType">
                        # XSD: <xs:element name="staining" type="xs:string" minOccurs="0"/>
                        stain_in = sp_in.get_staining()
                        if stain_in is not None:
                            self.check_set(stain_in.get_details, smpl_prep.set_staining)
                        # element 5 - <xs:complexType name="smplPrepType">
                        # XSD: <xs:element name="specimenSupportDetails" type="xs:string" minOccurs="0"/>
                        grid_in = sp_in.get_grid()
                        if grid_in is not None:
                            self.check_set(grid_in.get_details, smpl_prep.set_specimenSupportDetails)
                        # element 6 -8
                        # element 6 - <xs:complexType name="smplPrepType">
                        # XSD: <xs:element name="twoDCrystalParameters" type="twoDxtalParamType" minOccurs="0" maxOccurs="1"/>
                        # element 7 - <xs:complexType name="smplPrepType">
                        # XSD: <xs:element name="threeDCrystalParameters" type="threeDxtalParamType" minOccurs="0" maxOccurs="1"/>
                        # element 8 - <xs:complexType name="smplPrepType">
                        # XSD: <xs:element name="helicalParameters" type="helixParamType" minOccurs="0" maxOccurs="1"/>
                        imp_list_in = sd_in.get_image_processing()
                        if imp_list_in is not None and len(imp_list_in) > 0:
                            em_method = sd_in.get_method()
                            for imp_in in imp_list_in:
                                final_reconstruct_in = imp_in.get_final_reconstruction()
                                if em_method == const.EMM_EC:
                                    set_crystal_parameters(imp_in, spec_prep_1)
                                    set_helical_parameters(final_reconstruct_in, spec_prep_1)
                                elif em_method == const.EMM_HEL:
                                    set_crystal_parameters(imp_in, spec_prep_1)
                                    set_helical_parameters(final_reconstruct_in, spec_prep_1)
                                elif em_method == const.EMM_STOM:
                                    set_crystal_parameters(imp_in, spec_prep_1)
                                    set_helical_parameters(final_reconstruct_in, spec_prep_1)
                                elif em_method == const.EMM_TOM:
                                    set_crystal_parameters(imp_in, spec_prep_1)
                                    set_helical_parameters(final_reconstruct_in, spec_prep_1)
                                elif em_method == const.EMM_SP:
                                    set_helical_parameters(final_reconstruct_in, spec_prep_1)
                        # element 9 - <xs:complexType name="smplPrepType">
                        # XSD: <xs:element name="crystalGrowDetails" type="xs:string" minOccurs="0"/>
                        if sp_prep_type == 'crystallography_preparation':
                            cryst_form = sp_in.get_crystal_formation()
                            if cryst_form is not None:
                                self.check_set(cryst_form.get_details, smpl_prep.set_crystalGrowDetails)
                        else:
                            # if other preparations have crystalGrowDetails they are recorded as:
                            # crystalGrowDetails: text :crystalGrowDetails
                            det = sp_in.get_details()
                            if det is not None and det.find(' :crystalGrowDetails') != -1:
                                det2 = det.split(' :crystalGrowDetails')
                                det3 = det2
                                if len(det2) > 1:
                                    det3 = det2[0]
                                if det3.find('crystalGrowDetails: ') != -1:
                                    grow_det = det3.replace('crystalGrowDetails: ', '')
                                    smpl_prep.set_crystalGrowDetails(grow_det)
                        if spec_prep_1.hasContent_():
                            exp.set_specimenPreparation(spec_prep_1)
                    elif sp_in_id != 1 and vitr_in is None:
                        if sp_prep_type == 'crystallography_preparation':
                            cryst_form = sp_in.get_crystal_formation()
                            if cryst_form is not None:
                                # element 9 - <xs:complexType name="smplPrepType">
                                # XSD: <xs:element name="crystalGrowDetails" type="xs:string" minOccurs="0"/>
                                self.check_set(cryst_form.get_details, spec_prep_1.set_crystalGrowDetails)
                        else:
                            det = sp_in.get_details()
                            if det is not None and det.find(' :crystalGrowDetails') != -1:
                                det2 = det.split(' :crystalGrowDetails')
                                det3 = det2
                                if len(det2) > 1:
                                    det3 = det2[0]
                                if det3.find('crystalGrowDetails: ') != -1:
                                    grow_det = det3.replace('crystalGrowDetails: ', '')
                                    spec_prep_1.set_crystalGrowDetails(grow_det)
            xml_out.set_experiment(exp)

        # element 7 - <xs:complexType name="entryType">
        # XSD: <xs:element name="processing" type="processType" maxOccurs="1"/>
        # XSD: <xs:complexType name="processType"> has 2 elements + 1 that is a choice of 5 elements
        proc = emdb_19.processType()
        # processing is created from the structure determination in schema v3.0 strict
        imp_list_in = None
        if sd_in is not None:
            imp_list_in = sd_in.get_image_processing()
        if imp_list_in is not None and len(imp_list_in) > 0:
            em_method = sd_in.get_method()
            for imp_in in imp_list_in:
                # element 1 - <xs:complexType name="processType">
                # XSD: <xs:element name="method" type="methodType"/>
                if imp_in.get_image_processing_id() == 1:
                    if em_method in [const.EMM_SP, const.EMM_STOM, const.EMM_TOM, const.EMM_HEL]:
                        proc.set_method(em_method)
                    elif em_method == const.EMM_EC:
                        proc.set_method('twoDCrystal')
                # element 2 - <xs:complexType name="processType">
                # XSD: <xs:element name="reconstruction" type="reconsType" maxOccurs="unbounded"/>
                # XSD: <xs:complexType name="reconsType"> has 7 elements
                rec = emdb_19.reconsType()
                # reconstruction is taken from final reconstruction in image processing
                final_reconstruct_in = imp_in.get_final_reconstruction()
                if final_reconstruct_in is not None:
                    # element 1 - <xs:complexType name="reconsType">
                    # XSD: <xs:element name="algorithm" type="xs:string" minOccurs="0"/>
                    self.check_set(final_reconstruct_in.get_algorithm, rec.set_algorithm)
                    # element 2 - <xs:complexType name="reconsType">
                    # XSD: <xs:element name="software" type="xs:string" minOccurs="0"/>
                    soft_list_in = final_reconstruct_in.get_software_list()
                    if soft_list_in is not None:
                        soft_str = make_software_from_list(soft_list_in.get_software())
                        if soft_str is not None:
                            rec.set_software(soft_str)
                    # element 3 - <xs:complexType name="reconsType">
                    # XSD: <xs:element name="ctfCorrection" type="xs:string" minOccurs="0"/>
                    # set in the copy_ctf_and_euler_angles() call below
                    # element 4 - <xs:complexType name="reconsType">
                    # XSD: <xs:element name="resolutionByAuthor" type="xs:string" minOccurs="0"/>
                    res_in = final_reconstruct_in.get_resolution()
                    if res_in is not None:
                        rec.set_resolutionByAuthor(res_in.get_valueOf_())
                    # element 5 - <xs:complexType name="reconsType">
                    # XSD: <xs:element name="resolutionMethod" type="xs:string" minOccurs="0"/>
                    self.check_set(final_reconstruct_in.get_resolution_method, rec.set_resolutionMethod)
                    # element 6 - <xs:complexType name="reconsType">
                    # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    self.check_set(final_reconstruct_in.get_details, rec.set_details)
                    # element 7 - <xs:complexType name="reconsType">
                    # XSD: <xs:element name="eulerAnglesDetails" type="xs:string" minOccurs="0"/>
                    # set in the copy_ctf_and_euler_angles() call below
                proc.add_reconstruction(rec)
                # element 3 - <xs:complexType name="processType">
                proc_spec = None
                im_proc_in_id = imp_in.get_image_processing_id()
                if im_proc_in_id == 1:
                    # choice 1
                    # XSD: <xs:element name="twoDCrystal" type="xtal2DType" maxOccurs="1"/>
                    if em_method == const.EMM_EC:
                        # XSD <xs:complexType name="xtal2DType"> has 1 element
                        proc_spec = emdb_19.xtal2DType()
                        # set_crystal_parameters(imp_in, spec_prep_1)
                        # set_helical_symmetry(final_reconstruct_in, spec_prep_1)
                        # element 1 - <xs:complexType name="xtal2DType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        self.check_set(imp_in.get_details, proc_spec.set_details)
                        proc.set_twoDCrystal(proc_spec)
                    # choice 2
                    # XSD: <xs:element name="helical" type="helixType" maxOccurs="1"/>
                    elif em_method == const.EMM_HEL:
                        if final_reconstruct_in is not None:
                            alg_in = final_reconstruct_in.get_algorithm()
                            if alg_in is not None:
                                mtch = re.match(const.HEL_SP_PAT, alg_in)
                                if mtch is not None:
                                    match_groups = mtch.groups()
                                    hx_method_in = match_groups[1]
                                    alg_str = match_groups[0] + match_groups[2]
                                else:
                                    hx_method_in = const.EMM_HEL
                                    alg_str = alg_in
                                rec.set_algorithm(alg_str)
                            else:
                                hx_method_in = const.EMM_HEL
                            if hx_method_in == const.EMM_HEL:
                                # XSD: <xs:complexType name="helixType"> has 1 element
                                hel = emdb_19.helixType()
                                # element 1 - <xs:complexType name="helixType">
                                # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                                self.check_set(imp_in.get_details, hel.set_details)
                                # set_crystal_parameters(imp_in, spec_prep_1)
                                proc.set_helical(hel)
                            # set_helical_symmetry(final_reconstruct_in, spec_prep_1)
                    # choice 3
                    # XSD: <xs:element name="subtomogramAveraging" type="subTomType" maxOccurs="1"/>
                    elif em_method == const.EMM_STOM:
                        # XSD: <xs:complexType name="subTomType"> has 4 elements
                        proc_spec = emdb_19.subTomType()
                        # element 1 - <xs:complexType name="subTomType">
                        # XSD: <xs:element name="appliedSymmetry" type="xs:string" minOccurs="0"/>
                        if final_reconstruct_in is not None:
                            symm_in = final_reconstruct_in.get_applied_symmetry()
                            if symm_in is not None:
                                self.check_set(symm_in.get_point_group, proc_spec.set_appliedSymmetry)
                            # element 2 - <xs:complexType name="subTomType">
                            # XSD: <xs:element name="numSubtomograms" type="xs:positiveInteger" minOccurs="0"/>
                            self.check_set(final_reconstruct_in.get_number_subtomograms_used, proc_spec.set_numSubtomograms)
                        # element 3 - <xs:complexType name="subTomType">
                        # XSD: <xs:element name="numClassAverages" type="xs:positiveInteger" minOccurs="0"/>
                        sav_cls_in = imp_in.get_final_three_d_classification()
                        if sav_cls_in is not None:
                            self.check_set(sav_cls_in.get_number_classes, proc_spec.set_numClassAverages)
                        # element 4 - <xs:complexType name="subTomType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        self.check_set(imp_in.get_details, proc_spec.set_details)

                        # set_crystal_parameters(imp_in, spec_prep_1)
                        # set_helical_symmetry(final_reconstruct_in, spec_prep_1)
                        proc.set_subtomogramAveraging(proc_spec)
                    # choice 4
                    # XSD: <xs:element name="tomography" type="tomogrType" maxOccurs="1"/>
                    elif em_method == const.EMM_TOM:
                        # XSD: <xs:complexType name="tomogrType"> has 4 elements
                        proc_spec = emdb_19.tomogrType()
                        # element 1 - <xs:complexType name="tomogrType">
                        # XSD: <xs:element name="appliedSymmetry" type="pointGroupSymmetryType" minOccurs="0"/>
                        symm_in = final_reconstruct_in.get_applied_symmetry()
                        if symm_in is not None:
                            self.check_set(symm_in.get_point_group, proc_spec.set_appliedSymmetry)
                        # element 2 - <xs:complexType name="tomogrType">
                        # XSD: <xs:element name="tiltAngleIncrement" type="xs:string" minOccurs="0"/>
                        mic_list_in = sd_in.get_microscopy_list().get_microscopy()
                        tom_mic = mic_list_in[0]
                        tilt_series_in = tom_mic.get_tilt_series()
                        if len(tilt_series_in) > 0:
                            tilt_series = tilt_series_in[0]
                            axis_1_tilt_series = tilt_series.get_axis1()
                            tilt_inc = axis_1_tilt_series.get_angle_increment()
                            if tilt_inc is not None:
                                proc_spec.set_tiltAngleIncrement(tilt_inc.get_valueOf_())
                        # element 3 - <xs:complexType name="tomogrType">
                        # XSD: <xs:element name="numSections" type="xs:positiveInteger" minOccurs="0"/>
                        self.check_set(final_reconstruct_in.get_number_images_used, proc_spec.set_numSections)
                        # element 4 - <xs:complexType name="tomogrType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        self.check_set(imp_in.get_details, proc_spec.set_details)
                        # set_crystal_parameters(imp_in, spec_prep_1)
                        # set_helical_symmetry(final_reconstruct_in, spec_prep_1)
                        proc.set_tomography(proc_spec)
                    # choice 5
                    # XSD: <xs:element name="singleParticle" type="singPartType" maxOccurs="1"/>
                    elif em_method == const.EMM_SP:
                        # XSD: <xs:complexType name="singPartType"> has 4 elements
                        proc_spec = emdb_19.singPartType()
                        # set_el_single_particle(imp_in, final_reconstruct_in, proc_spec)
                        # element 1 - <xs:complexType name="singPartType">
                        # XSD: <xs:element name="appliedSymmetry" type="pointGroupSymmetryType" minOccurs="0"/>
                        symm_in = final_reconstruct_in.get_applied_symmetry()
                        if symm_in is not None:
                            self.check_set(symm_in.get_point_group, proc_spec.set_appliedSymmetry)
                        # element 2 - <xs:complexType name="singPartType">
                        # XSD: <xs:element name="numProjections" type="xs:positiveInteger" minOccurs="0"/
                        self.check_set(final_reconstruct_in.get_number_images_used, proc_spec.set_numProjections)
                        # element 3 - <xs:complexType name="singPartType">
                        # XSD: <xs:element name="numClassAverages" type="xs:positiveInteger" minOccurs="0"/>
                        sp_cls_in = imp_in.get_final_two_d_classification()
                        if sp_cls_in is not None:
                            self.check_set(sp_cls_in.get_number_classes, proc_spec.set_numClassAverages)
                        # element 4 - <xs:complexType name="singPartType">
                        # XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        self.check_set(imp_in.get_details, proc_spec.set_details)
#                         set_helical_symmetry(final_reconstruct_in, spec_prep_1)
                        proc.set_singleParticle(proc_spec)

                # Euler angles and ctf have to be set for all reconstruction objects
                copy_ctf_and_euler_angles(imp_in, rec, proc_spec)

#         if spec_prep_1.hasContent_():
#             exp.set_specimenPreparation(spec_prep_1)
        xml_out.set_processing(proc)
        # ---------------------------------
        # Write XML to file
        xml_v19_file = open(output_file, 'w') if output_file else sys.stdout
        xml_v19_file.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        # XSD: <xs:element name="emdEntry" type="entryType"/>
        xml_out.export(xml_v19_file, 0, name_='emdEntry')

        if xml_v19_file is not sys.stdout:
            xml_v19_file.close()

        # validate the output v1.9 file
        if self.validate_xml_out:
            print '**validation v1.9' + '*'*100
            self.validate(output_file, EMDBSettings.schema19)
            print '*'*117

    def validate_file(self, the_parser, xml_filename):
        """
        Method to validate any schema against any file
        """
        try:
            xml_file = open(xml_filename, 'r')
            try:
                etree.fromstring(xml_file.read(), the_parser)
            except etree.XMLSyntaxError:# as err:
                #self.logger.critical(self.Constants.VALIDATION_ERROR+" %s XMLSyntaxError for %s", err, xml_filename)
                return False
            except etree.XMLSchemaError:# as err:
                #self.logger.critical(self.Constants.VALIDATION_ERROR+" %s XMLSchemaError for %s", err, xml_filename)
                return False
            return True
        except IOError as exp:
            # self.logger.critical(self.Constants.VALIDATION_ERROR+"Error %s occured. Arguments %s.", exp.message, exp.args)
            return False
        finally:
            xml_file.close()

    def show_validation_errors(self, in_xml, in_schema_filename):
        """
        Called if the validation of the in_xml file against
        the schema in_schema_filename fails. Shows the list of validation errors.
        """
        try:
            xml_doc = etree.parse(in_xml)
            xsd = etree.parse(in_schema_filename)
            xml_schema = etree.XMLSchema(xsd)
            xml_schema.assertValid(xml_doc)
            print 'File %s validates' % in_xml
        except etree.XMLSyntaxError as exp:
            print 'PARSING ERROR %s' % exp
        except etree.DocumentInvalid as exp:
            i = 1
            for err in exp.error_log:
                print 'VALIDATION ERROR %d: %s' % (i, err)
                print ''
                i = i + 1

    def validate(self, in_xml, in_schema_filename):
        """
        Validate in_xml against in_schema
        """
        try:
            in_schema = open(in_schema_filename, 'r')
        except IOError as exp:
            print 'Validation error %s occurred. Arguments %s.' % (exp.message, exp.args)
            return False
        else:
            schema_doc = in_schema.read()
            schema_root = etree.XML(schema_doc)
            the_schema = etree.XMLSchema(schema_root)

            xml_parser = etree.XMLParser(schema=the_schema)
            validates = self.validate_file(xml_parser, in_xml)
            if not validates:
                # self.validation_logger_header(self.logger.critical, in_schema_filename)
                self.show_validation_errors(in_xml, in_schema_filename)
            return validates
        finally:
            in_schema.close()


def main():
    """
    Convert EMDB XML files from one schema version to another

    """

    # Handle command line options
    usage = """
            emdb_xml_translate.py [options] input_file
            Convert EMDB XML files from one schema version to another

            Examples:
            python emdb_xml_translate.py input_file

            Typical run:
            python emdb_xml_translate.py -f out.xml -i 1.9 -o 3.0 in.xml
            in.xml is assumed to be a EMDB 1.9 XML file and converted to
            an XML file following EMDB XML schema 3.0 and written out to out.xml
            """
    version = "0.19"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-i", "--in-schema", action="store", type="string", metavar="SCHEMA", dest="inputSchema", default="1.9", help="Schema version of output file - 1.9 or 3.0 [default: %default]")
    parser.add_option("-o", "--out-schema", action="store", type="string", metavar="SCHEMA", dest="outputSchema", default="3.0", help="Schema version of output file - 1.9 or 3.0 [default: %default]")
    parser.add_option("-f", "--out-file", action="store", type="string", metavar="FILE", dest="outputFile", help="Write output to FILE")
    parser.add_option("-w", "--warning-level", action="store", type="int", dest="warning_level", default=1, help="Level of warning output. 0 is none, 3 is max, default = 1")
    parser.add_option("-v", "--validate_ouput_xml", action="store_true", dest="validate", default=False, help="Validation flag. If true the output xml file will be validated.")
    parser.add_option("-r", "--use_relaxed_schema_v30", action="store_true", dest="relaxed", default=False, help="Schema v3.0 version flag. If true the schema is relaxed.")
    parser.add_option("-p", "--roundrip", action="store_true", dest="roundtrip", default=False, help="Roudtrip flag.If true the roundtrip files are created.")
    (options, args) = parser.parse_args()

    # Check for sensible/supported options
    if len(args) < 1:
        sys.exit("No input file specified!")
    else:
        input_file = args[0]
    if (options.inputSchema != "1.9" and options.outputSchema != "3.0") and (options.inputSchema != "3.0" and options.outputSchema != "1.9"):
        sys.exit("Conversion from version %s to %s not supported!" % (options.inputSchema, options.outputSchema))

    # Call appropriate conversion routine
    translator = EMDBXMLTranslator()
    translator.set_warning_level(options.warning_level)
    translator.set_validate(options.validate)
    translator.set_v30_schema(options.relaxed)
    translator.set_roundtrip(options.roundtrip)
    if options.inputSchema == "1.9" and options.outputSchema == "3.0":
        translator.translate_1_9_to_3_0(input_file, options.outputFile)
    if options.inputSchema == "1.9" and options.outputSchema == "1.9":
        translator.translate_1_9_to_1_9(input_file, options.outputFile)
    if options.inputSchema == "3.0" and options.outputSchema == "1.9":
        translator.translate_3_0_to_1_9(input_file, options.outputFile)


if __name__ == "__main__":
    main()
