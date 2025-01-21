#!/usr/bin/python
"""
relax_schema.py

Relaxes EMDB v3.x schemas.

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
import os
import xml.etree.ElementTree as ET


__author__ = 'Sanja Abbott'
__email__ = 'sanja@ebi.ac.uk'
__date__ = '2018-10-29'


def create_schema_filename(schema_filename_in):
    schema_filename, schema_extension = schema_filename_in.split(".")
    schema_filename_out = schema_filename + '_relaxed.' + schema_extension
    return schema_filename_out


def schema_file_from(schema_filename):
    if os.path.dirname(schema_filename) == '':
        schema_file = os.path.join(os.getcwd(), schema_filename)
    else:
        schema_file = schema_filename

    return schema_file


def child_elem_recursive(element, names):
    if len(names) == 1:
        return element, names
    else:
        for elem in element.iter():
            if elem.attrib.get('name') == names[0]:
                remain_names = names[1:]
                elem, remain_names = child_elem_recursive(elem, remain_names)
                return elem, remain_names


def relax_child(root, child_type, child_name, subs_names):
    for element in root.iter('{http://www.w3.org/2001/XMLSchema}'+child_type):
        if element.attrib.get('name') == child_name:
            last_elem, last_name = child_elem_recursive(element, subs_names)
            for elem_in_last_elem in last_elem.iter():
                if elem_in_last_elem.attrib.get('name') == last_name[0]:
                    print('relax %s in %s' % (elem_in_last_elem.attrib.get('name'), last_elem.attrib.get('name')))
                    if 'minOccurs' not in elem_in_last_elem.attrib.keys():
                        elem_in_last_elem.attrib['minOccurs'] = "0"
            break


def add_enum(root, child_type, child_name, new_enumeration, subchild_name=None):
    for element in root.iter('{http://www.w3.org/2001/XMLSchema}'+child_type):
        if element.attrib.get('name') == child_name:
            if subchild_name is not None:
                for elem in element.iter():
                    if elem.attrib.get('name') == subchild_name:
                        print('in <%s><%s> add %s' % (child_name, subchild_name, new_enumeration))
                        for el in elem.iter('{http://www.w3.org/2001/XMLSchema}restriction'):
                            new_enum = ET.SubElement(el, 'xs:enumeration')
                            new_enum.attrib['value'] = new_enumeration
            else:
                print('in <%s> add %s' % (child_name, new_enumeration))
                for el in element.iter('{http://www.w3.org/2001/XMLSchema}restriction'):
                    new_enum = ET.SubElement(el, 'xs:enumeration')
                    new_enum.attrib['value'] = new_enumeration


def comment_out_child(root, child_type, child_name, child_to_comment_out):
    for element in root.iter('{http://www.w3.org/2001/XMLSchema}'+child_type):
        if element.attrib.get('name') == child_name:
            for elem in element.iter():
                if elem.attrib.get('name') == child_to_comment_out:
                    print('in %s comment out %s' % (child_name, child_to_comment_out))
                    pattern_el = None
                    for el in elem.iter('{http://www.w3.org/2001/XMLSchema}pattern'):
                        pattern_el = el
                    restriction_el = None
                    for el in elem.iter('{http://www.w3.org/2001/XMLSchema}restriction'):
                        restriction_el = el
                    str_pattern_el = ET.tostring(pattern_el)
                    restriction_el.remove(pattern_el)
                    comment = ET.Comment(str_pattern_el)
                    restriction_el.append(comment)

def relax(schema_filename_in, schema_filename_out=None):
    """
    differences listed: https://www.ebi.ac.uk/seqdb/confluence/display/PDBE/Differences+between+relaxed+v3+schema+and+v3+schema
    :param schema_filename_in: Schema to relax
    :param schema_filename_out: Relaxed schema
    """
    schema_file_in = schema_file_from(schema_filename_in)
    if os.path.isfile(schema_file_in):
        print('schema in: %s' % schema_file_in)
        if schema_filename_out is None:
            schema_filename_out = create_schema_filename(schema_filename_in)
        schema_file_out = schema_file_from(schema_filename_out)
        print('schema out: %s' % schema_file_out)

        # Read the schema to releax
        tree = ET.parse(schema_file_in)
        root = tree.getroot()
        print('root is %s' % root)
        # 1. in <xs:complexType name="supersedes_type"> relax <xs:element name="date" type="xs:date"/>
        relax_child(root, 'complexType', 'supersedes_type', ('date',))
        # 2. in <xs:element name="journal_citation" substitutionGroup="citation_type"> relax <xs:element name="journal_abbreviation" type="xs:token"/>
        relax_child(root, 'element', 'journal_citation', ('journal_abbreviation',))
        # 3. in <xs:complexType name="auxiliary_link_type"> relax <xs:element name="type">
        relax_child(root, 'complexType', 'auxiliary_link_type', ('type',))
        # 4. in <xs:complexType name="base_supramolecule_type"> relax <xs:element name="parent" type="xs:nonNegativeInteger"/>
        relax_child(root, 'complexType', 'base_supramolecule_type', ('parent',))
        # 5. in <xs:complexType name="recombinant_source_type"> relax <xs:element name="recombinant_organism" type="organism_type"/>
        relax_child(root, 'complexType', 'recombinant_source_type', ('recombinant_organism',))
        # 6. in <xs:complexType name="virus_supramolecule_type"> <xs:extension base="base_supramolecule_type"> <xs:element name="virus_type"> add <xs:enumeration value="OTHER"/>
        add_enum(root, 'complexType', 'virus_supramolecule_type', 'OTHER', 'virus_type')
        # 7. in <xs:complexType name="base_macromolecule_type"> relax <xs:element name="name" type="sci_name_type"/>
        relax_child(root, 'complexType', 'base_macromolecule_type', ('name',))
        # 8. in <xs:complexType name="dna_macromolecule_type"><xs:extension base="base_macromolecule_type"> relax <xs:element name="sequence">
        relax_child(root, 'complexType', 'dna_macromolecule_type', ('sequence',))
        # 9. in <xs:complexType name="protein_or_peptide_macromolecule_type"><xs:extension base="base_macromolecule_type"> relax <xs:element name="enantiomer">
        relax_child(root, 'complexType', 'protein_or_peptide_macromolecule_type', ('enantiomer',))
        # 10. in <xs:complexType name="rna_macromolecule_type"><xs:extension base="base_macromolecule_type"> relax <xs:element name="sequence">
        relax_child(root, 'complexType', 'rna_macromolecule_type', ('sequence',))
        # 11. in <xs:complexType name="structure_determination_type"> relax <xs:element name="aggregation_state">
        relax_child(root, 'complexType', 'structure_determination_type', ('aggregation_state',))
        # 12. in <xs:complexType name="base_preparation_type"><xs:element name="staining" minOccurs="0"> relax <xs:element name="material" type="xs:token"/>
        relax_child(root, 'complexType', 'base_preparation_type', ('staining', 'material'))
        # 13. in <xs:complexType name="buffer_type"> relax <xs:element name="ph">
        relax_child(root, 'complexType', 'buffer_type', ('ph',))
        # 14. in <xs:complexType name="vitrification_type"> <xs:element name="cryogen_name"> add <xs:enumeration value="ETHANE-PROPANE MIXTURE"/>
        add_enum(root, 'complexType', 'vitrification_type', 'ETHANE-PROPANE MIXTURE', 'cryogen_name')
        # 14. in <xs:complexType name="vitrification_type"> <xs:element name="cryogen_name"> add <xs:enumeration value="NONE"/>
        add_enum(root, 'complexType', 'vitrification_type', 'NONE', 'cryogen_name')
        # 15. in <xs:complexType name="base_microscopy_type"> <xs:element name="microscope"> add <xs:enumeration value="OTHER"/>
        add_enum(root, 'complexType', 'base_microscopy_type', 'OTHER', 'microscope')
        # 16. in <xs:complexType name="base_microscopy_type"><xs:element name="alignment_procedure" minOccurs="0"><xs:element name="legacy"> relax <xs:element name="astigmatism" type="xs:string"/>
        relax_child(root, 'complexType', 'base_microscopy_type', ('alignment_procedure', 'legacy', 'astigmatism'))
        # 16. in <xs:complexType name="base_microscopy_type"><xs:element name="alignment_procedure" minOccurs="0"><xs:element name="legacy"> relax <xs:element name="electron_beam_tilt_params" type="xs:string"/>
        relax_child(root, 'complexType', 'base_microscopy_type', ('alignment_procedure', 'legacy', 'electron_beam_tilt_params'))
        # 17. in <xs:complexType name="base_microscopy_type"><xs:element name="image_recording_list"><xs:element name="image_recording"> relax <xs:element name="film_or_detector_model">
        relax_child(root, 'complexType', 'base_microscopy_type', ('film_or_detector_model',))
        # 18. in <xs:complexType name="crystallography_microscopy_type"><xs:extension base="base_microscopy_type"> relax <xs:element name="camera_length">
        relax_child(root, 'complexType', 'crystallography_microscopy_type', ('camera_length',))
        # 19. in <xs:complexType name="base_image_processing_type"> relax <xs:element name="image_recording_id" type="xs:positiveInteger"/>
        relax_child(root, 'complexType', 'base_image_processing_type', ('image_recording_id',))
        # 20. in <xs:complexType name="helical_parameters_type"> relax <xs:element name="delta_z">
        relax_child(root, 'complexType', 'helical_parameters_type', ('delta_z',))
        # 20. in <xs:complexType name="helical_parameters_type"> relax <xs:element name="delta_phi">
        relax_child(root, 'complexType', 'helical_parameters_type', ('delta_phi',))
        # 20. in <xs:complexType name="helical_parameters_type"> relax <xs:element name="axial_symmetry">
        relax_child(root, 'complexType', 'helical_parameters_type', ('axial_symmetry',))
        # 21. in <xs:complexType name="crystal_parameters_type"> relax <xs:element name="unit_cell" type="unit_cell_type"/>
        relax_child(root, 'complexType', 'crystal_parameters_type', ('unit_cell',))
        # 22. in <xs:complexType name="unit_cell_type"> relax <xs:element name="c" type="cell_type"/>
        relax_child(root, 'complexType', 'unit_cell_type', ('c',))
        # 23. in <xs:complexType name="angle_assignment_type"> relax <xs:element name="type">
        relax_child(root, 'complexType', 'angle_assignment_type', ('type',))
        # 24. in <xs:group name="subtomogram_averaging_proc_add_group"> relax <xs:element name="extraction">
        relax_child(root, 'group', 'subtomogram_averaging_proc_add_group', ('extraction',))
        # 25. in <xs:complexType name="map_type"> <xs:element name="file"> comment out/remove <xs:pattern value="emd_\d{4,}([A-Za-z0-9_]*)\.map(\.gz|)"/>
        comment_out_child(root, 'complexType', 'map_type', 'file')
        # 26. in <xs:simpleType name="reconstruction_algorithm_type"> add <xs:enumeration value="OTHER"/>
        add_enum(root, 'simpleType', 'reconstruction_algorithm_type', 'OTHER')
        # 27. in <xs:complexType name="base_microscopy_type"> relax <xs:element name="image_recording_list">
        relax_child(root, 'complexType', 'base_microscopy_type', ('image_recording_list',))
        # 28a. in <xs:group name="subtomogram_averaging_proc_add_group"><xs:element name="extraction"> relax <xs:element name="number_tomograms" type="xs:positiveInteger"/>
        relax_child(root, 'group', 'subtomogram_averaging_proc_add_group', ('extraction', 'number_tomograms'))
        # 28b. in <xs:group name="subtomogram_averaging_proc_add_group"><xs:element name="extraction"> relax <xs:element name="number_images_used" type="xs:positiveInteger"/>
        relax_child(root, 'group', 'subtomogram_averaging_proc_add_group', ('extraction', 'number_images_used'))
        # Write the relaxed schema
        if schema_file_out is not None:
            tree.write(schema_file_out, encoding="UTF-8", method="xml")
    else:
        print('Input schema file %s does not exist. It cannot be relaxed.' % schema_file_in)


def is_schema_name_correct(schema_name):
    if schema_name.endswith('.xsd'):
        return True
    else:
        print('Schema file name "%s" is not correct. It should have the ".xsd" extension' % schema_name)
        return False


def main():
    usage = """
    USAGE:
        python relax_schema.py 'myschema.xsd'
    where: 
        'myschema.xsd' is the input v3.x EMDB schema xsd file that will be relaxed
        the output file is 'myschema_relaxed.xsd'
        
        python relax_schema.py 'myschema.xsd' 'my_relaxed_schema_name.xsd'
    where:
        'myschema.xsd' is the input v3.x EMDB schema xsd file that will be relaxed
        the output file is given as 'my_relaxed_schema_name.xsd'
    """

    args = sys.argv[1:]
    if len(args) == 1:
        if is_schema_name_correct(args[0]):
            relax(args[0])
    elif len(args) == 2:
        if is_schema_name_correct(args[0]) and is_schema_name_correct(args[1]):
            relax(args[0], args[1])
    else:
        print('No input given. Please see the usage.')
        print(usage)


if __name__ == "__main__":
    main()

