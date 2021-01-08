#!/usr/bin/env python
# -*- mode: pymode; coding: latin1; -*-
"""
emdb_user_methods.py

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

import re

__author__ = 'Ardan Patwardhan, Sanja Abbott'
__email__ = 'ardan@ebi.ac.uk, sanja@ebi.ac.uk'
__date__ = '2017-06-14'
#
# You must include the following class definition at the top of
#   your method specification file.
#


class MethodSpec(object):
    """
    TO DO:
    """
    def __init__(self, name='', source='', class_names='', class_names_compiled=None):
        """
        MethodSpec -- A specification of a method.
        Member variables:
            name -- The method name
            source -- The source code for the method.  Must be
                indented to fit in a class definition.
            class_names -- A regular expression that must match the
                class names in which the method is to be inserted.
            class_names_compiled -- The compiled class names.
                generateDS.py will do this compile for you.
        """
        self.name = name
        self.source = source
        if class_names is None:
            self.class_names = ('.*', )
        else:
            self.class_names = class_names
        if class_names_compiled is None:
            self.class_names_compiled = re.compile(self.class_names)
        else:
            self.class_names_compiled = class_names_compiled

    def get_name(self):
        """
        TO DO
        """
        return self.name

    def set_name(self, name):
        """
        TO DO
        """
        self.name = name

    def get_source(self):
        """
        TO DO
        """
        return self.source

    def set_source(self, source):
        """
        TO DO
        """
        self.source = source

    def get_class_names(self):
        """
        TO DO
        """
        return self.class_names

    def set_class_names(self, class_names):
        """
        TO DO
        """
        self.class_names = class_names
        self.class_names_compiled = re.compile(class_names)

    def get_class_names_compiled(self):
        """
        TO DO
        """
        return self.class_names_compiled

    def set_class_names_compiled(self, class_names_compiled):
        """
        TO DO
        """
        self.class_names_compiled = class_names_compiled

    def match_name(self, class_name):
        """Match against the name of the class currently being generated.
        If this method returns True, the method will be inserted in
          the generated class.
        """
        return bool(self.class_names_compiled.search(class_name))

    def get_interpolated_source(self, values_dict):
        """Get the method source code, interpolating values from values_dict
        into it.  The source returned by this method is inserted into
        the generated class.
        """
        return self.source % values_dict

    def show(self):
        """
        TO DO
        """
        print 'specification:'
        print '    name: %s' % (self.name, )
        print self.source
        print '    class_names: %s' % (self.class_names, )
        print '    names pat  : %s' % (self.class_names_compiled.pattern, )


#
# Provide one or more method specification such as the following.
# Notes:
# - Each generated class contains a class variable _member_data_items.
#   This variable contains a list of instances of class _MemberSpec.
#   See the definition of class _MemberSpec near the top of the
#   generated superclass file and also section "User Methods" in
#   the documentation, as well as the examples below.

#
# Replace the following method specifications with your own.

#
# Sample method specification #1
#

method1 = MethodSpec(name='gds_format_float', source='''def gds_format_float(self,input_data, input_name="" ): return ("%%g" %% input_data)''',
                     class_names=r'^resolutionType$|^contourType$|^contourLevelType$|^chamber_humidityType$|^eWindowType$|^vitrifType$|^imgType$|^originType$|^limitType$',
                     )


#
# Provide a list of your method specifications.
#   This list of specifications must be named METHOD_SPECS.
#
METHOD_SPECS = (method1,)


def test():
    """
    TO DO:
    """
    for spec in METHOD_SPECS:
        spec.show()


def main():
    """
    TO DO:
    """
    test()


if __name__ == '__main__':
    main()
