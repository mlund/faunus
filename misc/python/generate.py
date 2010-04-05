#!/usr/bin/env python

#
# imports and preparation
#

import os
import sys
import fnmatch
from glob import glob

# get the location of the faunus source tree
faunus_base = os.environ['FAUNUS_BASE']

base = '../auxiliary/lib/python%s'
for version in ('2.5', '2.6'):
    if os.path.exists(base % version):
        path = (base + '/site-packages') % version
        sys.path.insert(0, path)
        print 'found path: %s' % path

# do we really want this?
#rm generated/*.cpp generated/*.hpp generated/*.h generated/*.txt

#from IPython.Shell import IPShellEmbed
#ipshell = IPShellEmbed([])

try:
    from pyplusplus.module_builder import module_builder_t
    from pygccxml.declarations.matchers import access_type_matcher_t
except:
    print "import from pyplusplus or pygccxml failed, aborting"
    sys.exit(1)

#
# helper code
#

class TemplateSpec(object):

    def __init__(self, name, namespace):
        self.name = name
        self.namespace = namespace
        self.arguments = []
        self.aliases = []

    def __str__(self):
        result = self.full_name()
        for class_name, alias in zip(self.class_names(), self.aliases):
            result += '\n    %s --- %s' % (class_name, alias)
        return result

    def add_instantiation(self, argument, alias):
        self.arguments.append(argument)
        self.aliases.append(alias)

    def full_name(self):
        return '%s::%s' % (self.namespace, self.name)

    def class_names(self, use_namespace=True):
        if use_namespace:
            full_name = self.full_name()
        else:
            full_name = self.name
        for argument in self.arguments:
            yield '%s< %s >' % (full_name, argument)

    def typedef_code(self, indent=0):
        typedef_lines = [indent*' ' + 'typedef %s %s;' % (class_name, alias)
                         for class_name, alias in zip(self.class_names(), self.aliases)
                        ]
        return '\n'.join(typedef_lines)

    def sizeof_code(self, indent=0):
        sizeof_lines = [indent*' ' + 'sizeof(%s);' % class_name
                        for class_name in self.class_names()
                       ]
        return '\n'.join(sizeof_lines)


class TemplateData(list):

    def __init__(self):
        pass

    def __str__(self):
        return '\n'.join([t.full_name() for t in self])

    def full_str(self):
        return '\n'.join([str(t) for t in self])

    def class_names(self, use_namespace=True):
        for t in self:
            for cls in t.class_names(use_namespace=use_namespace):
                yield cls

    def add(self, template_spec):
        self.append(template_spec)

    def typedef_code(self, indent=0):
        return '\n'.join([t.typedef_code(indent=indent) for t in self])

    def sizeof_code(self, indent=0):
        return '\n'.join([t.sizeof_code(indent=indent) for t in self])


def locate(pattern, root=os.curdir):
    '''
    Locate all files matching supplied filename pattern
    in and below supplied root directory.
    '''

    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


# exec open('filename').read() 

#
# settings
#

# list of header files that will be "seen" by the wrapper generator
header_files_list = ['faunus/faunus.h',
                     'faunus/common.h',
                     'faunus/point.h',
                     'faunus/particles.h',
                     'faunus/slump.h',
                     'faunus/inputfile.h',
                     'faunus/species.h',
                     'faunus/io.h',
                     'faunus/container.h',
                     'faunus/potentials/pot_coulomb.h',
                     'faunus/titrate.h',
                     'faunus/moves/base.h',
                     'faunus/moves/translational.h',
                     'faunus/moves/charge.h',
                     'faunus/group.h',
                     'faunus/analysis.h',
                     'faunus/widom.h',
                     'faunus/energy/base.h',
                     'faunus/ensemble.h',
                     'faunus/histogram.h',
                     'faunus/xytable.h',
                     'faunus/mcloop.h',
                     'faunus/average.h',
                    ]

header_files = header_files_list
#header_files = list(locate('*.h', '../../include/faunus'))

# class wrapper requests
classes = ['point',                             # point.h
           'particle',
           'spherical',
           'particles',                         # particles.h
           'random',                            # slump.h
           'randomDefault',
           'randomTwister',
           'inputfile',                         # inputfile.h
           'atoms',                             # species.h
           'container',                         # container.h
           'cell',
           'box',
           'slit',
           'clutch',
           'cylinder',
           'ensemble',                          # ensemble.h
           'canonical',
           'grandcanonical',
           'pot_coulomb',                       # potentials/pot_coulomb.h
           'pot_lj',                            # potentials/base.h
           'energybase',                        # energy/base.h
           'markovmove',                        # moves
           'saltmove',
           'titrate',
           'chargereg',
           'group',                             # group.h
           'macromolecule',
           'salt',
           'analysis',                          # analysis.h
           'virial',
           'systemenergy',
           'widom',                             # widom.h
           'widomSW',
           'mcloop',                            # mcloop.h
           'ioaam',                             # io.h
           'iopqr',
           'ioxyz',
           ]

other_stuff = ['faunus_splash']

# template instantiation requests
# dictionary with the fillowing structure:
# - key - C++ template
# - value - a tuple of of two-tuples, each of them contains template arguments and a python alias
# also include here any non-Faunus templates that you want instantiated and/or aliased (e. g. vector<double>)

template_data = TemplateData()

td = TemplateSpec('interaction', 'Faunus')
td.add_instantiation('Faunus::pot_coulomb', 'interaction_coulomb')
td.add_instantiation('Faunus::pot_hscoulomb', 'interaction_hscoulomb')
template_data.append(td)

td = TemplateSpec('average', 'Faunus')
td.add_instantiation('float', 'average_float')
td.add_instantiation('double', 'average_double')
template_data.append(td)

td = TemplateSpec('vector', 'std')
#td.add_instantiation('float', 'vector_float')
#td.add_instantiation('double', 'vector_double')
#template_data.append(td)

del td

print template_data.full_str()


#
# generator code
#

# documentation says it is safer to use absolute paths for the header files
header_files_abs = [os.path.abspath(p) for p in header_files]

generated_header_template = """%s

namespace pyplusplus{ namespace aliases{

%s

    inline void instantiate() {
%s
    }
} }
"""

def generate_header(header_files, template_data_list):
    """
    template_data is a dictionary as explained above
    """

    headers = '\n'.join(['#include <%s>' % header for header in header_files])

    return generated_header_template % (headers, template_data.typedef_code(indent=4), template_data.sizeof_code(indent=8))


# generate the header file to be passed to gccxml
header_code = generate_header(header_files, template_data)

# write the generated header file to disk
# (needed to get includes in the code generated by py++ right)
try:
    os.mkdir('generated')
except OSError:
    pass
out = file('generated/generated_header.h', 'w')
out.write(header_code)
out.close()

# parse the header file
mb = module_builder_t(files=['generated/generated_header.h'],
                      include_paths=[faunus_base + '/include/'],
                      gccxml_path='../auxiliary/bin/gccxml',
                      indexing_suite_version=1
                     )

# by default, do not expose anything
mb.decls().exclude()

# prepare custom lists
decls_to_include = []
decls_to_exclude = []

# include all the requested classes
decls_to_include.extend([mb.class_(cls) for cls in classes])

# include other requested stuff
decls_to_include.extend([mb.decl(decl) for decl in other_stuff])

# also include all the instantiated templates
decls_to_include.extend([mb.decl(template) for template in template_data.class_names(use_namespace=False)])

# print all the declaration to be included
print
print 'Declarations to be included:'
for d in decls_to_include:
    print d
print

# apply the requested includes and excludes
for decl in decls_to_include:
   decl.include()
for decl in decls_to_exclude:
    decl.exclude()

# do not expose private and protected members
mb.calldefs(access_type_matcher_t('private')).exclude()
mb.calldefs(access_type_matcher_t('protected')).exclude()

# run the generator
mb.build_code_creator(module_name='faunus')

# write source files to disk
mb.split_module('./generated')
