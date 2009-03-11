#!/usr/bin/env python

#from IPython.Shell import IPShellEmbed
#ipshell = IPShellEmbed([])

try:
    from pyplusplus.module_builder import module_builder_t, create_text_fc
    from pygccxml.declarations.matchers import access_type_matcher_t
except:
    import sys
    print "import from pyplusplus or pygccxml failed, aborting"
    sys.exit(1)

# TODO: generate this string for all templates we want instantiated
code = """
#include <faunus/energy.h>
#include <faunus/potentials/pot_coulomb.h>

namespace pyplusplus{ namespace aliases{

    typedef Faunus::interaction<Faunus::pot_coulomb> interaction_coulomb;

    inline void instantiate() {
        sizeof(Faunus::interaction<Faunus::pot_coulomb>);
    }

} }
"""

header_files = [create_text_fc(code),
                'faunus/common.h',
                'faunus/point.h',
                'faunus/particles.h',
                'faunus/slump.h',
                'faunus/inputfile.h',
                'faunus/species.h',
                'faunus/io.h',
                'faunus/container.h',
                'faunus/potentials/pot_coulomb.h',
                'faunus/moves/base.h',
                'faunus/moves/translational.h',
                'faunus/group.h',
                'faunus/analysis.h',
                'faunus/widom.h',
                'faunus/energy.h',
                'faunus/ensemble.h',
                'faunus/histogram.h',
                'faunus/xytable.h',
               ]

# parse the header files
mb = module_builder_t(files=header_files,
                      include_paths=['../../include/'],
                      #indexing_suite_version=1
                     )

#
# pick what gets exposed
#

# by default, do not expose anything
mb.decls().exclude()

# will that be it in the end?
#to_expose.append(mb.classes(header_dir='/home/andy/code/faunus/trunk/include'))

classes = [
# point.h
'point',
'particle',
'spherical',
# particles.h
'particles',
# slump.h
'random',
'randomDefault',
'randomTwister',
# inputfile.h
'inputfile',
# species.h
'atoms',
# container.h
'container',
'cell',
'box',
'slit',
'clutch',
'cylinder',
# ensemble.h
'ensemble',
'canonical',
# potentials
'pot_coulomb',
'pot_lj',
# energy.h
'energybase',
'interaction<Faunus::pot_coulomb>',
# # moves
'markovmove',
'saltmove',
# # group.h
'group',
# analysis.h
'analysis',
'systemenergy',
# widom.h
'widom',
]

decls_to_include = []
decls_to_exclude = []

decls_to_include.extend([mb.class_(cls) for cls in classes])

# TODO: do we need this?
#mb.decls(lambda decl: 'interaction' in decl.name).include()

for decl in decls_to_include:
   decl.include()

for decl in decls_to_exclude:
    decl.exclude()

# do not expose private and protected members
mb.calldefs( access_type_matcher_t( 'private' ) ).exclude()
mb.calldefs( access_type_matcher_t( 'protected' ) ).exclude()

# run the generator
mb.build_code_creator(module_name='faunus')

# write source files to disk
# TODO: switch to split files as soon as moves are split into .h and .cpp
#mb.split_module('./generated')
mb.write_module('./generated/main.cpp')
