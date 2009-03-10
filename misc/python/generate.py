#!/usr/bin/env python

from IPython.Shell import IPShellEmbed

ipshell = IPShellEmbed([])

try:
    from pyplusplus import module_builder
except:
    import sys
    print "pyplusplus import failed, aborting"
    sys.exit(1)

# TODO: get gccxml path from the configuration/cmake
# TODO: a proper way to treat include search directories?

# parse the header files
mb = module_builder.module_builder_t(files=['faunus_export.h'],
                                     include_paths=['/home/andy/code/faunus/trunk/include/'],
                                     indexing_suite_version=2
                                    )

# pick what gets exposed
#TODO: get classes from Faunus::, perhaps not all of them
#      get some stl? how much?

# by default, do not expose anything
mb.decls().exclude()

# will that be it in the end?
#to_expose.append(mb.classes(header_dir='/home/andy/code/faunus/trunk/include'))

ipshell()

# point.h
mb.class_('point').include()
mb.class_('particle').include()
mb.class_('spherical').include()

# particles.h
mb.class_('particles').include()

# slump.h
mb.class_('random').include()
mb.class_('randomDefault').include()
mb.class_('randomTwister').include()

# inputfile.h
mb.class_('inputfile').include()

# species.h
mb.class_('atoms').include()

# container.h
mb.class_('container').include()
mb.class_('cell').include()
mb.class_('box').include()
mb.class_('slit').include()
mb.class_('clutch').include()
mb.class_('cylinder').include()

# ensemble.h
mb.class_('ensemble').include()
mb.class_('canonical').include()

# potentials
mb.class_('pot_coulomb').include()
mb.class_('pot_lj').include()

# energy.h
mb.class_('energybase').include()
mb.class_('interaction<Faunus::pot_coulomb>').include()
mb.decls(lambda decl: 'interaction' in decl.name).include()

# moves
mb.class_('markovmove').include()
mb.class_('saltmove').include()

# group.h
mb.class_('group').include()

# analysis.h
mb.class_('analysis').include()
mb.class_('systemenergy').include()

# widom.h
mb.class_('widom').include()

# take care of other settings
# (nothing at the moment)

# generate the code
mb.build_code_creator(module_name='faunus')
mb.split_module('./generated')

