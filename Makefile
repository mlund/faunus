# Specify compiler (gnu,gnu686,intel,sun,pgi,pathscale,debug)
MODEL = gnu

# Set to yes if you need Gromacs xtc file support
# (requires a working Gromacs installation)
GROMACS = no

# Set to "yes" to enable parallel execution on multi-core
# CPU's. OpenMP must be supported by the compiler.
OPENMP = no

###########################################
#  Normally you would not want to modify  #
#  things beyond this point.              #
###########################################

CXX=g++
CLASSDIR=classes
INCDIR=-I$(CLASSDIR)
LDFLAGS=-L./lib -lfaunus

ifeq ($(GROMACS), yes)
  INCDIR:=${INCDIR} -I/usr/local/gromacs/include/gromacs
  LDFLAGS:=${LDFLAGS} -L/usr/local/gromacs/lib/ -lgmx
  EXTRA=-DGROMACS
endif

ifeq ($(MODEL), debug)
  CXXFLAGS = -O0 -W -Winline -Wno-sign-compare -g $(INCDIR) $(EXTRA)
endif

ifeq ($(MODEL), gnu)
  ifeq ($(OPENMP), yes)
    EXTRA:=$(EXTRA) -fopenmp
  endif
  CXXFLAGS = -O3 -w -funroll-loops $(INCDIR) $(EXTRA)
endif

ifeq ($(MODEL), gnu686)
  ifeq ($(OPENMP), yes)
    EXTRA:=$(EXTRA) -fopenmp
  endif
  CXXFLAGS = -mtune=i686 -msse3 -O3 -w -funroll-loops $(INCDIR) $(EXTRA)
endif

ifeq ($(MODEL), intel)
  ifeq ($(OPENMP), yes)
    EXTRA:=$(EXTRA) -openmp
  endif
  CXX=icc
  CXXFLAGS = -O3 -w $(INCDIR) $(EXTRA)
endif

ifeq ($(MODEL), pathscale)
  CXX=pathCC
  CXXFLAGS = -Ofast $(INCDIR) $(EXTRA)
endif

ifeq ($(MODEL), pgi)
  CXX=pgCC
  ifeq ($(OPENMP), yes)
    EXTRA:=$(EXTRA) -mp
  endif 
  CXXFLAGS = -O3 -w $(INCDIR) $(EXTRA)
endif

ifeq ($(MODEL), sun)
  ifeq ($(OPENMP), yes)
    EXTRA:=$(EXTRA) -xopenmp
  endif
  CXX=sunCC
  CXXFLAGS = -fast -w $(INCDIR) $(EXTRA)
endif
	    
OBJS=$(CLASSDIR)/inputfile.o \
     $(CLASSDIR)/io.o\
     $(CLASSDIR)/titrate.o\
     $(CLASSDIR)/point.o \
     $(CLASSDIR)/physconst.o\
     $(CLASSDIR)/slump.o\
     $(CLASSDIR)/container.o\
     $(CLASSDIR)/potentials.o\
     $(CLASSDIR)/hardsphere.o\
     $(CLASSDIR)/group.o \
     $(CLASSDIR)/particles.o \
     $(CLASSDIR)/analysis.o \
     $(CLASSDIR)/species.o 
all:	classes examples libfaunus

classes:	$(OBJS)
libfaunus:      $(OBJS)
	ar cr lib/libfaunus.a $(OBJS)

#$(CXX) $(CXXFLAGS) -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o lib/libfaunus.A.dylib -dynamiclib $(OBJS)

manual:
	doxygen doc/Doxyfile
manualul:
	scp -rC doc/html/* mikaellund@shell.sourceforge.net:/home/groups/f/fa/faunus/htdocs/
widom:	examples/widom/widom.C libfaunus
	$(CXX) $(CXXFLAGS) examples/widom/widom.C -o examples/widom/widom ${LDFLAGS}
ewald:		examples/ewald/ewald.C libfaunus
	$(CXX) $(CXXFLAGS) examples/ewald/ewald.C -o examples/ewald/ewald ${LDFLAGS}
twobody:	examples/twobody/twobody.C libfaunus
	$(CXX) $(CXXFLAGS) examples/twobody/twobody.C -o examples/twobody/twobody ${LDFLAGS}
twobody-hof:	examples/twobody-hofmeister/twobody-hof.C libfaunus
	$(CXX) $(CXXFLAGS) examples/twobody-hofmeister/twobody-hof.C -o examples/twobody-hofmeister/twobody-hof ${LDFLAGS}
manybody:	examples/manybody/manybody.C libfaunus 
	$(CXX) $(CXXFLAGS) examples/manybody/manybody.C -o examples/manybody/manybody ${LDFLAGS}
manybodyclust:	examples/manybody/manybodyclust.C libfaunus 
	$(CXX) $(CXXFLAGS) examples/manybody/manybodyclust.C -o examples/manybody/manybodyclust ${LDFLAGS}
isobaric:	examples/isobaric/isobaric.C libfaunus 
	$(CXX) $(CXXFLAGS) examples/isobaric/isobaric.C -o examples/isobaric/isobaric ${LDFLAGS}
tools:	examples/tools/printpotential.C examples/tools/aam2pqr.C libfaunus 
	$(CXX) $(CXXFLAGS) examples/tools/printpotential.C -o examples/tools/printpotential ${LDFLAGS}
	$(CXX) $(CXXFLAGS) examples/tools/aam2pqr.C -o examples/tools/aam2prq ${LDFLAGS}
pka:	examples/titration/pka.C libfaunus
	$(CXX) $(CXXFLAGS) examples/titration/pka.C -o examples/titration/pka $(LDFLAGS) 
GCpka:	examples/titration/GCpka.C libfaunus
	$(CXX) $(CXXFLAGS) examples/titration/GCpka.C -o examples/titration/GCpka $(LDFLAGS) 
binding:	examples/binding/binding.C libfaunus
	$(CXX) $(CXXFLAGS) examples/binding/binding.C -o examples/binding/binding $(LDFLAGS) 
undone:		undone/mikael/namespace.C libfaunus
	$(CXX) $(CXXFLAGS) \
	undone/mikael/namespace.C \
	-o undone/mikael/namespace ${LDFLAGS}
examples:	binding tools widom pka GCpka ewald twobody twobody-hof manybody isobaric
clean:
	rm -f $(OBJS) *.o \
	examples/titration/pka \
	examples/widom/widom \
	examples/ewald/ewald \
	examples/twobody/twobody \
	examples/manybody/manybody \
	examples/manybody/manybodyclust \
	examples/isobaric/isobaric \
	examples/twobody-hofmeister/twobody-hof \
	lib/libfaunus.a
docclean:
	rm -fR doc/html doc/latex
babel:
	#curl -L -o openbabel-2.1.1.tar.gz http://downloads.sourceforge.net/openbabel/openbabel-2.1.1.tar.gz
	#tar -zxf openbabel-2.1.1.tar.gz
	cd openbabel-2.1.1 ; ./configure ;make
	
