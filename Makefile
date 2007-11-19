# Specify compiler (gnu,intel,pgi,pathscale,debug)
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
CLASSDIR=./classes
INCDIR=-I$(CLASSDIR)

ifeq ($(GROMACS), yes)
  INCDIR=-I/usr/local/gromacs/include/gromacs/ -I$(CLASSDIR)
  LDFLAGS=-L/usr/local/gromacs/lib/ -lgmx
  GRO=-DGROMACS
endif

ifeq ($(MODEL), debug)
  CXXFLAGS = -O0 -Wextra -Winline -Wno-sign-compare -g $(INCDIR) $(GRO)
endif

ifeq ($(MODEL), gnu)
  ifeq ($(OPENMP), yes)
    OMP=-fopenmp
  endif
  CXXFLAGS = -O3 -w -funroll-loops $(INCDIR) $(GRO) $(OMP)
endif

ifeq ($(MODEL), intel)
  ifeq ($(OPENMP), yes)
    OMP=-openmp
  endif
  CXX=icc
  CXXFLAGS = -O3 -w $(INCDIR) $(GRO) $(OMP)
endif

ifeq ($(MODEL), pathscale)
  CXX=pathCC
  CXXFLAGS = -Ofast $(INCDIR) $(GRO)
endif

ifeq ($(MODEL), pgi)
  CXX=pgCC
  ifeq ($(OPENMP), yes)
    OMP=-mp
  endif 
  CXXFLAGS = -O3 $(INCDIR) $(GRO) $(OMP)
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
all:	classes examples

classes:	$(OBJS)
manual:
	doxygen doc/Doxyfile

widom:	examples/widom/widom.C $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) examples/widom/widom.C -o examples/widom/widom
	
ewald:	examples/ewald/ewald.C $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) examples/ewald/ewald.C -o examples/ewald/ewald

twobody:	examples/twobody/twobody.C $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) examples/twobody/twobody.C -o examples/twobody/twobody

pka:	examples/titration/pka.C $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(INCDIR) examples/titration/pka.C -o examples/titration/pka

examples:	widom pka ewald twobody

clean:
	rm -vf $(OBJS) examples/titration/pka examples/widom/widom examples/ewald/ewald

docclean:
	rm -vfR doc/html doc/latex

