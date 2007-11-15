  MODEL = generic
  GROMACS = no

  ##################################################################

CLASSDIR=./classes
INCDIR=-I$(CLASSDIR)

ifeq ($(GROMACS), yes)
  INCDIR=-I/usr/local/gromacs/include/gromacs/ -I$(CLASSDIR)
  LDFLAGS=-L/usr/local/gromacs/lib/ -lgmx
  GROFLAG=-DGROMACS
endif

ifeq ($(MODEL), debug)
  CXX=g++
  CXXFLAGS = -O0 -w -Winline  -g $(INCDIR) $(GROFLAG)
endif

ifeq ($(MODEL), generic)
  CXX=g++
  CXXFLAGS = -O3 -funroll-loops -w -Winline -g $(INCDIR) $(GROFLAG)
endif

OBJS=$(CLASSDIR)/inputfile.o \
     $(CLASSDIR)/io.o\
     $(CLASSDIR)/titrate.o\
     $(CLASSDIR)/point.o \
     $(CLASSDIR)/physconst.o\
     $(CLASSDIR)/potentials.o\
     $(CLASSDIR)/slump.o\
     $(CLASSDIR)/container.o\
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

