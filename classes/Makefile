  MODEL = generic

  ##################################################################

  ifeq ($(MODEL), generic)
    CXX=g++
    CXXFLAGS = -ftree-vectorize -O3 -funroll-loops -w -Winline
  endif

CLASSDIR=.
INCDIR=.

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
all:	classes 

classes:	$(OBJS)
doc:	
	doxygen doc/Doxyfile

widom:	examples/widom-example.C $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) examples/widom-example.C -o examples/widom-example

pka:	examples/protein-example.C $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) examples/protein-example.C -o examples/protein-example

examples:	widom pka

clean:
	rm -f $(OBJS) examples/protein-example examples/widom-example 

docclean:
	rm -fR doc/html doc/latex

