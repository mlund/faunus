  MODEL = generic

  ##################################################################

  ifeq ($(MODEL), generic)
    CXX=g++
    CXXFLAGS = -O3 -funroll-loops -w -Winline  -g
  endif

CLASSDIR=classes
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

pka:	pka.C $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) pka.C -o pka

clutch:	examples/clutch-example.C $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) examples/clutch-example.C -o examples/clutch-example

wall:   wall.C $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) wall.C -o wall

examples:	widom pka clutch

clean:
	rm -f $(OBJS) examples/protein-example examples/widom-example pka wall

docclean:
	rm -fR doc/html doc/latex

slumptest: slumptest.C $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) slumptest.C -o slumptest

