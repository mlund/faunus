  MODEL = generic

  ##################################################################

  ifeq ($(MODEL), generic)
    CXX=g++
    CXXFLAGS = -g -O3 -Winline
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
     $(CLASSDIR)/species.o
all:	classes 

classes:	$(OBJS)
doc:	
	doxygen doc/Doxyfile
examples:	$(OBJS) examples/markovmove-example.C
		$(CXX) $(CXXFLAGS) $(OBJS) examples/markovmove-example.C -o examples/markovmove-example

clean:
	rm -f $(OBJS) 

