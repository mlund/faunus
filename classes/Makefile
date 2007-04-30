  MODEL = generic

  ##################################################################

  ifeq ($(MODEL), generic)
    CXX=g++
    CXXFLAGS = -O0 -Winline
  endif

CLASSDIR=./
INCDIR=./

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
     $(CLASSDIR)/markovmove.o \
     $(CLASSDIR)/species.o
all:	classes 

classes:	$(OBJS)
doc:
	doxygen
clean:
	rm -f $(OBJS) 

