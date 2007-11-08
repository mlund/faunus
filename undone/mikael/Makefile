  #
  # CONFIGURATION
  #
  # Uncomment line to match your system:

  #POTENTIAL = POT_COULOMB

  #MODEL = intel
  #MODEL = ibmxl
  #MODEL = debug
  #MODEL = mac
  MODEL = generic

  ##################################################################

  INCFLAG = -I../classes

  CXX=g++
  ifeq ($(MODEL), mac)
    CXXFLAGS = -O3 -funroll-loops \
	       -Wall -Wno-sign-compare -Wno-unused \
	       -mtune=G5 -mcpu=G5 -mpowerpc64 \
	       -mpowerpc-gpopt -force_cpusubtype_ALL
  endif
  ifeq ($(MODEL), debug)
    CXX=g++
    CXXFLAGS = $(INCFLAG) -g -O0 -w -Wno-sign-compare -Wno-unused
  endif
  ifeq ($(MODEL), generic)
    CXX=g++
    CXXFLAGS = $(INCFLAG) -O0 -funroll-loops -Winline
  endif
  ifeq ($(MODEL), ibmxl)
    CXX=xlC
    CXXFLAGS = $(INCFLAG) -O4 -DXLC
  endif
  ifeq ($(MODEL), intel)
    CXX=icc
    CXXFLAGS = $(INCFLAG) -fast -w -openmp -parallel -D$(POTENTIAL)
endif

CLASSDIR=../classes
INCDIR=./

OBJS=$(CLASSDIR)/space.o \
     $(CLASSDIR)/group.o \
     $(CLASSDIR)/inputfile.o \
     $(CLASSDIR)/montecarlo.o \
     $(CLASSDIR)/io.o\
     $(CLASSDIR)/titrate.o\
     $(CLASSDIR)/point.o \
     $(CLASSDIR)/povray.o \
     $(CLASSDIR)/physconst.o\
     $(CLASSDIR)/potentials.o\
     $(CLASSDIR)/slump.o\
     $(CLASSDIR)/container.o\
     $(CLASSDIR)/hardsphere.o\
     $(CLASSDIR)/species.o
all:	prot_md 

prot_md:	$(OBJS) prot_md.C
		$(CXX) $(CXXFLAGS) $(OBJS) prot_md.C -o prot_md

volume:      $(OBJS) volume.C
	        $(CXX) $(CXXFLAGS) $(OBJS) volume.C -o volume

surface:	$(OBJS) surface.C
	$(CXX) $(CXXFLAGS) $(OBJS) surface.C -o surface

doc:
	doxygen

clean:
	rm -f $(OBJS) $(OBJSN) prot_md 

