  #
  # CONFIGURATION
  #
  # Uncomment line to match your system:

  #MODEL = pentium
  MODEL = intel
  #MODEL = macg5
  #MODEL = macg4
  #MODEL = ibmxl
  #MODEL = debug
  #MODEL = mac
  #MODEL = generic
 #$(CLASSDIR)/imdwrap.o\

  ##################################################################

  CXX=g++
  ifeq ($(MODEL), mac)
    CXXFLAGS = -O3 -funroll-loops \
	       -Wall -Wno-sign-compare -Wno-unused \
	       -mtune=G5 -mcpu=G5 -mpowerpc64 \
	       -mpowerpc-gpopt -force_cpusubtype_ALL
  endif
  ifeq ($(MODEL), macg5)
    CXXFLAGS = -O3 -DMACOSX -funroll-loops \
	       -Wall -Wno-sign-compare -Wno-unused \
	       -fsched-interblock -falign-loops=32 -falign-jumps=16 \
	       -falign-functions=16 -falign-jumps-max-skip=15 \
	       -falign-loops-max-skip=15 -mdynamic-no-pic -mpowerpc-gpopt \
	       -force_cpusubtype_ALL -mtune=G5 -mcpu=G5 -mpowerpc64
  endif
  ifeq ($(MODEL), macg4)
    CXXFLAGS = -O3 -DMACOSX -funroll-loops \
	       -Wall -Wno-sign-compare -Wno-unused \
	       -fsched-interblock -falign-loops=32 -falign-jumps=16 \
	       -falign-functions=16 -falign-jumps-max-skip=15 \
	       -falign-loops-max-skip=15 -mdynamic-no-pic -mpowerpc-gpopt \
	       -force_cpusubtype_ALL -mtune=7450 -mcpu=7450 -mpowerpc
  endif
  ifeq ($(MODEL), pentium)
    CXX=g++
    CXXFLAGS = -O3 -W -funroll-loops \
	       -mcpu=pentium4 -march=pentium4 \
	       -msse2
  endif
  ifeq ($(MODEL), debug)
    CXXFLAGS = -g -O2 -w -Wno-sign-compare -Wno-unused
  endif
  ifeq ($(MODEL), generic)
    CXX=g++
    CXXFLAGS = -O2 -w -funroll-loops
  endif
  ifeq ($(MODEL), ibmxl)
    CXX=xlC
    CXXFLAGS = -O4 -DXLC
  endif
  ifeq ($(MODEL), intel)
    CXX=icc
    CXXFLAGS = -O3 -w #-openmp #-parallel -openmp
endif

CLASSDIR=./
INCDIR=./

OBJS=$(CLASSDIR)/space.o \
 $(CLASSDIR)/vpython.o \
 $(CLASSDIR)/interact.o \
 $(CLASSDIR)/group.o \
 $(CLASSDIR)/point.o \
 $(CLASSDIR)/inputfile.o \
 $(CLASSDIR)/montecarlo.o \
 $(CLASSDIR)/histogram.o \
 $(CLASSDIR)/povray.o \
 $(CLASSDIR)/physconst.o\
 $(CLASSDIR)/slump.o\
 $(CLASSDIR)/cell.o\
 $(CLASSDIR)/chain.o\
 $(CLASSDIR)/average.o\
 $(CLASSDIR)/titrate.o\
 $(CLASSDIR)/peptide.o\
 $(CLASSDIR)/widommod2.o\
 $(CLASSDIR)/pdbclass.o

 #$(CLASSDIR)/vmd/imd.o\

all:	prot pka pka_chain pka_dielec prot_dh hairy_dh

prot:	$(OBJS) prot.C
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) prot.C -o prot

prot_dh:	$(OBJS) prot_dh.C    
		$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) prot_dh.C -o prot_dh

hairy_dh:       $(OBJS) hairy_dh.C
	        $(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) hairy_dh.C -o hairy_dh

pka:	$(OBJS) pka.C
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) pka.C -o pka

pka_chain:	$(OBJS) pka_chain.C
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) pka_chain.C -o pka_chain 

pka_dielec:      $(OBJS) pka_dielec.C
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) pka_dielec.C -o pka_dielec

volume:      $(OBJS) volume.C
	        $(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) volume.C -o volume

surface:	$(OBJS) surface.C
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) surface.C -o surface

rfield_3d:	$(OBJS) examples/rfield_3d.C
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) examples/rfield_3d.C -o rfield_3d

pkaGC:		$(OBJS) pkaGC.C
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) pkaGC.C -o pkaGCnewpka

pkaIT:		$(OBJS) pkaIT.C
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) pkaIT.C -o pkaIT
doc:
	doxygen

clean:
	rm -f $(OBJS) prot pka pka_dielec pka_chain

