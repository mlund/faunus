  #
  # CONFIGURATION
  #
  # Uncomment line to match your system:

  POTENTIAL = POT_COULOMB
  #POTENTIAL = POT_COULOMB

  #MODEL = pentium
  #MODEL = intel
  #MODEL = macg5
  #MODEL = macg4
  #MODEL = ibmxl
  #MODEL = debug
  #MODEL = mac
  MODEL = generic

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
    CXX=g++
    CXXFLAGS = -g -O0 -w -Wno-sign-compare -Wno-unused
  endif
  ifeq ($(MODEL), generic)
    CXX=g++
    CXXFLAGS = -g -O0 -w -funroll-loops -D$(POTENTIAL)
  endif
  ifeq ($(MODEL), ibmxl)
    CXX=xlC
    CXXFLAGS = -O4 -DXLC
  endif
  ifeq ($(MODEL), intel)
    CXX=icc
    CXXFLAGS = -fast -w -openmp -parallel -D$(POTENTIAL)
endif

CLASSDIR=./
INCDIR=./

OBJS=$(CLASSDIR)/space.o \
     $(CLASSDIR)/group.o \
     $(CLASSDIR)/point.o \
     $(CLASSDIR)/inputfile.o \
     $(CLASSDIR)/montecarlo.o \
     $(CLASSDIR)/povray.o \
     $(CLASSDIR)/physconst.o\
     $(CLASSDIR)/cell.o\
     $(CLASSDIR)/potentials.o\
     $(CLASSDIR)/slump.o\
     $(CLASSDIR)/species.o\
     $(CLASSDIR)/io.o\
     $(CLASSDIR)/titrate.o\
     $(CLASSDIR)/pdbclass.o\

all:	prot prot_md pka pka_chain pka_dielec prot_dh hairy_dh

prot:	$(OBJS) prot.C
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) prot.C -o prot

prot_md:	$(OBJS) prot_md.C
		$(CXX) $(CXXFLAGS) -I$(INCDIR) $(OBJS) prot_md.C -o prot_md

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

doc:
	doxygen

clean:
	rm -f $(OBJS) $(OBJSN) prot pka pka_dielec pka_chain

