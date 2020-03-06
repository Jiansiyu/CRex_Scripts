









#------------------------------------------------------------------------------
# Core library
SRC = setup.C

EXTRAHDR = 

CORE = CRexReplay

LINKDEF = $(CORE)_LinkDef.h

#------------------------------------------------------------------------------
# Compile debug version (for gdb)
export DEBUG = 1
# Compile extra code for printing verbose messages (enabled with fDebug)
export VERBOSE = 1
# Compile extra diagnostic code (extra computations and global variables)
export TESTCODE = 1
# Compile support code for MC input data
export MCDATA = 1

#export I387MATH = 1
export EXTRAWARN = 1

# Architecture to compile for
MACHINE := $(shell uname -s)
ARCH    := linux
SOSUF   := so
ifeq ($(MACHINE),Darwin)
  ARCH := macosx
  SOSUF := dylib
endif

#------------------------------------------------------------------------------
# Directory locations. All we need to know is INCDIRS.
# INCDIRS lists the location(s) of the C++ Analyzer header (.h) files

ifndef ANALYZER
  $(error $$ANALYZER environment variable not defined)
endif

ifndef TREESEARCH
  $(error $$TREESEARCH environment variable not defined)
endif

INCDIRS  = $(wildcard $(addprefix $(ANALYZER)/, include src hana_decode hana_scaler Podd HallA) $(TREESEARCH)/ )

ifdef EVIO_INCDIR
  INCDIRS += ${EVIO_INCDIR}
else ifdef EVIO
  INCDIRS += ${EVIO}/include
endif

#------------------------------------------------------------------------------
# Do not change anything  below here unless you know what you are doing

ifeq ($(strip $(INCDIRS)),)
  $(error No Analyzer header files found. Check $$ANALYZER)
endif

ROOTVERMAJOR := $(shell root-config --version | cut -d. -f1)
ROOTVERMINOR := $(shell root-config --version | cut -d. -f1 | cut -d/ -f1)
ROOTVERPATCH := $(shell root-config --version | cut -d. -f1 | cut -d/ -f2)
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTBIN      := $(shell root-config --bindir)
CXX          := $(shell root-config --cxx)
LD           := $(shell root-config --ld)

## ROOT5 and under used rootcint to make dictionaries
ifeq ($(ROOTVERMAJOR),5)
	ROOTDICT_CMD = $(ROOTBIN)/rootcint
	ROOTDICT_CMD_FLAGS = 
else
## ROOT6 uses rootcling
	ROOTDICT_CMD = $(ROOTBIN)/rootcling
	ROOTDICT_CMD_FLAGS = -rmf $(COREDICT).rootmap
endif

PKGINCLUDES  = $(addprefix -I, $(INCDIRS) ) -I$(shell pwd)
INCLUDES     = -I$(shell root-config --incdir) $(PKGINCLUDES)

CORELIB  = lib$(CORE).$(SOSUF)
COREDICT = $(CORE)Dict

LIBS          = 
GLIBS         = 

ifeq ($(ARCH),linux)
# Linux with gcc (RedHat)
ifdef DEBUG
  CXXFLAGS    = -g -O0
  LDFLAGS     = -g -O0
  DEFINES     =
else
  CXXFLAGS    = -O2 -g #-march=pentium4
  LDFLAGS     = -O -g
  DEFINES     = -DNDEBUG
endif
DEFINES      += -DLINUXVERS -DHAS_SSTREAM
CXXFLAGS     += -Wall -Woverloaded-virtual -fPIC
DICTCXXFLG   :=
ifdef EXTRAWARN
#FIXME: should be configure'd:
CXXVER       := $(shell g++ --version | head -1 | sed 's/.* \([0-9]\)\..*/\1/')
ifeq ($(CXXVER),4)
CXXFLAGS     += -Wextra -Wno-missing-field-initializers
DICTCXXFLG   := -Wno-strict-aliasing 
endif
endif
SOFLAGS       = -shared
ifdef I387MATH
CXXFLAGS     += -mfpmath=387
else
CXXFLAGS     += -march=core2 -mfpmath=sse
endif
endif

ifeq ($(ARCH),macosx)
# Mac OS X with gcc >= 3.x or clang++ >= 5
ifdef DEBUG
  CXXFLG     := -g -O0
  LDFLAGS    := -g -O0
  DEFINES    :=
else
  CXXFLG     := -O
  LDFLAGS    := -O
  DEFINES    := -DNDEBUG
endif
DEFINES      += -DMACVERS -DHAS_SSTREAM
CXXFLG       += -Wall -fPIC
CXXEXTFLG     =
LD           := $(CXX)
LDCONFIG     :=
SOFLAGS      := -shared -Wl,-undefined,dynamic_lookup
SONAME       := -Wl,-install_name,
ifeq ($(CXX),clang++)
CXXEXTFLG    += -Wextra -Wno-missing-field-initializers -Wno-unused-parameter
else
#FIXME: should be configure'd:
CXXVER       := $(shell g++ --version | head -1 | sed 's/.* \([0-9]\)\..*/\1/')
ifeq ($(CXXVER),4)
CXXEXTFLG    += -Wextra -Wno-missing-field-initializers
DICTCXXFLG   := -Wno-strict-aliasing
endif
endif
endif

ifdef VERBOSE
DEFINES      += -DVERBOSE
endif

ifdef TESTCODE
DEFINES      += -DTESTCODE
endif
ifdef MCDATA
DEFINES      += -DMCDATA
endif

CXXFLAGS     += $(DEFINES) $(ROOTCFLAGS) $(ROOTCFLAGS) $(PKGINCLUDES)
LIBS         += $(ROOTLIBS)
GLIBS        += $(ROOTGLIBS)

PHONY: test



test:
	@echo ${LIBS}
	@echo ${INCLUDES}