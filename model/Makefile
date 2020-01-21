CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)
CMSSW_INC = $(CMSSW_BASE)
REPO = $(shell git rev-parse --show-toplevel)/model
COMMON = $(shell git rev-parse --show-toplevel)/common

CPPFLAGS := $(shell root-config --cflags) -I$(REPO)/include -I$(COMMON)/include -I$(CMSSW_INC)/src
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)  -lRooFit -lRooFitCore -L$(CMSSW_INC)/lib/slc7_amd64_gcc700/ -lHiggsAnalysisCombinedLimit

CPPFLAGS += -g -std=c++1z

TARGET = MakeWorkspace
TARGET2 = MakeWorkspace_VR
TARGET3 = MakeWorkspace_DifferentYears

SRC = $(REPO)/build_ws.cc $(COMMON)/src/CommandLineInput.cc
SRC2 = $(REPO)/app/build_ws_validation_region.cc $(COMMON)/src/CommandLineInput.cc
SRC3 = $(REPO)/app/build_ws_combine_years.cc $(COMMON)/src/CommandLineInput.cc

OBJ = $(SRC:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)
OBJ3 = $(SRC3:.cc=.o)

all : $(TARGET) $(TARGET2) $(TARGET3)

$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET2) : $(OBJ2)
	$(LD) $(CPPFLAGS) -o $(TARGET2) $(OBJ2) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET3) : $(OBJ3)
	$(LD) $(CPPFLAGS) -o $(TARGET3) $(OBJ3) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@
	@echo $<

clean :
	rm -f *.o app/*.o src/*.o include/*.o $(TARGET) $(TARGET2) $(TARGET3) $(COMMON)/src/*.o *~
