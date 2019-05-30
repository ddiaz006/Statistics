CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)
CMSSW_INC = $(CMSSW_BASE)
REPO = $(shell git rev-parse --show-toplevel)/model
COMMON = $(shell git rev-parse --show-toplevel)/common

CPPFLAGS := $(shell root-config --cflags) -I$(REPO)/include -I$(COMMON)/include -I$(CMSSW_INC)/src
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)  -lRooFit -lRooFitCore -L$(CMSSW_INC)/lib/slc6_amd64_gcc530/ -lHiggsAnalysisCombinedLimit

CPPFLAGS += -g -std=c++1y

TARGET = MakeWorkspace

SRC = $(REPO)/build_ws.cc $(COMMON)/src/CommandLineInput.cc

OBJ = $(SRC:.cc=.o)

all : $(TARGET) 

$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@
	@echo $<
clean :
	rm -f *.o src/*.o include/*.o $(TARGET) $(COMMON)/src/*.o *~