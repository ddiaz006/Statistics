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
#TARGET3 = MakeWorkspace_DifferentYears
TARGET4 = MakeWorkspace_VR_SYS
TARGET5 = MakeMockupSystematics
TARGET6 = MakeWorkspace_VR_Merged_EMu
TARGET7 = MakeWorkspace_VR_SYS_ALL_TF

SRC = $(REPO)/build_ws.cc $(COMMON)/src/CommandLineInput.cc
SRC2 = $(REPO)/app/build_ws_validation_region.cc $(COMMON)/src/CommandLineInput.cc
#SRC3 = $(REPO)/app/build_ws_combine_years.cc $(COMMON)/src/CommandLineInput.cc
SRC4 = $(REPO)/app/build_ws_validation_region_systematics.cc $(COMMON)/src/CommandLineInput.cc
SRC5 = $(REPO)/app/make_mockup_systematics.cc $(COMMON)/src/CommandLineInput.cc
SRC6 = $(REPO)/app/build_ws_validation_region_merged_emu.cc $(COMMON)/src/CommandLineInput.cc
SRC7 = $(REPO)/app/build_ws_validation_region_systematics_all_tf.cc $(COMMON)/src/CommandLineInput.cc

OBJ = $(SRC:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)
#OBJ3 = $(SRC3:.cc=.o)
OBJ4 = $(SRC4:.cc=.o)
OBJ5 = $(SRC5:.cc=.o)
OBJ6 = $(SRC6:.cc=.o)
OBJ7 = $(SRC7:.cc=.o)

all : $(TARGET) $(TARGET2) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) #$(TARGET3)

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

#$(TARGET3) : $(OBJ3)
#	$(LD) $(CPPFLAGS) -o $(TARGET3) $(OBJ3) $(LDFLAGS)
#	@echo $@
#	@echo $<
#	@echo $^


$(TARGET4) : $(OBJ4)
	$(LD) $(CPPFLAGS) -o $(TARGET4) $(OBJ4) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET5) : $(OBJ5)
	$(LD) $(CPPFLAGS) -o $(TARGET5) $(OBJ5) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET6) : $(OBJ6)
	$(LD) $(CPPFLAGS) -o $(TARGET6) $(OBJ6) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET7) : $(OBJ7)
	$(LD) $(CPPFLAGS) -o $(TARGET7) $(OBJ7) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@
	@echo $<

clean :
	rm -f *.o app/*.o src/*.o include/*.o $(TARGET) $(TARGET2) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(COMMON)/src/*.o *~ #$(TARGET3)
