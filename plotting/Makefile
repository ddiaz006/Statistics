CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)
CMSSW_INC = $(CMSSW_BASE)
REPO = $(shell git rev-parse --show-toplevel)/plotting
COMMON = $(shell git rev-parse --show-toplevel)/common

CPPFLAGS := $(shell root-config --cflags) -I$(REPO)/include -I$(COMMON)/include -I$(CMSSW_INC)/src
#LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)  -lRooFit -lRooFitCore -L$(CMSSW_INC)/lib/slc6_amd64_gcc530/ -lHiggsAnalysisCombinedLimit
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)  -lRooFit -lRooFitCore -L$(CMSSW_INC)/lib/slc6_amd64_gcc530/

CPPFLAGS += -g -std=c++1y -fsanitize=address -Wall -Wextra -Wno-sign-compare

TARGET = PlotLimits
TARGET2 = PlotVR
TARGET6 = PlotVR_LL
TARGET3 = GetExtraSysUnc
TARGET4 = GetZ_TF_InclusiveCorrection
TARGET5 = Plot_TF

SRC = $(REPO)/app/plot_limits.cc $(COMMON)/src/CommandLineInput.cc
SRC2 = $(REPO)/app/plot_validation_region.cc $(REPO)/src/helper_functions.cc $(COMMON)/src/CommandLineInput.cc
SRC6 = $(REPO)/app/plot_validation_region_ll.cc $(REPO)/src/helper_functions.cc $(COMMON)/src/CommandLineInput.cc
SRC3 = $(REPO)/app/get_extra_sys_unc.cc $(REPO)/src/helper_functions.cc $(COMMON)/src/CommandLineInput.cc
SRC4 = $(REPO)/app/get_Z_tf_inclusive_correction.cc $(REPO)/src/helper_functions.cc $(COMMON)/src/CommandLineInput.cc
SRC5 = $(REPO)/app/plot_tf.cc $(REPO)/src/helper_functions.cc $(COMMON)/src/CommandLineInput.cc

OBJ = $(SRC:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)
OBJ3 = $(SRC3:.cc=.o)
OBJ4 = $(SRC4:.cc=.o)
OBJ5 = $(SRC5:.cc=.o)
OBJ6 = $(SRC6:.cc=.o)

all : $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6)

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

%.o : %.cc
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@
	@echo $<

clean :
	rm -f *.o app/*.o src/*.o include/*.o $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) *~
