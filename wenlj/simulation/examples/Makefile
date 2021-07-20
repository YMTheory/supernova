all: etSpec2dVis
CXX:=`root-config --cxx`

SNSW=$SNsim/simualtion
CPPFLAGS=-D_PGTRACK_  -I$(SNSW)/include -I$(ROOTSYS)/include 
CXXFLAGS= -ggdb  $(CPPFLAGS)  

EXTRALIBS = -lSNsim  
EXTRALIBS += `root-config --libs`

ifeq "$(CXX)" "g++"
	EXTRALIBS += -lgfortran
endif


%:%.cxx 
	$(CXX) -g -o $@ $@.cxx -L$(SNSW)/lib  $(CXXFLAGS)  $(EXTRALIBS)
clean:
	@rm -rvf etSpec2dVis
