all: 
CXX:=`root-config --cxx`
SNSW=/junofs/users/lihl/neutrino/superN/SNsim/simulation

CPPFLAGS=-D_PGTRACK_  -I$(SNSW)/include -I$(ROOUNFOLD)/src -isystem$(ROOTSYS)/include -I$(ROOTSYS)/include 
CXXFLAGS= -ggdb  $(CPPFLAGS)  

EXTRALIBS = -lSNsim -lRooUnfold -lMathMore -lMinuit -lHistPainter 
EXTRALIBS += `root-config --libs` -lc 

ifeq "$(CXX)" "g++"
	EXTRALIBS += -lgfortran
endif


%:%.cxx 
	$(CXX) -Wall -o $@ $@.cxx -L/lib64 -L$(SNSW)/lib -L$(ROOUNFOLD) $(CXXFLAGS)  $(EXTRALIBS)
clean:
	@rm -rvf unfoldtest
