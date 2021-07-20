#-----------------------------------
#Created by Huiling, 2016.08.18
#-----------------------------------



#================COMPILER================#
CXX:=`root-config --cxx`

#=================FILES===================#
DSRC=$(CURDIR)/src/
DINC=$(CURDIR)/include/
DOBJ=$(CURDIR)/obj/
DLIB=$(CURDIR)/lib/
DICTDIR=$(CURDIR)/rootdict/
SRC=$(wildcard $(DSRC)*.cxx)
HEADERS=$(wildcard $(DINC)*.hh)
OBJS=$(patsubst $(DSRC)%.cxx, $(DOBJ)%.o, $(SRC))
MYLIB=libSNsim
LIBSO=$(DLIB)$(MYLIB).so
LIBA=$(DLIB)$(MYLIB).a
ROOTDICT=$(MYLIB)_dict
#SRC+=$(DSRC)$(ROOTDICT).cxx

#====FIXME: Set the ROOTSYS environment=============#
ROOTSYS=/home/liangjian/SoftWare/bin/ROOT
#====FIXME: NO NEED to set ROOTSYS if use JUNO environment=====#

#=============COMPILING FLAG=============#
CPPFLAGS=-D_PGTRACK_ -I$(DINC) -L/lib64 -isystem$(ROOTSYS)/include -I$(ROOTSYS)/include
CXXFLAGS := $(CPPFLAGS) -ggdb    

EXTRALIBS = -lMathMore 
EXTRALIBS += `root-config --libs`
ifeq "$(CXX)" "g++"
	EXTRALIBS += -lgfortran
endif

#===============================COMPILING RULES=====================#
all: $(LIBSO) 

$(LIBSO): $(OBJS) $(DOBJ)$(ROOTDICT).o
	@echo Compiling $@
	@echo $(OBJS)
	@if [ -d $(DLIB) ]; then \
		echo $(DLIB); \
	 else \
	    mkdir $(DLIB);\
	 fi;
	$(CXX) -shared -o $@ $^ $(CPPFLAGS) $(EXTRALIBS)
	@cp $(DICTDIR)$(MYLIB)_dict_rdict.pcm $(DLIB)
	@rm -r $(DOBJ)


$(DOBJ)%.o: $(DSRC)%.cxx
	@echo Compiling $@ ...
	@if [ -d $(DOBJ) ]; then \
		echo $(DOBJ); \
	 else \
	    mkdir $(DOBJ);\
	 fi;
#	$(CXX) -o $@ -c -std=c++0x -Wall -fPIC $^ $(CPPFLAGS) $(EXTRALIBS)
	$(CXX) -o $@ -c -Wall -fPIC $^ $(CPPFLAGS) $(EXTRALIBS)


$(DOBJ)$(ROOTDICT).o: $(DICTDIR)$(ROOTDICT).cxx
	@echo Compiling $@...
	$(CXX) -o $@ -c -p -Wall -fPIC $^ $(CPPFLAGS) $(EXTRALIBS)
#	$(CXX) -o $@ -c -std=c++1z -p -Wall -fPIC $^ $(CPPFLAGS) $(EXTRALIBS)

$(DICTDIR)$(ROOTDICT).cxx: $(HEADERS) $(DICTDIR)LinkDef.h
	@echo $@
	rootcint -f $@ $(CPPFLAGS) $^

#==================================Clean=============================#
.PHONY: clean

clean:
	@echo Remove objects and lib.....
	@rm -rf $(DOBJ)* $(ROOTDICT).* $(LIBSO) $(ROOTDICT).* $(DICTDIR)/*.pcm $(DLIB)/*.pcm 

