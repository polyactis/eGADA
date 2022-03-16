
BoostPythonTarget	= eGADA.so
BoostPythonLib  = -lboost_python38
# "-L ~/lib64/boost/lib" removed from BoostPythonLib. "-lboost_python-mt" changed to "-lboost_python"

SharedLibFlags	= -shared -fPIC
IncludeDirs	= `python3-config --includes` -I`readlink -f ..`
# "-I ~/lib64/boost/include/" removed from IncludeDirs

BoostLib	= -lboost_program_options -lboost_iostreams

CFLAGS	= -O3 
CXXFLAGS	= -fPIC $(CFLAGS) -std=c++0x
CXXCOMPILER	= g++
CXXLDFLAGS	= $(BoostLib)

all:	depend eGADA $(BoostPythonTarget)

# for boost.python target(shared libraries)
$(BoostPythonTarget):	%.so:	eGADA_shared.o BaseGADA_shared.o
	# -l... should be put in the last. If they are put in front of input files, they will be ignored!
	$(CXXCOMPILER) $< BaseGADA_shared.o -DGADALib -o $@ $(SharedLibFlags) -lm $(BoostPythonLib)

eGADA_shared.o:	eGADA.cc
	$(CXXCOMPILER) $< -c -o $@ $(CXXFLAGS) -DGADALib $(IncludeDirs) $(SharedLibFlags)
BaseGADA_shared.o:	BaseGADA.cc
	$(CXXCOMPILER) $< -c -o $@ $(CXXFLAGS) -DGADALib $(IncludeDirs) $(SharedLibFlags)

.cc.o:
	$(CXXCOMPILER) $< -c -o $@ $(CXXFLAGS) $(IncludeDirs) -DGADABIN $(SharedLibFlags)

.c.o:
	$(CXXCOMPILER) $< -c -o $@ $(CXXFLAGS) $(IncludeDirs) -DGADABIN $(SharedLibFlags)


eGADA:	eGADA.o BaseGADA.o BaseGADA.h consoledefines.h
	# `python3-config --libs` and $(BoostPythonLib) would be needed if boost.python is used. However, -DGADABIN renders those code off. 
	$(CXXCOMPILER) eGADA.o BaseGADA.o -o eGADA $(CXXFLAGS) $(IncludeDirs) -DGADABIN -lm $(CXXLDFLAGS)

depend: .depend

.depend: eGADA.cc BaseGADA.cc
	rm -f ./.depend
	-$(CXXCOMPILER) $^ $(CXXFLAGS) $(IncludeDirs) -MM >>./.depend;

include .depend

clean:
	-rm eGADA *.o $(BoostPythonTarget)
