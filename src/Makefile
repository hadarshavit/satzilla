###################################################################
#
# You should be able to just fill in the vars here for any project.
# You will need dirs "obj" and "debugobj" each containing
# Makefiles containing only the following:
#
# VPATH = ..
#
# all : ${OBJ}
#
###################################################################

COMMONOBJ = extractor_core.o SATinstance.o stopwatch.o BinSolver.o BuildSolvers.o

FEATOBJ = features.o $(COMMONOBJ)
SHAREDNAME = libsatzilla_features.so

EXECNAME = features

INCLUDE =

LIB = -l:liblpk5.a -lvarsat -ldl  -lz

LIBDIRS = -Llp_solve_5.0 -LVARSAT 
SHAREDLIBDIRS = -Llp_solve_5.0

RELEASEFLAGS = -O3 -DDEBUG -g -std=c++17
DEBUGFLAGS = -O0 -DDEBUG -g3 -std=c++17 -Wall -Wextra -Wpedantic

SOLVERLIB =  

release:
	cd obj; $(MAKE) all "OBJ = ${FEATOBJ}" "CPPFLAGS = -DFEATURES ${RELEASEFLAGS} ${INCLUDE}";
	$(MAKE) allfeat "VPATH = obj" "OBJ = $(FEATOBJ)" "FLAGS = -o features ${RELEASEFLAGS}"

libubcsat.a:
	cd UBCSAT; $(MAKE);
	
libvarsat.a:
	cd VARSAT; $(MAKE);

liblpk.a:
	cd lp_solve_5.0/lpsolve5; ./ccc; cp liblpk5.a ../liblpk5.a

liblpk.so:
	cd lp_solve_5.0/lpsolve5; ./ccc; cp liblpk5.so ../liblpk5.so; cp liblpk5.so ../../liblpk5.so


allfeat : $(OBJ)
	+$(MAKE) libvarsat.a; $(MAKE) liblpk.a; g++  ${LIBDIRS} $^ ${FLAGS} $(LIB)

feature:
	cd obj; $(MAKE) all "OBJ = ${FEATOBJ}" "CPPFLAGS = -DFEATURES ${RELEASEFLAGS} ${INCLUDE}";
	$(MAKE) allfeat "VPATH = obj" "OBJ = $(FEATOBJ)" "FLAGS = -o features ${RELEASEFLAGS}"

shared:
	$(MAKE) liblpk.so
	g++ -shared -fPIC ${RELEASEFLAGS} -Wl,-rpath,'$$ORIGIN' -o $(SHAREDNAME) \
		extractor_core.cc satzilla_c_api.cc SATinstance.cc stopwatch.cpp BinSolver.cc BuildSolvers.cc \
		VARSAT/Main.C VARSAT/Solver.C $(SHAREDLIBDIRS) -llpk5 -ldl -lz -pthread


.PHONY: clean depend again

test:
	$(MAKE) feature
	$(MAKE) shared
	python3 tests/test_package.py
	python3 tests/test_python_interface.py
	python3 tests/test_package_install.py

clean-all:
	rm -f obj/*.o features $(SHAREDNAME);  cd VARSAT; $(MAKE) clean ; cd ../lp_solve_5.0;  $(MAKE) clean

clean:
	rm -f obj/*.o features $(SHAREDNAME) lp_solve_5.0/liblpk5.so liblpk5.so

depend:
	makedepend -Y -fobj/Makefile *.cc *.cpp >& /dev/null
	makedepend -Y -fdebugobj/Makefile *.cc *.cpp >& /dev/null

again:
	$(MAKE) clean; $(MAKE) depend; $(MAKE)
