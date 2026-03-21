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

COMMONOBJ = SATinstance.o  stopwatch.o  BinSolver.o BuildSolvers.o

FEATOBJ = features.o $(COMMONOBJ)

EXECNAME = features

INCLUDE =

LIB = -l:liblpk5.a -lvarsat -ldl  -lz

LIBDIRS = -Llp_solve_5.0 -LVARSAT 

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


allfeat : $(OBJ)
	+$(MAKE) libvarsat.a; $(MAKE) liblpk.a; g++  ${LIBDIRS} $^ ${FLAGS} $(LIB)

feature:
	cd obj; $(MAKE) all "OBJ = ${FEATOBJ}" "CPPFLAGS = -DFEATURES ${RELEASEFLAGS} ${INCLUDE}";
	$(MAKE) allfeat "VPATH = obj" "OBJ = $(FEATOBJ)" "FLAGS = -o features ${RELEASEFLAGS}"


.PHONY: clean depend again

clean-all:
	rm -f obj/*.o features;  cd VARSAT; $(MAKE) clean ; cd ../lp_solve_5.0;  $(MAKE) clean

clean:
	rm -f obj/*.o features

depend:
	makedepend -Y -fobj/Makefile *.cc *.cpp >& /dev/null
	makedepend -Y -fdebugobj/Makefile *.cc *.cpp >& /dev/null

again:
	$(MAKE) clean; $(MAKE) depend; $(MAKE)
