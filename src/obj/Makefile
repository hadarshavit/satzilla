VPATH = ..

all : ${OBJ}
# DO NOT DELETE

BinSolver.o: BinSolver.h global.h stopwatch.h
BuildSolvers.o: BinSolver.h global.h stopwatch.h
features.o: SATinstance.h BinSolver.o BuildSolvers.o
SATinstance.o: SATinstance.h global.h stopwatch.h
stopwatch.o: stopwatch.h

