CPP  = g++
CC   = gcc
CXXFLAGS = $(CXXINCS) 

BIN  = libvarsat.a
LINKOBJ=Main.o Solver.o
RM = rm -f	

.PHONY: all clean clean-custom

all: libvarsat.a

clean: clean-custom
	${RM} $(LINKOBJ) $(BIN)

$(BIN): $(LINKOBJ)
	ar r $(BIN) $(LINKOBJ)
	ranlib $(BIN)
	
Main.o: Main.C
	$(CPP) -c Main.C -o Main.o $(CXXFLAGS)
	
Solver.o: Solver.C
	$(CPP) -c Solver.C -o Solver.o $(CXXFLAGS)

