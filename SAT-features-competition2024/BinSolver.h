#ifndef __BIN_SOLVER_H__
#define __BIN_SOLVER_H__

class BinSolver
{
public:
  const char* name; //Solver name
  int argc; // number of arguments
  const char** argv;
  int inputFileParam;

  char outFileName[512];
  char solverName[512];
  
  bool outFileCreated;

  BinSolver(const char* _name, int _argc, int _inputFileParam);
  ~BinSolver();

  int spawnBinary(const char* binFile, const char* const argv[], const char* outFileName, int timeout);

  virtual int execute(const char* inputFile, int timeout);
  virtual void cleanup();
};



extern BinSolver* SolverSaps;
extern BinSolver* SolverGsat;
extern BinSolver* SolverZchaff;
extern BinSolver* SolverSatelite;


void BuildSolvers(const char* strseed, const char* outfile);


#endif
