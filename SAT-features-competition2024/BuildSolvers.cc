#include "BinSolver.h"
#include "global.h"

#include "model.h"

#include <stdio.h>
#include <string.h>
#include <initializer_list>

BinSolver *SolverSatelite;
BinSolver *SolverZchaff;
BinSolver *SolverSaps;
BinSolver *SolverGsat;

namespace {
void setSolverArgs(BinSolver *solver, std::initializer_list<const char *> args)
{
  int idx = 1;
  for (const char *arg : args)
    solver->argv[idx++] = arg;
  solver->argv[idx] = nullptr;
}
}

void BuildSolvers(const char *strseed, const char *outfile)
{
  DestroySolvers();

  (void)strseed;

  // --vallst
  SolverSatelite = new BinSolver("sbva", 7, 2);
  setSolverArgs(SolverSatelite, {"-i", nullptr, "-o", outfile, "-t", "30"});

  // -- zchaff07 for compute features
  SolverZchaff = new BinSolver("cadical2023", 3, 1);
  setSolverArgs(SolverZchaff, {nullptr, "--plain"});

  SolverSaps = new BinSolver("ubcsat2006", 18, 2);
  setSolverArgs(SolverSaps, {
                                 "-inst", nullptr, "-alg", "sparrow", "-noimprove", "0.1n", "-r", "stats", outfile,
                                 "best[mean+cv],firstlmstep[mean+median+cv+q10+q90],bestavgimpr[mean+cv],firstlmratio[mean+cv],estacl,numsolve",
                                 "-runs", UBCSAT_NUM_RUNS, "-gtimeout", "5", "-solve", "-v", "sat11"});

  SolverGsat = new BinSolver("ubcsat2006", 16, 2);
  setSolverArgs(SolverGsat, {
                                 "-inst", nullptr, "-alg", "gsat", "-noimprove", "0.5n", "-r", "stats", outfile,
                                 "best[mean+cv],firstlmstep[mean+median+cv+q10+q90],bestavgimpr[mean+cv],firstlmratio[mean+cv],estacl,numsolve",
                                 "-runs", UBCSAT_NUM_RUNS, "-gtimeout", "5", "-solve"});
}

void DestroySolvers()
{
  delete SolverSatelite;
  delete SolverZchaff;
  delete SolverSaps;
  delete SolverGsat;
  SolverSatelite = nullptr;
  SolverZchaff = nullptr;
  SolverSaps = nullptr;
  SolverGsat = nullptr;
}
