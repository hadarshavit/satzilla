#include "SATinstance.h"
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <array>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include "global.h"
using namespace std;

namespace {
constexpr int kDefaultTimeoutSeconds = 2419200; // 4 weeks
constexpr const char *kSolverSeed = "123456";
constexpr const char *kOutputTemplate = "/outputXXXXXX";
constexpr const char *kTmpConvertedTemplate = "/tmp_wcnf_to_cnf_XXXXXX";
constexpr int kPreprocessorTimeoutSeconds = 35;
constexpr int kSolvedSatCode = 10;
constexpr int kSolvedUnsatCode = 20;
constexpr int kKilledBySignalCode = 137;

bool createTempPath(char *buffer, size_t bufferSize, const char *suffixTemplate)
{
  if (snprintf(buffer, bufferSize, "%s%s", P_tmpdir, suffixTemplate) >= static_cast<int>(bufferSize))
    return false;
  int fd = mkstemp(buffer);
  if (fd == -1)
    return false;
  close(fd);
  return true;
}

struct FeatureOptions {
  bool doBase = false;
  bool doUnitProbe = false;
  bool doLSProbe = false;
  bool doCl = false;
  bool doDia = false;
  bool doSp = false;
  bool doLobjois = false;
  bool doLP = false;
  const char *inputFile = nullptr;
  const char *outputFile = nullptr;
};

enum class ParseResult {
  Ok,
  Help,
  Error,
};

struct PreparedInput {
  bool weightedCnf = false;
  const char *selectedInputFile = nullptr;
  std::array<char, 512> convertedTmpFile{};
};

void enableAllFeatures(FeatureOptions &opts)
{
  opts.doBase = true;
  opts.doUnitProbe = true;
  opts.doLSProbe = true;
  opts.doLobjois = true;
  opts.doCl = true;
  opts.doDia = true;
  opts.doSp = true;
  opts.doLP = true;
}

void printUsage(const char *progName)
{
  cerr << "Usage: " << progName
       << " [options] infile [outfile]\n"
       << "Options:\n"
       << "  -all       Enable all feature groups\n"
       << "  -base      Base features\n"
       << "  -sp        Survey propagation\n"
       << "  -dia       Diameter\n"
       << "  -cl        Clause learning\n"
       << "  -lp        Linear programming\n"
       << "  -unit      Unit propagation probe\n"
       << "  -ls        Local search probes (GSAT + SAPS)\n"
       << "  -lobjois   Lobjois probe\n"
       << "  -h, --help Show this help\n";
}

ParseResult parseArgs(int argc, char **argv, FeatureOptions &opts)
{
  bool sawFeatureFlag = false;
  int index = 1;
  while (index < argc && argv[index][0] == '-')
  {
    string flag = argv[index];
    string lowerFlag = flag;
    transform(lowerFlag.begin(), lowerFlag.end(), lowerFlag.begin(), [](unsigned char c) {
      return static_cast<char>(tolower(c));
    });

    if (lowerFlag == "-h" || lowerFlag == "--help")
    {
      printUsage(argv[0]);
      return ParseResult::Help;
    }
    if (lowerFlag == "-all")
    {
      enableAllFeatures(opts);
      sawFeatureFlag = true;
    }
    else if (lowerFlag == "-base")
    {
      opts.doBase = true;
      sawFeatureFlag = true;
    }
    else if (lowerFlag == "-unit")
    {
      opts.doUnitProbe = true;
      sawFeatureFlag = true;
    }
    else if (lowerFlag == "-lp")
    {
      opts.doLP = true;
      sawFeatureFlag = true;
    }
    else if (lowerFlag == "-sp")
    {
      opts.doSp = true;
      sawFeatureFlag = true;
    }
    else if (lowerFlag == "-dia")
    {
      opts.doDia = true;
      sawFeatureFlag = true;
    }
    else if (lowerFlag == "-cl")
    {
      opts.doCl = true;
      sawFeatureFlag = true;
    }
    else if (lowerFlag == "-ls")
    {
      opts.doLSProbe = true;
      sawFeatureFlag = true;
    }
    else if (lowerFlag == "-lobjois")
    {
      opts.doLobjois = true;
      sawFeatureFlag = true;
    }
    else
    {
      cerr << "Unknown option: " << flag << "\n";
      printUsage(argv[0]);
      return ParseResult::Error;
    }
    index++;
  }

  if (index >= argc)
  {
    cerr << "Missing input file.\n";
    printUsage(argv[0]);
    return ParseResult::Error;
  }

  opts.inputFile = argv[index++];

  if (index < argc)
    opts.outputFile = argv[index++];

  if (index < argc)
  {
    cerr << "Too many arguments. Expected at most one output file.\n";
    printUsage(argv[0]);
    return ParseResult::Error;
  }

  if (!sawFeatureFlag)
    opts.doBase = true;

  return ParseResult::Ok;
}

std::string resolveExecutableDir(const char *argv0)
{
  const string fullPath(argv0);
  const size_t lastSlashPos = fullPath.find_last_of('/');
  if (lastSlashPos == string::npos)
    return ".";
  return fullPath.substr(0, lastSlashPos);
}

bool prepareInputFile(const char *originalInputFile, PreparedInput &prepared)
{
  prepared.selectedInputFile = originalInputFile;

  ifstream infile(originalInputFile);
  if (!infile)
  {
    fprintf(stderr, "c Error: Could not read from input file %s.\n", originalInputFile);
    return false;
  }

  char strbuf[1024];
  char wordbuf[256];

  while (infile >> wordbuf)
  {
    if (wordbuf[0] == 'p')
    {
      infile >> strbuf;
      infile >> OrigNumVars >> OrigNumClauses;
      break;
    }

    if (!strstr(wordbuf, "c{"))
      continue;

    prepared.weightedCnf = true;
    if (!createTempPath(prepared.convertedTmpFile.data(), prepared.convertedTmpFile.size(), kTmpConvertedTemplate))
    {
      fprintf(stderr, "c Error: Could not create temporary CNF file.\n");
      return false;
    }

    LogDebug("c Processing WCNF file to CNF... \n");
    preTime = gSW.TotalLap() - myTime;

    while (infile >> wordbuf)
    {
      if (strstr(wordbuf, "c}"))
        break;
      if (OrigNumVars == -1 && strstr(wordbuf, "nvars\":"))
        infile >> OrigNumVars;
      else if (OrigNumClauses == -1 && strstr(wordbuf, "ncls\":"))
        infile >> OrigNumClauses;
    }

    std::ofstream tmpOutFile(prepared.convertedTmpFile.data());
    tmpOutFile << "p cnf " << OrigNumVars << " " << OrigNumClauses << endl;

    for (std::string line; getline(infile, line);)
    {
      if (line.empty() || line[0] == 'w' || line[0] == 'c' || line.find_first_not_of(' ') == std::string::npos)
        continue;

      const size_t separator = line.find_first_of(" \t");
      if (separator == std::string::npos)
        continue;

      line = line.substr(separator + 1);
      tmpOutFile << line << endl;
    }

    prepared.selectedInputFile = prepared.convertedTmpFile.data();
    LogDebug("c WCNF processing time: %f \n", gSW.TotalLap() - preTime);
    break;
  }

  if (OrigNumVars == -1 || OrigNumClauses == -1)
  {
    fprintf(stderr, "\nc ERROR: Initial read did not find number of vars and clauses, reached Premature EOF in %s ", prepared.selectedInputFile);
    return false;
  }

  return true;
}
}

Stopwatch gSW;
int gTimeOut;
double preTime;
int OrigNumVars = -1, OrigNumClauses = -1;
double myTime = 0.0;
const char *mypath;

// int numFeats=64;
int getTimeOut(void)
{
  char *value = getenv("SATTIMEOUT");
  if (value == NULL)
    return kDefaultTimeoutSeconds;

  return atoi(value);
}

int main(int argc, char **argv)
{
  FeatureOptions options;
  const ParseResult parseResult = parseArgs(argc, argv, options);
  if (parseResult == ParseResult::Help)
    return 0;
  if (parseResult != ParseResult::Ok)
    return 1;

  bool doComp = true;
  const string myownpath = resolveExecutableDir(argv[0]);
  mypath = myownpath.c_str();

  const bool doBase = options.doBase;
  const bool doUnitProbe = options.doUnitProbe;
  const bool doLSProbe = options.doLSProbe;
  const bool doCl = options.doCl;
  const bool doDia = options.doDia;
  const bool doSp = options.doSp;
  const bool doLobjois = options.doLobjois;
  const bool doLP = options.doLP;

  bool letsgo = true;
  bool prep_successful = false;

  std::array<char, 512> outfile;
  if (!createTempPath(outfile.data(), outfile.size(), kOutputTemplate))
  {
    fprintf(stderr, "c Error: Could not create temp output file.\n");
    return 1;
  }

  gTimeOut = getTimeOut();
  BuildSolvers(kSolverSeed, outfile.data());
  gSW.Start();

  PreparedInput prepared;
  if (!prepareInputFile(options.inputFile, prepared))
    return 1;

  LogDebug("c Orignal number of variables is %d, number of clauses is %d \n", OrigNumVars, OrigNumClauses);
  bool solved = false;
  LogDebug("c run SatELite as pre-processor ... \n");
  int returnVal = 0;
  LogDebug("c Input file is: %s. Output file is %s\n", prepared.selectedInputFile, outfile.data());

  // We can perhaps disable the pre-processor for WCNF to speed up feature computation? -> Looks currently to only be a small slowdown
  SATinstance *sat;
  if (!prepared.weightedCnf){
    preTime = gSW.TotalLap() - myTime;
    returnVal = SolverSatelite->execute(prepared.selectedInputFile, kPreprocessorTimeoutSeconds);
    myTime = gSW.TotalLap();
    LogDebug("c SatELite pre-process time is %f second\n", myTime - preTime);
    
    if (returnVal == kSolvedSatCode || returnVal == kSolvedUnsatCode)
    {
        LogDebug("c This instance is solved by pre-processor with %d!\n", returnVal);
        solved = true;
        doComp = false;
    }

    SolverSatelite->cleanup();
    if (returnVal == kKilledBySignalCode)
    {
        sat = new SATinstance(prepared.selectedInputFile, doComp);
    }
    else
    {
        sat = new SATinstance(outfile.data(), doComp);
        prep_successful = true;
    }
  }
  else {
    sat = new SATinstance(prepared.selectedInputFile, doComp);
    //prep_sucessfull = true;
  }

  if (!prepared.weightedCnf){  // Pre processing mainly makes sense for regular SAT, although the potential formula simplification is usefull.
    preTime = gSW.TotalLap() - myTime;
    sat->start_computation(solved, preTime);

    myTime = gSW.TotalLap();
    LogDebug("c Pre-process time is %f second\n", preTime);
  }

  if (doBase && letsgo)
  {
    preTime = gSW.TotalLap() - myTime;
    returnVal = sat->computeFeatures(doComp, !prepared.weightedCnf);
    double feature_time = gSW.TotalLap() - preTime;
    LogDebug("c Base features time is %f second\n", feature_time);
    if (sat->getNumVals() == 0 || sat->getNumClaus() == 0)
    {
      doComp = false;
      LogDebug("c Instance can be solved by unit propodation alone!!!\n");
    }
  }

  bool timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
  if (doDia && letsgo && !timeout && returnVal != VCG_TIMEOUT_CODE)
    sat->init_diameter(doComp);

  timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
  if (doCl && letsgo && !timeout)
  {
    if (prep_successful)
        sat->cl_prob(outfile.data(), doComp);
    else
      sat->cl_prob(prepared.selectedInputFile, doComp);

  }

  timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
  if (doSp && letsgo && !timeout)
    sat->sp(doComp);

  timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
  if (doUnitProbe && letsgo && !timeout)
    sat->unitPropProbe(false, doComp);

  timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
  if (doLP && letsgo && !timeout)
    sat->compute_lp(doComp);

  timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
  if (doLSProbe && letsgo && !timeout)
  {
    if (prep_successful){
    sat->localSearchProbeSaps(outfile.data(), doComp);
    sat->localSearchProbeGsat(outfile.data(), doComp);
    }
    else{
      sat->localSearchProbeSaps(prepared.selectedInputFile, doComp);
      sat->localSearchProbeGsat(prepared.selectedInputFile, doComp);
    }
  }

  timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
  if (doLobjois && letsgo && !timeout)
    sat->lobjoisProbe(false, doComp);

  sat->finish_computation();
  if (letsgo)
    if (options.outputFile != nullptr)
    {
      sat->writeFeatNamesToFile(options.outputFile);
      sat->writeFeaturesToFile(options.outputFile);
    }
    else
    {
      sat->writeFeatNamesToFile(stdout);
      sat->writeFeaturesToFile(stdout);
    }

  delete sat;
  remove(outfile.data());
  if (prepared.weightedCnf){
    remove(prepared.convertedTmpFile.data());
  }

  return 0;
}
