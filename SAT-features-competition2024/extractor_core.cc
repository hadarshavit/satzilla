#include "extractor_core.h"

#include "SATinstance.h"
#include "global.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {
constexpr int kDefaultTimeoutSeconds = 2419200; // 4 weeks
constexpr const char *kSolverSeed = "123456";
constexpr const char *kOutputTemplate = "/outputXXXXXX";
constexpr const char *kTmpConvertedTemplate = "/tmp_wcnf_to_cnf_XXXXXX";
constexpr int kPreprocessorTimeoutSeconds = 35;
constexpr int kSolvedSatCode = 10;
constexpr int kSolvedUnsatCode = 20;
constexpr int kKilledBySignalCode = 137;

struct PreparedInput {
  bool weightedCnf = false;
  const char *selectedInputFile = nullptr;
  std::array<char, 512> convertedTmpFile{};
};

std::string gRuntimeRootStorage;

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

int determineTimeout(const FeatureOptions &options)
{
  if (options.timeoutSeconds >= 0)
    return options.timeoutSeconds;

  char *value = getenv("SATTIMEOUT");
  if (value == nullptr)
    return kDefaultTimeoutSeconds;
  return atoi(value);
}

void resetGlobalState(const FeatureOptions &options)
{
  gSW = Stopwatch();
  gTimeOut = determineTimeout(options);
  preTime = 0.0;
  OrigNumVars = -1;
  OrigNumClauses = -1;
  myTime = 0.0;
  gRuntimeRootStorage = options.solverRoot != nullptr ? options.solverRoot : ".";
  mypath = gRuntimeRootStorage.c_str();
}

bool prepareInputFile(const char *originalInputFile, PreparedInput &prepared, std::string &errorMessage)
{
  prepared.selectedInputFile = originalInputFile;

  std::ifstream infile(originalInputFile);
  if (!infile)
  {
    errorMessage = std::string("Could not read from input file ") + originalInputFile;
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
      errorMessage = "Could not create temporary CNF file";
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
    tmpOutFile << "p cnf " << OrigNumVars << " " << OrigNumClauses << '\n';

    for (std::string line; getline(infile, line);)
    {
      if (line.empty() || line[0] == 'w' || line[0] == 'c' || line.find_first_not_of(' ') == std::string::npos)
        continue;

      const size_t separator = line.find_first_of(" \t");
      if (separator == std::string::npos)
        continue;

      tmpOutFile << line.substr(separator + 1) << '\n';
    }

    prepared.selectedInputFile = prepared.convertedTmpFile.data();
    LogDebug("c WCNF processing time: %f \n", gSW.TotalLap() - preTime);
    break;
  }

  if (OrigNumVars == -1 || OrigNumClauses == -1)
  {
    errorMessage = std::string("Initial read did not find number of vars and clauses in ") + prepared.selectedInputFile;
    return false;
  }

  return true;
}

void cleanupPreparedInput(const PreparedInput &prepared)
{
  if (prepared.weightedCnf && prepared.convertedTmpFile[0] != '\0')
    remove(prepared.convertedTmpFile.data());
}
}

Stopwatch gSW;
int gTimeOut;
double preTime;
int OrigNumVars = -1, OrigNumClauses = -1;
double myTime = 0.0;
const char *mypath = nullptr;

void enableAllFeatures(FeatureOptions &opts)
{
  opts.doBase = true;
  opts.doStructure = true;
  opts.doNcnfGraphs = true;
  opts.doNcnfConstraints = true;
  opts.doNcnfRwh = true;
  opts.doUnitProbe = true;
  opts.doLSProbe = true;
  opts.doLobjois = true;
  opts.doCl = true;
  opts.doDia = true;
  opts.doSp = true;
  opts.doLP = true;
}

bool hasSelectedFeatureGroup(const FeatureOptions &opts)
{
  return opts.doBase || opts.doStructure || opts.doNcnfGraphs || opts.doNcnfConstraints ||
         opts.doNcnfRwh || opts.doUnitProbe || opts.doLSProbe || opts.doCl || opts.doDia ||
         opts.doSp || opts.doLobjois || opts.doLP;
}

std::string resolveExecutableDir(const char *path)
{
  const std::string fullPath(path != nullptr ? path : "");
  const size_t lastSlashPos = fullPath.find_last_of('/');
  if (lastSlashPos == std::string::npos)
    return ".";
  return fullPath.substr(0, lastSlashPos);
}

int runFeatureExtraction(const char *inputFile,
                         const FeatureOptions &rawOptions,
                         FeatureRunResult &result,
                         std::string &errorMessage)
{
  FeatureOptions options = rawOptions;
  if (!hasSelectedFeatureGroup(options))
    options.doBase = true;

  result.featureNames.clear();
  result.featureValues.clear();
  errorMessage.clear();

  PreparedInput prepared;
  std::array<char, 512> solverOutput{};
  SATinstance *sat = nullptr;

  try
  {
    if (!createTempPath(solverOutput.data(), solverOutput.size(), kOutputTemplate))
    {
      errorMessage = "Could not create temporary solver output file";
      return 1;
    }

    resetGlobalState(options);
    BuildSolvers(kSolverSeed, solverOutput.data());
    gSW.Start();

    if (!prepareInputFile(inputFile, prepared, errorMessage))
    {
      DestroySolvers();
      remove(solverOutput.data());
      cleanupPreparedInput(prepared);
      return 1;
    }

    LogDebug("c Orignal number of variables is %d, number of clauses is %d \n", OrigNumVars, OrigNumClauses);
    LogDebug("c run SatELite as pre-processor ... \n");
    LogDebug("c Input file is: %s. Output file is %s\n", prepared.selectedInputFile, solverOutput.data());

    bool doComp = true;
    bool prepSuccessful = false;
    bool solved = false;
    int returnVal = 0;

    if (!prepared.weightedCnf)
    {
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
        sat = new SATinstance(solverOutput.data(), doComp);
        prepSuccessful = true;
      }
    }
    else
    {
      sat = new SATinstance(prepared.selectedInputFile, doComp);
    }

    if (!prepared.weightedCnf)
    {
      preTime = gSW.TotalLap() - myTime;
      sat->start_computation(solved, preTime);
      myTime = gSW.TotalLap();
      LogDebug("c Pre-process time is %f second\n", preTime);
    }

    if (options.doBase)
    {
      preTime = gSW.TotalLap() - myTime;
      returnVal = sat->computeFeatures(doComp, !prepared.weightedCnf);
      LogDebug("c Base features time is %f second\n", gSW.TotalLap() - preTime);
      if (sat->getNumVals() == 0 || sat->getNumClaus() == 0)
      {
        doComp = false;
        LogDebug("c Instance can be solved by unit propodation alone!!!\n");
      }
    }

    bool timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doStructure && !timeout)
      sat->structureFeatures(doComp);

    timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doNcnfGraphs && !timeout)
      sat->newCnfGraphFeatures(doComp);

    timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doNcnfConstraints && !timeout)
      sat->newCnfConstraintFeatures(doComp);

    timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doNcnfRwh && !timeout)
      sat->newCnfRwhFeatures(doComp);

    timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doDia && !timeout && returnVal != VCG_TIMEOUT_CODE)
      sat->init_diameter(doComp);

    timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doCl && !timeout)
      sat->cl_prob(prepSuccessful ? solverOutput.data() : prepared.selectedInputFile, doComp);

    timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doSp && !timeout)
      sat->sp(doComp);

    timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doUnitProbe && !timeout)
      sat->unitPropProbe(false, doComp);

    timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doLP && !timeout)
      sat->compute_lp(doComp);

    timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doLSProbe && !timeout)
    {
      const char *lsInput = prepSuccessful ? solverOutput.data() : prepared.selectedInputFile;
      sat->localSearchProbeSaps(lsInput, doComp);
      sat->localSearchProbeGsat(lsInput, doComp);
    }

    timeout = gSW.TotalLap() > TOTAL_TIMEOUT;
    if (options.doLobjois && !timeout)
      sat->lobjoisProbe(false, doComp);

    sat->finish_computation();
    for (int i = 0; i < sat->getNumFeatures(); ++i)
    {
      result.featureNames.push_back(sat->getFeatureNames()[i]);
      result.featureValues.push_back(sat->getFeatureVals()[i]);
    }
  }
  catch (const std::exception &exc)
  {
    errorMessage = exc.what();
  }
  catch (...)
  {
    errorMessage = "Unknown extraction failure";
  }

  delete sat;
  DestroySolvers();
  if (solverOutput[0] != '\0')
    remove(solverOutput.data());
  cleanupPreparedInput(prepared);

  return errorMessage.empty() ? 0 : 1;
}

bool writeFeatureRunResultCsv(FILE *stream,
                              const FeatureRunResult &result,
                              std::string &errorMessage)
{
  if (stream == nullptr)
  {
    errorMessage = "Output stream is null";
    return false;
  }
  if (result.featureNames.size() != result.featureValues.size())
  {
    errorMessage = "Feature name/value size mismatch";
    return false;
  }

  for (size_t i = 0; i < result.featureNames.size(); ++i)
  {
    fprintf(stream, "%s%s", result.featureNames[i].c_str(),
            (i + 1 == result.featureNames.size()) ? "\n" : ",");
  }
  for (size_t i = 0; i < result.featureValues.size(); ++i)
  {
    fprintf(stream, "%.9lf%s", result.featureValues[i],
            (i + 1 == result.featureValues.size()) ? "\n" : ",");
  }
  return true;
}

bool writeFeatureRunResultCsv(const char *path,
                              const FeatureRunResult &result,
                              std::string &errorMessage)
{
  FILE *stream = fopen(path, "w");
  if (stream == nullptr)
  {
    errorMessage = std::string("Could not open output file ") + path;
    return false;
  }
  const bool ok = writeFeatureRunResultCsv(stream, result, errorMessage);
  fclose(stream);
  return ok;
}
