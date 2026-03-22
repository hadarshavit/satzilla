#include "extractor_core.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <iostream>
#include <string>

using namespace std;

namespace {
enum class ParseResult {
  Ok,
  Help,
  Error,
};

void printUsage(const char *progName)
{
  cerr << "Usage: " << progName
       << " [options] infile [outfile]\n"
       << "Options:\n"
       << "  -all       Enable all feature groups\n"
       << "  -base      Base features\n"
       << "  -structure Structure features (Ansotegui)\n"
       << "  -ncnf-graphs      New-CNF graph features\n"
       << "  -ncnf-constraints New-CNF constraint features\n"
       << "  -ncnf-rwh         New-CNF recursive weight heuristic features\n"
       << "  -sp        Survey propagation\n"
       << "  -dia       Diameter\n"
       << "  -cl        Clause learning\n"
       << "  -lp        Linear programming\n"
       << "  -unit      Unit propagation probe\n"
       << "  -ls        Local search probes (GSAT + SAPS)\n"
       << "  -lobjois   Lobjois probe\n"
       << "  -h, --help Show this help\n";
}

ParseResult parseArgs(int argc,
                      char **argv,
                      FeatureOptions &opts,
                      const char *&inputFile,
                      const char *&outputFile)
{
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
      enableAllFeatures(opts);
    else if (lowerFlag == "-base")
      opts.doBase = true;
    else if (lowerFlag == "-structure")
      opts.doStructure = true;
    else if (lowerFlag == "-ncnf-graphs")
      opts.doNcnfGraphs = true;
    else if (lowerFlag == "-ncnf-constraints")
      opts.doNcnfConstraints = true;
    else if (lowerFlag == "-ncnf-rwh")
      opts.doNcnfRwh = true;
    else if (lowerFlag == "-unit")
      opts.doUnitProbe = true;
    else if (lowerFlag == "-lp")
      opts.doLP = true;
    else if (lowerFlag == "-sp")
      opts.doSp = true;
    else if (lowerFlag == "-dia")
      opts.doDia = true;
    else if (lowerFlag == "-cl")
      opts.doCl = true;
    else if (lowerFlag == "-ls")
      opts.doLSProbe = true;
    else if (lowerFlag == "-lobjois")
      opts.doLobjois = true;
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

  inputFile = argv[index++];
  outputFile = (index < argc) ? argv[index++] : nullptr;

  if (index < argc)
  {
    cerr << "Too many arguments. Expected at most one output file.\n";
    printUsage(argv[0]);
    return ParseResult::Error;
  }

  return ParseResult::Ok;
}
}

int main(int argc, char **argv)
{
  FeatureOptions options;
  const char *inputFile = nullptr;
  const char *outputFile = nullptr;
  const ParseResult parseResult = parseArgs(argc, argv, options, inputFile, outputFile);
  if (parseResult == ParseResult::Help)
    return 0;
  if (parseResult != ParseResult::Ok)
    return 1;

  const string runtimeRoot = resolveExecutableDir(argv[0]);
  options.solverRoot = runtimeRoot.c_str();

  FeatureRunResult result;
  string errorMessage;
  const int status = runFeatureExtraction(inputFile, options, result, errorMessage);
  if (status != 0)
  {
    cerr << "c Error: " << errorMessage << "\n";
    return status;
  }

  bool wrote = false;
  if (outputFile != nullptr)
    wrote = writeFeatureRunResultCsv(outputFile, result, errorMessage);
  else
    wrote = writeFeatureRunResultCsv(stdout, result, errorMessage);

  if (!wrote)
  {
    cerr << "c Error: " << errorMessage << "\n";
    return 1;
  }

  return 0;
}
