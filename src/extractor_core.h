#ifndef SATZILLA_EXTRACTOR_CORE_H
#define SATZILLA_EXTRACTOR_CORE_H

#include <cstdio>
#include <string>
#include <vector>

struct FeatureOptions {
  bool doBase = false;
  bool doStructure = false;
  bool doNcnfGraphs = false;
  bool doNcnfConstraints = false;
  bool doNcnfRwh = false;
  bool doUnitProbe = false;
  bool doLSProbe = false;
  bool doCl = false;
  bool doDia = false;
  bool doSp = false;
  bool doLobjois = false;
  bool doLP = false;
  int timeoutSeconds = -1;
  int groupTimeoutSeconds = -1;
  int preprocessTimeoutSeconds = -1;
  const char *solverRoot = nullptr;
};

struct FeatureRunResult {
  std::vector<std::string> featureNames;
  std::vector<double> featureValues;
};

void enableAllFeatures(FeatureOptions &opts);
bool hasSelectedFeatureGroup(const FeatureOptions &opts);
std::string resolveExecutableDir(const char *path);

int runFeatureExtraction(const char *inputFile,
                         const FeatureOptions &options,
                         FeatureRunResult &result,
                         std::string &errorMessage);

bool writeFeatureRunResultCsv(FILE *stream,
                              const FeatureRunResult &result,
                              std::string &errorMessage);

bool writeFeatureRunResultCsv(const char *path,
                              const FeatureRunResult &result,
                              std::string &errorMessage);

#endif
