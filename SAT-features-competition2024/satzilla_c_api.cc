#include "satzilla_c_api.h"

#include "extractor_core.h"

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <mutex>
#include <string>
#include <unistd.h>
#include <vector>
#include <fcntl.h>

namespace {
std::mutex gExtractionMutex;

class StdoutSilencer
{
public:
  StdoutSilencer()
      : savedFd_(-1), nullFd_(-1)
  {
    fflush(stdout);
    savedFd_ = dup(STDOUT_FILENO);
    nullFd_ = open("/dev/null", O_WRONLY);
    if (savedFd_ != -1 && nullFd_ != -1)
      dup2(nullFd_, STDOUT_FILENO);
  }

  ~StdoutSilencer()
  {
    fflush(stdout);
    if (savedFd_ != -1)
    {
      dup2(savedFd_, STDOUT_FILENO);
      close(savedFd_);
    }
    if (nullFd_ != -1)
      close(nullFd_);
  }

private:
  int savedFd_;
  int nullFd_;
};

FeatureOptions toFeatureOptions(const satzilla_options_t *options)
{
  FeatureOptions converted;
  if (options == nullptr)
    return converted;

  converted.doBase = options->do_base != 0;
  converted.doStructure = options->do_structure != 0;
  converted.doNcnfGraphs = options->do_ncnf_graphs != 0;
  converted.doNcnfConstraints = options->do_ncnf_constraints != 0;
  converted.doNcnfRwh = options->do_ncnf_rwh != 0;
  converted.doUnitProbe = options->do_unit_probe != 0;
  converted.doLSProbe = options->do_ls_probe != 0;
  converted.doCl = options->do_cl != 0;
  converted.doDia = options->do_dia != 0;
  converted.doSp = options->do_sp != 0;
  converted.doLobjois = options->do_lobjois != 0;
  converted.doLP = options->do_lp != 0;
  converted.timeoutSeconds = options->timeout_seconds;
  converted.solverRoot = options->solver_root;
  return converted;
}

void resetResult(satzilla_result_t *result)
{
  if (result == nullptr)
    return;
  result->status_code = 0;
  result->feature_count = 0;
  result->feature_names = nullptr;
  result->feature_values = nullptr;
  result->error_message = nullptr;
}

void assignError(satzilla_result_t *result, int statusCode, const std::string &message)
{
  resetResult(result);
  if (result == nullptr)
    return;
  result->status_code = statusCode;
  result->error_message = strdup(message.c_str());
}

int copyFeatureResult(const FeatureRunResult &input, satzilla_result_t *result)
{
  resetResult(result);
  result->feature_count = input.featureNames.size();
  result->status_code = 0;

  if (input.featureNames.empty())
    return 0;

  result->feature_names = static_cast<const char **>(calloc(input.featureNames.size(), sizeof(char *)));
  result->feature_values = static_cast<double *>(calloc(input.featureValues.size(), sizeof(double)));
  if (result->feature_names == nullptr || result->feature_values == nullptr)
  {
    assignError(result, 2, "Out of memory while copying feature result");
    return 2;
  }

  for (size_t i = 0; i < input.featureNames.size(); ++i)
  {
    result->feature_names[i] = strdup(input.featureNames[i].c_str());
    if (result->feature_names[i] == nullptr)
    {
      assignError(result, 2, "Out of memory while copying feature names");
      return 2;
    }
    result->feature_values[i] = input.featureValues[i];
  }

  return 0;
}
}

extern "C" int satzilla_extract_from_dimacs_path(const char *input_path,
                                                  const satzilla_options_t *options,
                                                  satzilla_result_t *result)
{
  std::lock_guard<std::mutex> lock(gExtractionMutex);
  if (result == nullptr)
    return 2;
  resetResult(result);

  if (input_path == nullptr)
  {
    assignError(result, 2, "Input path is null");
    return result->status_code;
  }

  FeatureRunResult featureResult;
  std::string errorMessage;
  StdoutSilencer silencer;
  const int status = runFeatureExtraction(input_path, toFeatureOptions(options), featureResult, errorMessage);
  if (status != 0)
  {
    assignError(result, status, errorMessage);
    return status;
  }

  return copyFeatureResult(featureResult, result);
}

extern "C" int satzilla_extract_from_clauses(int num_vars,
                                              size_t num_clauses,
                                              const int *literals,
                                              const size_t *clause_offsets,
                                              const satzilla_options_t *options,
                                              satzilla_result_t *result)
{
  std::lock_guard<std::mutex> lock(gExtractionMutex);
  if (result == nullptr)
    return 2;
  resetResult(result);

  if (num_vars < 0)
  {
    assignError(result, 2, "Number of variables must be non-negative");
    return result->status_code;
  }
  if (num_clauses > 0 && (literals == nullptr || clause_offsets == nullptr))
  {
    assignError(result, 2, "Clause buffers are null");
    return result->status_code;
  }

  char tempPath[] = "/tmp/satzilla_clauses_XXXXXX.cnf";
  int fd = mkstemps(tempPath, 4);
  if (fd == -1)
  {
    assignError(result, 2, "Could not create temporary CNF file");
    return result->status_code;
  }

  FILE *stream = fdopen(fd, "w");
  if (stream == nullptr)
  {
    close(fd);
    unlink(tempPath);
    assignError(result, 2, "Could not open temporary CNF file");
    return result->status_code;
  }

  fprintf(stream, "p cnf %d %zu\n", num_vars, num_clauses);
  for (size_t clauseIdx = 0; clauseIdx < num_clauses; ++clauseIdx)
  {
    const size_t start = clause_offsets[clauseIdx];
    const size_t end = clause_offsets[clauseIdx + 1];
    if (end < start)
    {
      fclose(stream);
      unlink(tempPath);
      assignError(result, 2, "Clause offsets are not monotonically increasing");
      return result->status_code;
    }

    for (size_t litIdx = start; litIdx < end; ++litIdx)
      fprintf(stream, "%d ", literals[litIdx]);
    fprintf(stream, "0\n");
  }
  fclose(stream);

  FeatureRunResult featureResult;
  std::string errorMessage;
  StdoutSilencer silencer;
  const int status = runFeatureExtraction(tempPath, toFeatureOptions(options), featureResult, errorMessage);
  unlink(tempPath);
  if (status != 0)
  {
    assignError(result, status, errorMessage);
    return status;
  }

  return copyFeatureResult(featureResult, result);
}

extern "C" void satzilla_free_result(satzilla_result_t *result)
{
  if (result == nullptr)
    return;

  if (result->feature_names != nullptr)
  {
    for (size_t i = 0; i < result->feature_count; ++i)
      free(const_cast<char *>(result->feature_names[i]));
    free(const_cast<char **>(result->feature_names));
  }

  free(result->feature_values);
  free(result->error_message);
  resetResult(result);
}
