#ifndef SATZILLA_C_API_H
#define SATZILLA_C_API_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct satzilla_options_t {
  int do_base;
  int do_structure;
  int do_ncnf_graphs;
  int do_ncnf_constraints;
  int do_ncnf_rwh;
  int do_unit_probe;
  int do_ls_probe;
  int do_cl;
  int do_dia;
  int do_sp;
  int do_lobjois;
  int do_lp;
  int timeout_seconds;
  const char *solver_root;
} satzilla_options_t;

typedef struct satzilla_result_t {
  int status_code;
  size_t feature_count;
  const char **feature_names;
  double *feature_values;
  char *error_message;
} satzilla_result_t;

int satzilla_extract_from_dimacs_path(const char *input_path,
                                      const satzilla_options_t *options,
                                      satzilla_result_t *result);

int satzilla_extract_from_clauses(int num_vars,
                                  size_t num_clauses,
                                  const int *literals,
                                  const size_t *clause_offsets,
                                  const satzilla_options_t *options,
                                  satzilla_result_t *result);

void satzilla_free_result(satzilla_result_t *result);

#ifdef __cplusplus
}
#endif

#endif
