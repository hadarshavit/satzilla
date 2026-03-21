#ifndef _SAT_INSTANCE_H
#define _SAT_INSTANCE_H

#include "BinSolver.h"
#include "global.h"

#include <stdio.h>

#include <map>
#include <set>
#include <string>
#include <stack>
#include <vector>
using namespace std;

class SATinstance
{

private:
  static const char *badFeatNames[];
  static const int numBadFeats;

  typedef enum
  {
    TRUE_VAL,
    FALSE_VAL,
    UNASSIGNED,
    IRRELEVANT
  } varState;
  typedef enum
  {
    ACTIVE,
    PASSIVE
  } clauseState;

  /* These never change after constructor. */
  int numVars, numClauses;
  int **clauses;

  vector<int> *negClausesWithVar;
  vector<int> *posClausesWithVar;

  /* These change as tree is explored. */
  varState *varStates;
  clauseState *clauseStates;
  int numActiveVars, numActiveClauses;
  int *numActiveClausesWithVar;
  int *numBinClausesWithVar;
  int *clauseLengths;

  stack<int> unitClauses;

  stack<int> reducedClauses;
  stack<int> numReducedClauses;

  stack<int> reducedVars;
  stack<int> numReducedVars;

  inline vector<int> &clausesWithLit(int lit);
  bool reduceClauses(int lit, int &numClausesReduced, int &numVarsReduced);
  bool setVarAndProp(int var, bool val);
  bool unitprop(int &numClausesReduced, int &numVarsReduced);
  // backtrack undoes one call of setVar *or* unitprop
  void backtrack();

  void retardedSearch();
  void outputAssignment(); // for testing.

  /* OUR ACCESSORS */

  bool *active_clause;

  int currentClause;
  inline int *firstClause();
  int *nextClause();

  // two lit iterators are provided so that two
  // clauses can be iterated over simultaneously

  int currentLit, currentLit2;
  int *currentClauseForLitIter, *currentClauseForLitIter2;
  inline int firstLitInClause(int *clause);
  inline int nextLitInClause();
  inline int firstLitInClause2(int *clause);
  inline int nextLitInClause2();

  int currentClauseWithVar;
  int currentVarForClauseIter;
  bool posClauses;
  inline int *firstClauseWithVar(int var, bool positive);
  int *nextClauseWithVar();

  void print();
  void setVar(int var);
  int varsInMostBinaryClauses(int *array);
  double getAvePropsFromVars(int *toSet, bool *vals, int &cumulProps, int depth);

  int computer_diameter(const int varnr, int nrofliterals);
  int *dia_clause_table, *dstack, *var_depth, *clause_stamps;
  int CCS;
  int **dia_clause_list;

  /* feature helpers */

  bool isVarInClause(int *clause, int var);
  inline double square(double num) { return num * num; }

  inline double mean(int *array, int num);
  inline double mean(double *array, int num);
  inline double stdev(int *array, int num, double mean);
  inline double stdev(double *array, int num, double mean);
  inline int array_min(int *array, int num);
  inline double array_min(double *array, int num);
  inline int array_max(int *array, int num);
  inline double array_max(double *array, int num);
  inline double array_entropy(double *array, int num, int vals, int maxval);
  inline double array_entropy(int *array, int num, int vals);

  // this is true only for active vars that occur
  // in some active clause
  bool isActiveVariable(int i);

  // write out node stats
  // could these stats all be computed in one pass?
  void writeStats(int *array, int num, const char *name);
  void writeStatsQ(double *array, int num, const char *name);
  void writeStats(double *array, int num, const char *name);
  void writeStatsSTDEV(int *array, int num, const char *name);
  void writeStatsSTDEV(double *array, int num, const char *name);

  char *featureNames[MAX_FEATURES]; // list of all feature names
  int indexCount;                   // count of how many features we currently have (built up using writeFeature() )
  double featureVals[MAX_FEATURES]; // the array of actual feature  values
  //   char featurePrefix[128];

  /* these three for the clause graph */
  bool cgTimeout;
  void neighborClauses(int *clause, set<int *> *neighbor);
  bool conflicted(int *clause1, int *clause2);

  int *test_flag;
  int lp_return_val;

  void mkVarTranslation(map<int, int> *trans_for, map<int, int> *trans_back);

  const char *inputFileName;

public:
  SATinstance(const char *filename, bool doComp, long seed = 0);
  ~SATinstance();

  int computeFeatures(bool doComp, bool doClauseGraphFeatures);
  void clauseGraphFeatures(bool realCC);
  int sp(bool doComp);
  int compute_lp(bool doComp);
  int cl_prob(const char *outfile, bool doComp);
  int init_diameter(bool doComp);
  int unitPropProbe(bool haltOnAssignment, bool doComp);
  int localSearchProbeSaps(const char *outfile, bool doComp);
  int localSearchProbeGsat(const char *outfile, bool doComp);
  int lobjoisProbe(bool haltOnAssignment, bool doComp);

  inline double *getFeatureVals() { return featureVals; }
  inline char **getFeatureNames() { return featureNames; }

  inline int getNumVals() { return numActiveVars; }
  inline int getNumClaus() { return numActiveClauses; }

  inline int getNumFeatures() { return indexCount; }

  void writeFeature(const char *name, double val);

  void writeFeaturesToFile(const char *name); // write comma separated features to passed file
  void writeFeaturesToFile(FILE *f);

  void writeFeatNamesToFile(const char *name); // write comma separated feature names to passed file
  void writeFeatNamesToFile(FILE *f);

  void outputActiveFeat(bool *active);

  void testAPI();
  void testBackTrack();
  bool stupidSearch();

  void start_computation(bool preprocessor_solved, float pre_time);
  void finish_computation();
  bool ignoreBadFeats;
  char *vlineFilename;
  map<string, int> nameToIndex; // maps feature names to index in featureVals array

  long seed; // - random seed
  bool solved;
};

#endif
