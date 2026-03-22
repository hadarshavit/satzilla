#include "SATinstance.h"
#include "global.h"

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include <sys/times.h>
#include <limits.h>
#include <stdlib.h>

#include <set>
#include <map>
#include <queue>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <cassert>
#include <stdexcept>
#include "lp_solve_5.0/lpkit.h"
#include "lp_solve_5.0/patchlevel.h"
#include "stopwatch.h"

#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X, Y) ((X) > (Y) ? (Y) : (X))
#define ABS(X) ((X) > 0 ? (X) : -(X))

#define positive(X) ((X) > 0)
#define negative(X) ((X) < 0)

#define RESERVED_VALUE (-512)

const char *SATinstance::badFeatNames[] = {"KLB-featuretime", "lpTIME", "CG-featuretime", "lobjois-featuretime", "UnitProbeTIME", "POSNEG-RATIO-CLAUSE-max", "UNARY", "VCG-VAR-mean", "VCG-VAR-coeff-variation", "VCG-CLAUSE-min", "POSNEG-RATIO-CLAUSE-min", "POSNEG-RATIO-VAR-mean", "POSNEG-RATIO-VAR-min", "POSNEG-RATIO-VAR-entropy", "HORNY-VAR-mean", "VG-max", "CG-mean", "CG-max", "CG-entropy", "vars-reduced-depth-256", "gsat_BestStep_Median", "gsat_BestStep_Q.10", "gsat_BestStep_Q.90", "gsat_FirstLMRatio_CoeffVariance", "gsat_FirstLMRatio_Mean", "saps_BestStep_Median", "saps_BestStep_Q.10", "saps_BestStep_Q.90", "gsat_totaltime", "saps_totaltime", "cluster-coeff-min", "LPSLack-min", "LPSLack-max", "gsat_AvgImproveToBest_Mean"};
const int SATinstance::numBadFeats = 34;

// == THIS IS HORRIBLE STYLE! But the easiest way to interface
// == with ubcsat for now
SATinstance *currInstanceForUBC = 0;

namespace {
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

struct TranslatedFormula
{
	std::vector<std::vector<int> > clauses;
	std::vector<int> clauseLiteralOccurrences;
	std::vector<int> literalOccurrences;
	int numVars = 0;
};

struct SequenceStats
{
	double min = 0.0;
	double max = 0.0;
	double mode = 0.0;
	double mean = 0.0;
	double stddev = 0.0;
	double zeros = 0.0;
	double entropy = 0.0;
	double q1 = 0.0;
	double q2 = 0.0;
	double q3 = 0.0;
	double valRate = 0.0;
};

struct EdgeKey
{
	int a;
	int b;

	EdgeKey(int lhs = 0, int rhs = 0)
	{
		if (lhs <= rhs)
		{
			a = lhs;
			b = rhs;
		}
		else
		{
			a = rhs;
			b = lhs;
		}
	}

	bool operator<(const EdgeKey &other) const
	{
		if (a != other.a)
			return a < other.a;
		return b < other.b;
	}
};

struct DirectedEdgeKey
{
	int from;
	int to;

	DirectedEdgeKey(int lhs = 0, int rhs = 0)
		: from(lhs), to(rhs)
	{
	}

	bool operator<(const DirectedEdgeKey &other) const
	{
		if (from != other.from)
			return from < other.from;
		return to < other.to;
	}
};

double quantileCunnane(const std::vector<double> &sortedValues, double probability)
{
	if (sortedValues.empty())
		return 0.0;
	if (sortedValues.size() == 1)
		return sortedValues[0];

	const double alphap = 0.4;
	const double betap = 0.4;
	const double aleph = sortedValues.size() * probability + alphap + probability * (1.0 - alphap - betap);
	const long j = static_cast<long>(floor(aleph));
	const double gamma = aleph - static_cast<double>(j);

	if (j <= 0)
		return sortedValues.front();
	if (static_cast<size_t>(j) >= sortedValues.size())
		return sortedValues.back();

	return (1.0 - gamma) * sortedValues[j - 1] + gamma * sortedValues[j];
}

SequenceStats computeSequenceStats(const std::vector<double> &values)
{
	SequenceStats stats;
	stats.zeros = 0.0;

	std::vector<double> normalized;
	normalized.reserve(values.size());
	for (size_t i = 0; i < values.size(); ++i)
	{
		if (fabs(values[i]) < EPSILON)
			stats.zeros += 1.0;
		else
			normalized.push_back(values[i]);
	}

	if (normalized.empty())
		normalized.push_back(0.0);

	sort(normalized.begin(), normalized.end());

	stats.min = normalized.front();
	stats.max = normalized.back();
	stats.mean = std::accumulate(normalized.begin(), normalized.end(), 0.0) / static_cast<double>(normalized.size());

	double variance = 0.0;
	for (size_t i = 0; i < normalized.size(); ++i)
		variance += (normalized[i] - stats.mean) * (normalized[i] - stats.mean);
	stats.stddev = sqrt(variance / static_cast<double>(normalized.size()));

	size_t bestCount = 0;
	size_t distinctCount = 0;
	size_t runCount = 0;
	double currentValue = normalized.front();
	for (size_t i = 0; i < normalized.size(); ++i)
	{
		if (i == 0 || fabs(normalized[i] - currentValue) < EPSILON)
		{
			runCount++;
		}
		else
		{
			distinctCount++;
			if (runCount > bestCount)
			{
				bestCount = runCount;
				stats.mode = currentValue;
			}
			currentValue = normalized[i];
			runCount = 1;
		}
	}
	distinctCount++;
	if (runCount > bestCount)
	{
		bestCount = runCount;
		stats.mode = currentValue;
	}

	double sumLogCounts = 0.0;
	runCount = 0;
	currentValue = normalized.front();
	for (size_t i = 0; i < normalized.size(); ++i)
	{
		if (i == 0 || fabs(normalized[i] - currentValue) < EPSILON)
		{
			runCount++;
		}
		else
		{
			sumLogCounts += static_cast<double>(runCount) * log(static_cast<double>(runCount));
			currentValue = normalized[i];
			runCount = 1;
		}
	}
	sumLogCounts += static_cast<double>(runCount) * log(static_cast<double>(runCount));

	const double sequenceSize = static_cast<double>(normalized.size());
	stats.entropy = log(sequenceSize) - (sumLogCounts / sequenceSize);
	stats.q1 = quantileCunnane(normalized, 0.25);
	stats.q2 = quantileCunnane(normalized, 0.50);
	stats.q3 = quantileCunnane(normalized, 0.75);
	stats.valRate = static_cast<double>(distinctCount) / sequenceSize;
	return stats;
}

void writeSequenceStats(SATinstance *sat, const std::string &prefix, const std::vector<double> &values)
{
	const SequenceStats stats = computeSequenceStats(values);
	sat->writeFeature((prefix + "min").c_str(), stats.min);
	sat->writeFeature((prefix + "max").c_str(), stats.max);
	sat->writeFeature((prefix + "mode").c_str(), stats.mode);
	sat->writeFeature((prefix + "mean").c_str(), stats.mean);
	sat->writeFeature((prefix + "std").c_str(), stats.stddev);
	sat->writeFeature((prefix + "zeros").c_str(), stats.zeros);
	sat->writeFeature((prefix + "entropy").c_str(), stats.entropy);
	sat->writeFeature((prefix + "q1").c_str(), stats.q1);
	sat->writeFeature((prefix + "q2").c_str(), stats.q2);
	sat->writeFeature((prefix + "q3").c_str(), stats.q3);
	sat->writeFeature((prefix + "val_rate").c_str(), stats.valRate);
}

void writeReservedSequenceStats(SATinstance *sat, const std::string &prefix)
{
	static const char *kSuffixes[] = {
		"min", "max", "mode", "mean", "std", "zeros", "entropy", "q1", "q2", "q3", "val_rate"};
	for (size_t i = 0; i < sizeof(kSuffixes) / sizeof(kSuffixes[0]); ++i)
		sat->writeFeature((prefix + kSuffixes[i]).c_str(), RESERVED_VALUE);
}

void writeMantheyGraphFeatures(SATinstance *sat,
							   const std::string &graphPrefix,
							   const std::vector<double> &nodeDegrees,
							   const std::vector<double> &weights)
{
	writeSequenceStats(sat, graphPrefix + "node_", nodeDegrees);
	writeSequenceStats(sat, graphPrefix + "weights_", weights.empty() ? std::vector<double>(1, 0.0) : weights);
}

void writeReservedMantheyGraphFeatures(SATinstance *sat, const std::string &graphPrefix)
{
	writeReservedSequenceStats(sat, graphPrefix + "node_");
	writeReservedSequenceStats(sat, graphPrefix + "weights_");
}

int literalNodeIndex(int literal, int numVars)
{
	if (literal > 0)
		return literal;
	return numVars + (-literal);
}

std::vector<std::vector<int> > initUndirectedAdjacency(int nodeCount)
{
	return std::vector<std::vector<int> >(nodeCount + 1);
}

void finalizeUndirectedAdjacency(std::vector<std::vector<int> > &adjacency)
{
	for (size_t i = 1; i < adjacency.size(); ++i)
	{
		sort(adjacency[i].begin(), adjacency[i].end());
		adjacency[i].erase(unique(adjacency[i].begin(), adjacency[i].end()), adjacency[i].end());
	}
}

std::vector<double> nodeDegreesFromUndirected(const std::vector<std::vector<int> > &adjacency, int nodeCount)
{
	std::vector<double> degrees(nodeCount, 0.0);
	for (int node = 1; node <= nodeCount; ++node)
		degrees[node - 1] = static_cast<double>(adjacency[node].size());
	return degrees;
}

std::vector<double> doubledUndirectedWeights(const std::map<EdgeKey, double> &weights)
{
	std::vector<double> values;
	values.reserve(weights.size() * 2);
	for (std::map<EdgeKey, double>::const_iterator it = weights.begin(); it != weights.end(); ++it)
	{
		values.push_back(it->second);
		values.push_back(it->second);
	}
	return values;
}

std::vector<double> directedOutDegrees(const std::vector<std::vector<int> > &adjacency, int nodeCount)
{
	std::vector<double> degrees(nodeCount, 0.0);
	for (int node = 1; node <= nodeCount; ++node)
		degrees[node - 1] = static_cast<double>(adjacency[node].size());
	return degrees;
}

std::vector<double> directedWeights(const std::map<DirectedEdgeKey, double> &weights)
{
	std::vector<double> values;
	values.reserve(weights.size());
	for (std::map<DirectedEdgeKey, double>::const_iterator it = weights.begin(); it != weights.end(); ++it)
		values.push_back(it->second);
	return values;
}

std::vector<int> bfsWithinRadius(const std::vector<std::vector<int> > &adjacency, int start, int radius)
{
	std::vector<int> nodes;
	std::vector<int> distance(adjacency.size(), -1);
	std::queue<int> pending;
	pending.push(start);
	distance[start] = 0;

	while (!pending.empty())
	{
		const int node = pending.front();
		pending.pop();
		nodes.push_back(node);
		if (distance[node] == radius)
			continue;
		for (size_t i = 0; i < adjacency[node].size(); ++i)
		{
			const int neighbor = adjacency[node][i];
			if (distance[neighbor] != -1)
				continue;
			distance[neighbor] = distance[node] + 1;
			pending.push(neighbor);
		}
	}

	return nodes;
}

double regressionSlope(const std::vector<double> &xs, const std::vector<double> &ys)
{
	if (xs.empty() || xs.size() != ys.size())
		return 1.0;

	const double sx = std::accumulate(xs.begin(), xs.end(), 0.0);
	const double sy = std::accumulate(ys.begin(), ys.end(), 0.0);
	double sxx = 0.0;
	double sxy = 0.0;
	for (size_t i = 0; i < xs.size(); ++i)
	{
		sxx += xs[i] * xs[i];
		sxy += xs[i] * ys[i];
	}

	const double denom = sx * sx - static_cast<double>(xs.size()) * sxx;
	if (fabs(denom) < EPSILON)
		return 1.0;
	return (sx * sy - static_cast<double>(xs.size()) * sxy) / denom;
}

double powerLawC(int x, int xmin, double alpha)
{
	const int maxIterations = 10000;
	double num = 0.0;
	double den = 0.0;
	int i = xmin;

	if (xmin < 25)
	{
		while (i < x)
		{
			den += pow(static_cast<double>(i), alpha);
			++i;
		}

		double pOld = -2.0;
		double p = -1.0;
		int iterations = 0;
		while (fabs(p - pOld) > 1e-8 && iterations < maxIterations)
		{
			den += pow(static_cast<double>(i), alpha);
			num += pow(static_cast<double>(i), alpha);
			++i;
			++iterations;
			pOld = p;
			p = num / den;
		}

		if (iterations < maxIterations)
			return p;
	}

	return pow(static_cast<double>(x) / static_cast<double>(xmin), alpha + 1.0);
}

double mostLikelyAlpha(const std::vector<int> &x,
					   const std::vector<double> &y,
					   const std::vector<double> &syLogX,
					   const std::vector<double> &syX,
					   int maxXmin = 10)
{
	double bestAlpha = 0.0;
	double bestDiff = 1.0;
	const int n = static_cast<int>(x.size());

	for (int ind = 1; ind <= maxXmin; ++ind)
	{
		if (ind >= n - 3)
			continue;

		const int xmin = x[ind];
		const double alpha = -1.0 - (1.0 / ((syLogX[ind] / y[ind]) - log(static_cast<double>(xmin) - 0.5)));
		double worstDiff = -1.0;

		for (int j = ind + 1; j < n; ++j)
		{
			const double diff = fabs(y[j] / y[ind] - powerLawC(x[j], xmin, alpha));
			if (diff >= bestDiff)
			{
				worstDiff = diff;
				break;
			}
			worstDiff = MAX(worstDiff, diff);
		}

		for (int j = ind; j < n - 1; ++j)
		{
			if (x[j] + 1 >= x[j + 1])
				continue;
			const double diff = fabs(y[j + 1] / y[ind] - powerLawC(x[j] + 1, xmin, alpha));
			if (diff >= bestDiff)
			{
				worstDiff = diff;
				break;
			}
			worstDiff = MAX(worstDiff, diff);
		}

		if (worstDiff < bestDiff)
		{
			bestDiff = worstDiff;
			bestAlpha = alpha;
		}
	}

	return -bestAlpha;
}

int countConnectedComponents(const std::vector<std::vector<int> > &adjacency, int nodeCount)
{
	std::vector<bool> visited(nodeCount + 1, false);
	int componentCount = 0;
	for (int node = 1; node <= nodeCount; ++node)
	{
		if (visited[node])
			continue;
		componentCount++;
		std::queue<int> pending;
		pending.push(node);
		visited[node] = true;
		while (!pending.empty())
		{
			const int current = pending.front();
			pending.pop();
			for (size_t i = 0; i < adjacency[current].size(); ++i)
			{
				const int neighbor = adjacency[current][i];
				if (visited[neighbor])
					continue;
				visited[neighbor] = true;
				pending.push(neighbor);
			}
		}
	}
	return componentCount;
}

std::vector<int> burningByNodeDegree(const std::vector<std::vector<int> > &adjacency, int nodeCount)
{
	std::vector<int> coverCounts(nodeCount, 0);
	if (nodeCount <= 1)
		return coverCounts;

	coverCounts[1] = nodeCount;
	std::vector<std::pair<int, int> > nodeDegrees;
	nodeDegrees.reserve(nodeCount);
	for (int node = 1; node <= nodeCount; ++node)
		nodeDegrees.push_back(std::make_pair(node, static_cast<int>(adjacency[node].size())));
	sort(nodeDegrees.begin(), nodeDegrees.end(), [](const std::pair<int, int> &lhs, const std::pair<int, int> &rhs) {
		if (lhs.second != rhs.second)
			return lhs.second > rhs.second;
		return lhs.first < rhs.first;
	});

	const int componentCount = countConnectedComponents(adjacency, nodeCount);
	const int maxRadius = MIN(16, nodeCount - 1);
	for (int radius = 1; radius <= maxRadius; ++radius)
	{
		if (coverCounts[radius - 1] <= componentCount)
			continue;

		std::vector<bool> burned(nodeCount + 1, false);
		int circles = 0;
		while (true)
		{
			int centre = 0;
			for (size_t i = 0; i < nodeDegrees.size(); ++i)
			{
				if (!burned[nodeDegrees[i].first])
				{
					centre = nodeDegrees[i].first;
					break;
				}
			}
			if (centre == 0)
				break;

			const std::vector<int> covered = bfsWithinRadius(adjacency, centre, radius - 1);
			for (size_t i = 0; i < covered.size(); ++i)
				burned[covered[i]] = true;
			circles++;
		}
		coverCounts[radius] = circles;
	}

	return coverCounts;
}

double fractalSlope(const std::vector<int> &coverCounts)
{
	std::vector<double> xs;
	std::vector<double> ys;
	for (size_t i = 0; i < coverCounts.size(); ++i)
	{
		if (coverCounts[i] <= 0)
			continue;
		xs.push_back(log(static_cast<double>(xs.size() + 1)));
		ys.push_back(log(static_cast<double>(coverCounts[i])));
	}
	if (xs.empty())
		return 0.0;
	return -regressionSlope(xs, ys);
}

double modularityFromPartition(const std::vector<std::map<int, double> > &adjacency,
							   const std::vector<int> &community,
							   const std::vector<double> &nodeWeights,
							   double totalWeightTwice)
{
	if (totalWeightTwice <= 0.0)
		return 0.0;

	std::map<int, double> internal;
	std::map<int, double> totals;
	for (size_t node = 1; node < adjacency.size(); ++node)
	{
		const int comm = community[node];
		totals[comm] += nodeWeights[node];
		for (std::map<int, double>::const_iterator it = adjacency[node].begin(); it != adjacency[node].end(); ++it)
		{
			if (community[it->first] == comm)
				internal[comm] += it->second;
		}
	}

	double modularity = 0.0;
	for (std::map<int, double>::const_iterator it = totals.begin(); it != totals.end(); ++it)
	{
		const double inWeight = internal[it->first] / totalWeightTwice;
		const double total = it->second / totalWeightTwice;
		modularity += inWeight - total * total;
	}
	return modularity;
}

double approximateLouvainModularity(const std::vector<std::map<int, double> > &adjacency, int nodeCount)
{
	std::vector<double> nodeWeights(nodeCount + 1, 0.0);
	double totalWeightTwice = 0.0;
	for (int node = 1; node <= nodeCount; ++node)
	{
		for (std::map<int, double>::const_iterator it = adjacency[node].begin(); it != adjacency[node].end(); ++it)
			nodeWeights[node] += it->second;
		totalWeightTwice += nodeWeights[node];
	}
	if (totalWeightTwice <= 0.0)
		return 0.0;

	std::vector<int> community(nodeCount + 1, 0);
	std::vector<double> communityTotals(nodeCount + 1, 0.0);
	for (int node = 1; node <= nodeCount; ++node)
	{
		community[node] = node;
		communityTotals[node] = nodeWeights[node];
	}

	bool changed = true;
	while (changed)
	{
		changed = false;
		for (int node = 1; node <= nodeCount; ++node)
		{
			std::map<int, double> weightsToCommunity;
			for (std::map<int, double>::const_iterator it = adjacency[node].begin(); it != adjacency[node].end(); ++it)
				weightsToCommunity[community[it->first]] += it->second;

			const int currentCommunity = community[node];
			communityTotals[currentCommunity] -= nodeWeights[node];

			int bestCommunity = currentCommunity;
			double bestGain = 0.0;
			for (std::map<int, double>::const_iterator it = weightsToCommunity.begin(); it != weightsToCommunity.end(); ++it)
			{
				const double gain = it->second - (nodeWeights[node] * communityTotals[it->first] / totalWeightTwice);
				if (gain > bestGain + 1e-12)
				{
					bestGain = gain;
					bestCommunity = it->first;
				}
			}

			community[node] = bestCommunity;
			communityTotals[bestCommunity] += nodeWeights[node];
			if (bestCommunity != currentCommunity)
				changed = true;
		}
	}

	return modularityFromPartition(adjacency, community, nodeWeights, totalWeightTwice);
}
}

void writeFeature(const char *name, double val)
{
	currInstanceForUBC->writeFeature(name, val);
}

SATinstance::SATinstance(const char *filename, bool doComp, long _seed)
	: numVars(0),
	  numClauses(0),
	  clauses(nullptr),
	  negClausesWithVar(nullptr),
	  posClausesWithVar(nullptr),
	  varStates(nullptr),
	  clauseStates(nullptr),
	  numActiveVars(0),
	  numActiveClauses(0),
	  numActiveClausesWithVar(nullptr),
	  numBinClausesWithVar(nullptr),
	  clauseLengths(nullptr),
	  active_clause(nullptr),
	  currentClause(0),
	  currentLit(0),
	  currentLit2(0),
	  currentClauseForLitIter(nullptr),
	  currentClauseForLitIter2(nullptr),
	  currentClauseWithVar(0),
	  currentVarForClauseIter(0),
	  posClauses(false),
	  dia_clause_table(nullptr),
	  dstack(nullptr),
	  var_depth(nullptr),
	  clause_stamps(nullptr),
	  CCS(0),
	  dia_clause_list(nullptr),
	  featureNames{},
	  indexCount(0),
	  featureVals{},
	  cgTimeout(false),
	  test_flag(nullptr),
	  lp_return_val(0),
	  inputFileName(filename),
	  ignoreBadFeats(false),
	  vlineFilename(nullptr),
	  seed(_seed),
	  solved(false)
{
	currInstanceForUBC = this;

	if (doComp)
	{
		ifstream infile(filename);
		if (!infile)
		{
			fprintf(stderr, "c Error: Could not read from input file %s.\n", filename);
			throw std::runtime_error(std::string("Could not read from input file ") + filename);
		}

		inputFileName = filename;
		char chbuf;
		char strbuf[1024];
		infile.get(chbuf);
		while (chbuf != 'p')
		{
			//    infile.getline(strbuf, 100);
			infile.ignore(1000, '\n');
			infile.get(chbuf);
			if (!infile)
			{
				fprintf(stderr, "c ERROR: Premature EOF reached in %s\n", filename);
				throw std::runtime_error(std::string("Premature EOF reached in ") + filename);
			}
		}

		infile >> strbuf; // "cnf"
		if (strcmp(strbuf, "cnf") != 0)
		{
			fprintf(stderr, "c Error: Can only understand cnf format!\n");
			throw std::runtime_error("Can only understand cnf format");
		}
		// == TODO: During parsing should really skip comment lines....
		// == TODO: Parser should really check for EOF

		infile >> numVars >> numClauses;

		clauses = new int *[numClauses];
		clauseStates = new clauseState[numClauses];
		clauseLengths = new int[numClauses];

		negClausesWithVar = new vector<int>[numVars + 1];
		posClausesWithVar = new vector<int>[numVars + 1];

		numActiveClausesWithVar = new int[numVars + 1];
		numBinClausesWithVar = new int[numVars + 1];

		for (int i = 1; i <= numVars; i++)
		{
			numActiveClausesWithVar[i] = 0;
			numBinClausesWithVar[i] = 0;
		}

		int *lits = new int[numVars + 1];

		// read stuff into data structure.
		// take care of all data structure
		// clauseLengths,  unitClauses, numBinClausesWithVar
		// negClausesWithVar, posClausesWithVar, numActiveClausesWithVar
		for (int clauseNum = 0; clauseNum < numClauses; clauseNum++)
		{
			int numLits = 0; // not including 0 terminator
			if (!infile)
			{
				fprintf(stderr, "c ERROR: Premature EOF reached in %s\n", filename);
				throw std::runtime_error(std::string("Premature EOF reached in ") + filename);
			}

			infile >> lits[numLits];
			while (lits[numLits] != 0)
			{
				infile >> lits[++numLits];
			}

			/* test if some literals are redundant and sort the clause */
			bool tautology = false;
			for (int i = 0; i < numLits - 1; i++)
			{
				int tempLit = lits[i];
				for (int j = i + 1; j < numLits; j++)
				{
					if (ABS(tempLit) > ABS(lits[j]))
					{
						// this is sorting the literals
						int temp = lits[j];
						lits[j] = tempLit;
						tempLit = temp;
					}
					else if (tempLit == lits[j])
					{
						lits[j--] = lits[--numLits];
						printf("c literal %d is redundant in clause %d\n", tempLit, clauseNum);
					}
					else if (ABS(tempLit) == ABS(lits[j]))
					{
						tautology = true;
						//	  printf("c Clause %d is tautological.\n", clauseNum);
						//	  break;
					}
				}
				if (tautology)
					break;
				else
					lits[i] = tempLit;
			}

			if (!tautology)
			{
				clauseLengths[clauseNum] = numLits;
				clauses[clauseNum] = new int[numLits + 1];
				clauseStates[clauseNum] = ACTIVE;

				if (numLits == 1)
					unitClauses.push(clauseNum);
				else if (numLits == 2)
					for (int i = 0; i < numLits; i++)
						numBinClausesWithVar[ABS(lits[i])]++;

				for (int litNum = 0; litNum < numLits; litNum++)
				{
					if (lits[litNum] < 0)
						negClausesWithVar[ABS(lits[litNum])].push_back(clauseNum);
					else
						posClausesWithVar[lits[litNum]].push_back(clauseNum);
					numActiveClausesWithVar[ABS(lits[litNum])]++;
					clauses[clauseNum][litNum] = lits[litNum];
				}
				clauses[clauseNum][numLits] = 0;
			}
			else
			{
				clauseNum--;
				numClauses--;
			}
			//    printf("%d, %d, %d, %d %d\n", clauses[clauseNum][0],clauses[clauseNum][1],clauses[clauseNum][2],clauseNum, numClauses);
		}
		delete[] lits;
		numActiveClauses = numClauses;

		// remove some redandant variables
		// prepar data sturcuture: varStates
		varStates = new varState[numVars + 1];
		numActiveVars = numVars;
		for (int i = 1; i <= numVars; i++)
		{
			if (numActiveClausesWithVar[i] == 0)
			{
				varStates[i] = IRRELEVANT;
				numActiveVars--;
			}
			else
				varStates[i] = UNASSIGNED;
		}

		// before doing anything first do a round of unit propogation to remove all the
		// unit clasues
		int dummy1, dummy2;
		unitprop(dummy1, dummy2);

		test_flag = new int[numVars + 1];
		indexCount = 0;

		if (seed == 0)
			seed = (long)time(NULL);

		srand(seed);
		if (DEB)
			printf("c Number of variabe is: %d, Number of clause is : %d \n", numActiveVars, numActiveClauses);
	}
}

SATinstance::~SATinstance()
{
	delete[] varStates;
	for (int i = 0; i < numClauses; i++)
		delete[] clauses[i];
	delete[] clauses;
	delete[] clauseStates;
	delete[] clauseLengths;

	delete[] negClausesWithVar;
	delete[] posClausesWithVar;

	delete[] numActiveClausesWithVar;
	delete[] numBinClausesWithVar;

	delete[] test_flag;
}

inline vector<int> &SATinstance::clausesWithLit(int lit)
{
	if (positive(lit))
		return posClausesWithVar[lit];
	else
		return negClausesWithVar[-lit];
}

bool SATinstance::reduceClauses(int lit, int &numClausesReduced, int &numVarsReduced)
{

	// "remove" vars from inconsistent clauses
	for (int i = 0; i < (int)clausesWithLit(-lit).size(); i++)
	{
		int clause = clausesWithLit(-lit)[i];
		if (clauseStates[clause] == ACTIVE)
		{
			reducedClauses.push(clause);
			numClausesReduced++;

			clauseLengths[clause]--;
			if (clauseLengths[clause] == 2)
				for (int i = 0; clauses[clause][i] != 0; i++)
					numBinClausesWithVar[ABS(clauses[clause][i])]++;
			else if (clauseLengths[clause] == 1)
			{
				for (int i = 0; clauses[clause][i] != 0; i++)
					numBinClausesWithVar[ABS(clauses[clause][i])]--;
				unitClauses.push(clause);
			}
			else if (clauseLengths[clause] == 0)
				return false;
		}
	}

	// satisfy consistent clauses
	for (int i = 0; i < (int)clausesWithLit(lit).size(); i++)
	{
		int clause = clausesWithLit(lit)[i];
		if (clauseStates[clause] == ACTIVE)
		{

			clauseStates[clause] = PASSIVE;
			reducedClauses.push(clause);
			numActiveClauses--;

			int j = 0;
			int otherVarInClause = ABS(clauses[clause][j]);
			while (otherVarInClause != 0)
			{
				numActiveClausesWithVar[otherVarInClause]--;
				if (clauseLengths[clause] == 2)
					numBinClausesWithVar[otherVarInClause]--;

				// is the var now irrelevant (active, but existing in no clauses)?
				if (numActiveClausesWithVar[otherVarInClause] == 0 &&
					varStates[otherVarInClause] == UNASSIGNED)
				{
					varStates[otherVarInClause] = IRRELEVANT;
					reducedVars.push(otherVarInClause);
					numActiveVars--;

					numVarsReduced++;
				}

				j++;
				otherVarInClause = ABS(clauses[clause][j]);
			}
			numClausesReduced++;
		}
	}
	return true;
}

bool SATinstance::setVarAndProp(int var, bool val)
{
	int numClausesReduced = 0;
	int numVarsReduced = 1;

	assert(varStates[var] == UNASSIGNED);
	varStates[var] = val ? TRUE_VAL : FALSE_VAL;
	reducedVars.push(var);
	numActiveVars--;

	int lit = val ? var : -var;
	bool consistent = reduceClauses(lit, numClausesReduced, numVarsReduced);

	if (consistent)
		consistent = unitprop(numClausesReduced, numVarsReduced);

	numReducedClauses.push(numClausesReduced);
	numReducedVars.push(numVarsReduced);

	return consistent;
}

bool SATinstance::unitprop(int &numClausesReduced, int &numVarsReduced)
{

	bool consistent = true;

	while (!unitClauses.empty() && consistent)
	{
		int clauseNum = unitClauses.top();
		unitClauses.pop();

		if (clauseStates[clauseNum] != ACTIVE)
			continue;

		int litNum = 0;
		while (varStates[ABS(clauses[clauseNum][litNum])] != UNASSIGNED)
		{
			litNum++;
		}

		// assertions are our friends!
		assert(clauseLengths[clauseNum] == 1);

		int lit = clauses[clauseNum][litNum];

		varStates[ABS(lit)] = positive(lit) ? TRUE_VAL : FALSE_VAL;
		reducedVars.push(ABS(lit));
		numActiveVars--;
		numVarsReduced++;

		consistent &= reduceClauses(lit, numClausesReduced, numVarsReduced);
	}

	return consistent;
}

void SATinstance::backtrack()
{
	int numVarsReduced = numReducedVars.top();
	numReducedVars.pop();
	for (int i = 0; i < numVarsReduced; i++)
	{
		int var = reducedVars.top();
		reducedVars.pop();
		varStates[var] = UNASSIGNED;
		numActiveVars++;
	}

	int numClausesReduced = numReducedClauses.top();
	numReducedClauses.pop();
	for (int i = 0; i < numClausesReduced; i++)
	{
		int clause = reducedClauses.top();
		reducedClauses.pop();

		if (clauseStates[clause] != ACTIVE)
		{
			numActiveClauses++;
			clauseStates[clause] = ACTIVE;

			if (clauseLengths[clause] == 2)
				for (int j = 0; clauses[clause][j] != 0; j++)
				{
					numActiveClausesWithVar[ABS(clauses[clause][j])]++;
					numBinClausesWithVar[ABS(clauses[clause][j])]++;
				}
			else
				for (int j = 0; clauses[clause][j] != 0; j++)
					numActiveClausesWithVar[ABS(clauses[clause][j])]++;
		}
		else
		{
			clauseLengths[clause]++;
			if (clauseLengths[clause] == 2)
				for (int j = 0; clauses[clause][j] != 0; j++)
					numBinClausesWithVar[ABS(clauses[clause][j])]++;

			else if (clauseLengths[clause] == 3)
				for (int j = 0; clauses[clause][j] != 0; j++)
					numBinClausesWithVar[ABS(clauses[clause][j])]--;
		}
	}

	while (!unitClauses.empty())
		unitClauses.pop();
}

void SATinstance::outputAssignment()
{
	FILE *vlineFile;
	vlineFilename = new char[512];
	if (!createTempPath(vlineFilename, 512, "/XXXXXX"))
	{
		fprintf(stderr, "c Couldn't create temp file path\n");
		delete[] vlineFilename;
		return;
	}
	printf("Assignment file name: %s", vlineFilename);
	if ((vlineFile = fopen(vlineFilename, "w")) == NULL)
	{
		fprintf(stderr, "c Couldn't open temp file\n");
		delete[] vlineFilename;
		return;
	}

	assert(numActiveClauses == 0 && numActiveVars == 0);
	fprintf(vlineFile, "s SATISFIABLE\n");
	fprintf(vlineFile, "v ");
	for (int i = 1; i <= numVars; i++)
	{
		switch (varStates[i])
		{
		case TRUE_VAL:
			fprintf(vlineFile, "%d ", i);
			break;
		case FALSE_VAL:
			fprintf(vlineFile, "-%d ", i);
			break;
		case IRRELEVANT:
			// fprintf(vlineFile, "%d ", i);
			break;
		default:
			fprintf(stderr, "c Error: outputAssignment called before solution reached.  Var %d is unassigned.\n", i);
		}
	}
	fprintf(vlineFile, "0\n");
	fclose(vlineFile);
	// delete[] vlineFilename;  *****done later
}

bool SATinstance::stupidSearch()
{
	int i = 1;

	while (varStates[i] != UNASSIGNED && i <= numVars)
		i++;

	assert(i <= numVars);
	assert(unitClauses.size() == 0);

	print();

	if (setVarAndProp(i, true))
	{
		assert(unitClauses.size() == 0);
		if (numActiveVars == 0 ||
			stupidSearch())
			return true;
	}

	print();

	backtrack();

	print();

	assert(unitClauses.size() == 0);

	if (setVarAndProp(i, false))
	{
		assert(unitClauses.size() == 0);
		if (numActiveVars == 0 ||
			stupidSearch())
			return true;
	}

	print();

	backtrack();

	print();

	assert(unitClauses.size() == 0);

	return false;
}

void SATinstance::retardedSearch()
{
	int i = 1;
	while (varStates[i] != UNASSIGNED && i <= numVars)
		i++;
	assert(i <= numVars);
	assert(unitClauses.size() == 0);

	//  printf("Entered Retard\n");
	print();

	printf("c Setting %d to true\n", i);

	if (setVarAndProp(i, true) && numActiveVars > 0)
	{
		assert(unitClauses.size() == 0);
		retardedSearch();
	}

	print();

	backtrack();

	printf("c Backtracked\n");
	print();

	assert(unitClauses.size() == 0);

	printf("c Setting %d to false\n", i);

	if (setVarAndProp(i, false) && numActiveVars > 0)
	{
		assert(unitClauses.size() == 0);
		retardedSearch();
	}

	print();

	backtrack();

	printf("c Backtracked\n");
	print();

	assert(unitClauses.size() == 0);
}

inline void p(const char *in)
{
	LogDebug("%s\n", in);
	fflush(stdout);
}

void SATinstance::start_computation(bool preprocessor_solved, float pre_time)
{
	solved = -512;
}

int SATinstance::computeFeatures(bool doComp, bool doClauseGraphFeatures)
{
	//  testBackTrack();
	// Stopwatch sw;
	if (DEB)
		p("c Initializing...");
	// sw.Start();
	if (!doComp)
	{
		double dummy[] = {RESERVED_VALUE, RESERVED_VALUE};
		writeFeature("nvarsOrig", (double)OrigNumVars);
		writeFeature("nclausesOrig", (double)OrigNumClauses);
		writeFeature("nvars", 0);
		writeFeature("nclauses", 0);
		writeFeature("reducedVars", RESERVED_VALUE);
		writeFeature("reducedClauses", RESERVED_VALUE);
		writeFeature("Pre-featuretime", preTime);
		writeFeature("vars-clauses-ratio", RESERVED_VALUE);
		writeStats(dummy, 2, "POSNEG-RATIO-CLAUSE");
		writeFeature("POSNEG-RATIO-CLAUSE-entropy", RESERVED_VALUE);
		writeStats(dummy, 2, "VCG-CLAUSE");
		writeFeature("VCG-CLAUSE-entropy", RESERVED_VALUE);
		writeFeature("UNARY", RESERVED_VALUE);
		writeFeature("BINARY+", RESERVED_VALUE);
		writeFeature("TRINARY+", RESERVED_VALUE);
		writeFeature("Basic-featuretime", RESERVED_VALUE);
		writeStats(dummy, 2, "VCG-VAR");
		writeFeature("VCG-VAR-entropy", RESERVED_VALUE);
		writeStatsSTDEV(dummy, 2, "POSNEG-RATIO-VAR");
		writeFeature("POSNEG-RATIO-VAR-entropy", RESERVED_VALUE);
		writeStats(dummy, 2, "HORNY-VAR");
		writeFeature("HORNY-VAR-entropy", RESERVED_VALUE);
		writeFeature("horn-clauses-fraction", RESERVED_VALUE);
		writeStats(dummy, 2, "VG");
		writeFeature("KLB-featuretime", RESERVED_VALUE);
		writeStats(dummy, 2, "CG");
		writeFeature("CG-entropy", RESERVED_VALUE);
		writeStats(dummy, 2, "cluster-coeff");
		writeFeature("cluster-coeff-entropy", RESERVED_VALUE);
		writeFeature("CG-featuretime", RESERVED_VALUE);
		return 0;
	}

	// fprintf(stderr, "Computing features. Prefix %s endpref\n", featurePrefix);
	//  node degree stats for var-clause graph
	int *var_array = new int[numVars + 1];
	int *var_graph = new int[numVars + 1];
	bool *var_graph_found = new bool[numVars + 1];
	double *var_graph_norm = new double[numVars + 1];
	int *horny_var = new int[numVars + 1];
	double *horny_var_norm = new double[numVars + 1];
	double *var_array_norm = new double[numVars + 1];
	int *clause_array = new int[numClauses];
	double *clause_array_norm = new double[numClauses];
	int *pos_in_clause = new int[numClauses];
	int *neg_in_clause = new int[numClauses];
	double *pos_frac_in_clause = new double[numClauses];
	int *pos_var = new int[numVars + 1];
	int *neg_var = new int[numVars + 1];
	double *pos_frac_per_var = new double[numVars + 1];
	int unary = 0, binary = 0, trinary = 0;
	int horn_clauses = 0;
	int t, tt;

	// initialize
	for (t = 1; t <= numVars; t++)
	{
		var_array[t] = 0;
		pos_var[t] = 0;
		neg_var[t] = 0;
		horny_var[t] = 0;
		var_array_norm[t] = RESERVED_VALUE;
		pos_frac_per_var[t] = RESERVED_VALUE;
	}

	for (t = 0; t < numClauses; t++)
	{
		clause_array[t] = (int)RESERVED_VALUE;
		clause_array_norm[t] = RESERVED_VALUE;
		pos_in_clause[t] = (int)RESERVED_VALUE;
		neg_in_clause[t] = (int)RESERVED_VALUE;
		pos_frac_in_clause[t] = RESERVED_VALUE;
	}
	if (DEB)
		p("c Go through clauses...");
	writeFeature("nvarsOrig", (double)OrigNumVars);
	writeFeature("nclausesOrig", (double)OrigNumClauses);
	writeFeature("nvars", (double)numActiveVars);
	writeFeature("nclauses", (double)numActiveClauses);
	if ((double)numActiveVars == 0)
	{
		writeFeature("reducedVars", RESERVED_VALUE);
		writeFeature("reducedClauses", RESERVED_VALUE);
		writeFeature("Pre-featuretime", preTime);
		writeFeature("vars-clauses-ratio", RESERVED_VALUE);
	}
	else
	{
		writeFeature("reducedVars", ((double)OrigNumVars - (double)numActiveVars) / (double)numActiveVars);
		writeFeature("reducedClauses", ((double)OrigNumClauses - (double)numActiveClauses) / (double)numActiveClauses);
		writeFeature("Pre-featuretime", preTime);
		writeFeature("vars-clauses-ratio", ((double)numActiveVars) / (double)numActiveClauses);
	}
	if (numActiveVars == 0 || numActiveClauses == 0)
	{
		solved = 1;
	}
	// go through all the clauses
	// What we get from here is
	// clause_array : number of lierals
	// pos_in_clause/neg_in_clause
	// var_array: number of cluses contain this variable
	// pos_var/neg_var
	int *clause, lit;
	t = 0;
	for (clause = firstClause(); clause != NULL; clause = nextClause())
	{
		// initialize
		clause_array[t] = 0;
		pos_in_clause[t] = 0;
		neg_in_clause[t] = 0;

		for (lit = firstLitInClause(clause); lit != 0; lit = nextLitInClause())
		{
			clause_array[t]++;
			var_array[ABS(lit)]++;

			if (positive(lit))
			{
				pos_in_clause[t]++;
				pos_var[ABS(lit)]++;
			}
			else
			{
				neg_in_clause[t]++;
				neg_var[ABS(lit)]++;
			}
		}

		// may be this is a bad name for this.
		// basically, it compute the bias for the assignment
		// do we need say anything for cluase_array[t]=0
		// this should not happened
		if (clause_array[t] != 0)
			pos_frac_in_clause[t] = 2.0 * fabs(0.5 - (double)pos_in_clause[t] / ((double)pos_in_clause[t] + (double)neg_in_clause[t]));
		else
		{
			pos_frac_in_clause[t] = RESERVED_VALUE;
			//	  fprintf(stderr, "L %d clause %d empty\n", featureLevel, t);
		}

		// cardinality
		switch (clause_array[t])
		{
		case 1:
			unary++;
			break;
		case 2:
			binary++;
			break;
		case 3:
			trinary++;
			break;
		}

		// NOTE: isn't neg_in_clause <= 1 also horny? GMA
		// this is really not make sense. by switching pos/neg, you can get different horn clause
		// horn clause
		if (pos_in_clause[t] <= 1)
		{
			for (lit = firstLitInClause(clause); lit != 0; lit = nextLitInClause())
				horny_var[ABS(lit)]++;
			horn_clauses++;
		}
		// normalize
		clause_array_norm[t] = (double)clause_array[t] / (double)numActiveVars;
		// increment clause index
		t++;
	}
	//  fprintf(stderr, "Level %d: Went through %d clauses\n", featureLevel, t);

	// positive ratio in clauses
	writeStats(pos_frac_in_clause, numClauses, "POSNEG-RATIO-CLAUSE");
	writeFeature("POSNEG-RATIO-CLAUSE-entropy", array_entropy(pos_frac_in_clause, numClauses, 100, 1));
	// clause side in the bipartite graph
	writeStats(clause_array_norm, numClauses, "VCG-CLAUSE");
	writeFeature("VCG-CLAUSE-entropy", array_entropy(clause_array, numClauses, numActiveVars + 1));
	// cardinality of clauses
	if ((double)numActiveVars == 0)
	{
		writeFeature("UNARY", RESERVED_VALUE);
		writeFeature("BINARY+", RESERVED_VALUE);
		writeFeature("TRINARY+", RESERVED_VALUE);
	}
	else
	{
		writeFeature("UNARY", (double)unary / (double)numActiveClauses);
		writeFeature("BINARY+", (double)(unary + binary) / (double)numActiveClauses);
		writeFeature("TRINARY+", (double)(unary + binary + trinary) / (double)numActiveClauses);
	}

	writeFeature("Basic-featuretime", gSW.TotalLap() - myTime);
	if (gSW.TotalLap() -myTime > TOTAL_TIMEOUT)
	{
		double dummy[] = {RESERVED_VALUE, RESERVED_VALUE};
		writeStats(dummy, 2, "VCG-VAR");
		writeFeature("VCG-VAR-entropy", RESERVED_VALUE);
		writeStatsSTDEV(dummy, 2, "POSNEG-RATIO-VAR");
		writeFeature("POSNEG-RATIO-VAR-entropy", RESERVED_VALUE);
		writeStats(dummy, 2, "HORNY-VAR");
		writeFeature("HORNY-VAR-entropy", RESERVED_VALUE);
		writeFeature("horn-clauses-fraction", RESERVED_VALUE);
		writeStats(dummy, 2, "VG");
		writeFeature("KLB-featuretime", 0);
		writeStats(dummy, 2, "CG");
		writeFeature("CG-entropy", RESERVED_VALUE);
		writeStats(dummy, 2, "cluster-coeff");
		writeFeature("cluster-coeff-entropy", RESERVED_VALUE);
		writeFeature("CG-featuretime", 0);
		return TOTAL_TIMEOUT_CODE;
	}

	myTime = gSW.TotalLap();

	if (DEB)
		p("c Go through variables...");

	double cur_time;
	// Go through the variables
	for (t = 1; t <= numVars; t++)
	{
		if (varStates[t] != UNASSIGNED || var_array[t] == 0) // do we still want the second part?
		{
			var_graph[t] = (int)RESERVED_VALUE;
			var_array_norm[t] = RESERVED_VALUE;
			var_graph_norm[t] = RESERVED_VALUE;
			horny_var[t] = (int)RESERVED_VALUE;
			horny_var_norm[t] = RESERVED_VALUE;
			var_array[t] = (int)RESERVED_VALUE;
			pos_var[t] = (int)RESERVED_VALUE;
			neg_var[t] = (int)RESERVED_VALUE;
			pos_frac_per_var[t] = RESERVED_VALUE;
			continue;
		}

		for (tt = 1; tt <= numVars; tt++)
			var_graph_found[tt] = false;

		// now do the variable graph
		for (clause = firstClauseWithVar(t, false); clause != NULL; clause = nextClauseWithVar())
		{
			// fprintf(stderr, "Var %d false: clause %xd\n", t, clause);
			for (lit = firstLitInClause(clause); lit != 0; lit = nextLitInClause())
				var_graph_found[ABS(lit)] = true;
		}
		for (clause = firstClauseWithVar(t, true); clause != NULL; clause = nextClauseWithVar())
		{
			// fprintf(stderr, "Var %d truee: clause %xd\n", t, clause);
			for (lit = firstLitInClause(clause); lit != 0; lit = nextLitInClause())
				var_graph_found[ABS(lit)] = true;
		}

		var_graph[t] = -1; // counting self
		for (tt = 1; tt <= numVars; tt++)
			if (var_graph_found[tt])
				var_graph[t]++;

		// calculate and normalize
		pos_frac_per_var[t] = 2.0 * fabs(0.5 - (double)pos_var[t] / ((double)pos_var[t] + (double)neg_var[t]));
		var_array_norm[t] = (double)var_array[t] / (double)numActiveClauses;
		var_graph_norm[t] = (double)var_graph[t] / (double)numActiveClauses;
		horny_var_norm[t] = (double)horny_var[t] / (double)numActiveClauses;

		if (t % 100 == 0){
			if (gSW.TotalLap() - myTime > TOTAL_TIMEOUT)
			{
				double dummy[] = {RESERVED_VALUE, RESERVED_VALUE};
				writeStats(dummy, 2, "VCG-VAR");
				writeFeature("VCG-VAR-entropy", RESERVED_VALUE);
				writeStatsSTDEV(dummy, 2, "POSNEG-RATIO-VAR");
				writeFeature("POSNEG-RATIO-VAR-entropy", RESERVED_VALUE);
				writeStats(dummy, 2, "HORNY-VAR");
				writeFeature("HORNY-VAR-entropy", RESERVED_VALUE);
				writeFeature("horn-clauses-fraction", RESERVED_VALUE);
				writeStats(dummy, 2, "VG");
				writeFeature("KLB-featuretime", gSW.TotalLap() - myTime);
				writeStats(dummy, 2, "CG");
				writeFeature("CG-entropy", RESERVED_VALUE);
				writeStats(dummy, 2, "cluster-coeff");
				writeFeature("cluster-coeff-entropy", RESERVED_VALUE);
				writeFeature("CG-featuretime", 0);
				return VCG_TIMEOUT_CODE;
			}
			// if (cur_time > TOTAL_TIMEOUT){
			// 	writeFeature("KLB-featuretime", gSW.TotalLap() - myTime);
			// 	return TOTAL_TIMEOUT_CODE;
			// }
		}
	}

	// variable side in the bipartite graph
	writeStats(var_array_norm + 1, numActiveVars, "VCG-VAR");
	writeFeature("VCG-VAR-entropy", array_entropy(var_array + 1, numActiveVars, numActiveClauses + 1));

	/* == DEBUG:
	fprintf(stderr, "c L %d: %lf %lf %lf %lf\n", featureLevel, array_min(clause_array_norm, NB_CLAUSE), array_max(clause_array_norm, NB_CLAUSE), mean(clause_array_norm, NB_CLAUSE), stdev(clause_array_norm, NB_CLAUSE, mean(clause_array_norm, NB_CLAUSE)));
	for(t=0; t<NB_CLAUSE; t++)
	{
	fprintf(stderr, "c L %d clause[%d]:\t", featureLevel, t);
	if(clause_array_norm[t]==RESERVED_VALUE) fprintf(stderr, "RESERVED\n");
	else fprintf(stderr, "c %lf\n", clause_array_norm[t]);
	}
	*/

	// positive ratio in variables
	writeStatsSTDEV(pos_frac_per_var + 1, numActiveVars, "POSNEG-RATIO-VAR");
	writeFeature("POSNEG-RATIO-VAR-entropy", array_entropy(pos_frac_per_var + 1, numActiveVars, 100, 1));

	// horn clauses
	writeStats(horny_var_norm + 1, numActiveVars, "HORNY-VAR");
	writeFeature("HORNY-VAR-entropy", array_entropy(horny_var + 1, numActiveVars, numActiveClauses + 1));
	if ((double)numActiveVars == 0)
		writeFeature("horn-clauses-fraction", RESERVED_VALUE);
	else
		writeFeature("horn-clauses-fraction", (double)horn_clauses / (double)numActiveClauses);

	// variable graph
	writeStats(var_graph_norm + 1, numActiveVars, "VG");

	// clean up after yourself, you pig!
	delete[] var_array;
	delete[] var_graph;
	delete[] var_graph_norm;
	delete[] horny_var;
	delete[] horny_var_norm;
	delete[] var_array_norm;
	delete[] clause_array;
	delete[] clause_array_norm;
	delete[] pos_in_clause;
	delete[] neg_in_clause;
	delete[] pos_frac_in_clause;
	delete[] pos_var;
	delete[] neg_var;
	delete[] pos_frac_per_var;
	delete[] var_graph_found;

	writeFeature("KLB-featuretime", gSW.TotalLap() - myTime);
	if (gSW.TotalLap() - myTime > TOTAL_TIMEOUT) {
		double dummy[] = {RESERVED_VALUE, RESERVED_VALUE};
		writeStats(dummy, 2, "CG");
		writeFeature("CG-entropy", RESERVED_VALUE);
		writeStats(dummy, 2, "cluster-coeff");
		writeFeature("cluster-coeff-entropy", RESERVED_VALUE);
		writeFeature("CG-featuretime", 0);
		return TOTAL_TIMEOUT_CODE;

	}

	myTime = gSW.TotalLap();
    if (doClauseGraphFeatures){
        if (DEB)
            p("c Clause graph...");
        clauseGraphFeatures(false);
    }
	if (DEB)
		p("c Done with base features");

	return FEAT_OK;
}

// compute clause graph features

void SATinstance::clauseGraphFeatures(bool realCC)
{
	int *clause;
	double *nodeDegrees = new double[numActiveClauses];
	int *nodeDegreesUnnormalized = new int[numActiveClauses];
	double *clusterCoeffs = new double[numActiveClauses];
	int count = 0;
	int inner_count = 0;
	//  set<int*> neighbors; // = new set<int *>;
	int nsize, confCount;
	set<int *>::iterator it, it2;
	myTime = gSW.TotalLap();
	cgTimeout = false;

	for (clause = firstClause(); clause != NULL; clause = nextClause())
	{
		set<int *> neighbors; // = new set<int *>;
		neighborClauses(clause, &neighbors);
		nsize = neighbors.size();
		nodeDegrees[count] = double(nsize) / double(numActiveClauses);
		nodeDegreesUnnormalized[count] = nsize;

		if (!realCC)
		{
			neighbors.insert(clause);
			nsize++;
		}

		//    fprintf(stderr, "Clause %d neighbours %d\n", count, nsize-1);

		confCount = 0;
		inner_count = 0;
		for (it = neighbors.begin(); it != neighbors.end(); it++)
		{

			// for(it2=neighbors->begin(); it2!=neighbors->end(); it2++) {
			for (it2 = it; it2 != neighbors.end(); it2++)
			{
				if (*it < *it2)
				{
					if (conflicted(*it, *it2))
					{
						confCount++;
					}
				}
			}
			printf("NEI%f, %d\n", gSW.TotalLap(), inner_count);

			if (inner_count++ % 100 == 0 && gSW.TotalLap() - myTime > TOTAL_TIMEOUT)
			{
				
				for (int i = 0; i < numActiveClauses; i++)
				{
					nodeDegreesUnnormalized[i] = -1;
					nodeDegrees[i] = -1;
					clusterCoeffs[i] = -1;
				}

				cgTimeout = true;
				break;
			}
		}

		if (cgTimeout)
			break;
		//    fprintf(stderr, "Clause %d conflicts %d\n", count, confCount);

		if (nsize > 1)
			clusterCoeffs[count] = 2.0 * (double)confCount / (double)(nsize * (nsize - 1));
		else
			clusterCoeffs[count] = 0;

		//    fprintf(stderr, "Clause %d coeff %lf\n", count, clusterCoeffs[count]);

		//    fprintf(stderr, "Time so far: %f\n", sw.Lap());
		if (count % 10000 == 0 && gSW.TotalLap() - myTime > TOTAL_TIMEOUT)
		{
			
			for (int i = 0; i < numActiveClauses; i++)
			{
				nodeDegreesUnnormalized[i] = -1;
				nodeDegrees[i] = -1;
				clusterCoeffs[i] = -1;
			}

			cgTimeout = true;
			break;
		}

		count++;
	}

	writeStats(nodeDegrees, count, "CG");
	if (!cgTimeout)
		writeFeature("CG-entropy", array_entropy(nodeDegreesUnnormalized, numActiveClauses, numActiveClauses + 1));
	else
		writeFeature("CG-entropy", RESERVED_VALUE);

	writeStats(clusterCoeffs, count, "cluster-coeff");
	writeFeature("cluster-coeff-entropy", array_entropy(clusterCoeffs, numActiveClauses, 100, 1));

	//  delete neighbors;
	delete[] nodeDegrees;
	delete[] nodeDegreesUnnormalized;
	delete[] clusterCoeffs;

	writeFeature("CG-featuretime", gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();
}

// compute LP
int SATinstance::compute_lp(bool doComp)
{

	if (!doComp)
	{
		// -- just write dummy values
		// == NOTE: These should be in the same order as before
		writeFeature("LP_OBJ", RESERVED_VALUE);
		double dummy[] = {RESERVED_VALUE, RESERVED_VALUE};
		writeStats(dummy, 2, "LPSLack");
		writeFeature("lpIntRatio", RESERVED_VALUE);
		writeFeature("lpTIME", 0);
		lp_return_val = LP_TIMEOUT;
		return lp_return_val;
	}

	REAL *solvedInstanceVarPtr;
	if (DEB)
		p("c Starting lp...");

	Stopwatch sw;

	sw.Start();

	int clauseCount = 0;
	for (int *clause = firstClause(); clause != NULL; clause = nextClause())
		clauseCount++;

	lprec *lp = make_lp(clauseCount, numActiveVars);
	if (lp == NULL)
	{
		fprintf(stderr, "c ERROR: compute_lp\n");
		return (lp_return_val = LP_ERROR);
	}

	/* Since we're dealing only with ACTIVE vars, this should
	be changed to index by the active var's position */

	REAL* obj = new REAL[numActiveVars + 1];
	REAL* cons = new REAL[numActiveVars + 1];

	// -- set up bounds
	for (int i = 1; i <= numActiveVars; i++)
	{
		obj[i] = 0;
		set_upbo(lp, i, 1);
		set_lowbo(lp, i, 0);
	}

	map<int, int> trans_for, trans_back;

	mkVarTranslation(&trans_for, &trans_back);

	int res;
	// -- add constraints
	// -- and compute objective
	for (int *clause = firstClause(); clause != NULL; clause = nextClause())
	{
		for (int i = 1; i <= numActiveVars; i++)
			cons[i] = 0;
		REAL rhs = 1;

		for (int lit = firstLitInClause(clause); lit != 0; lit = nextLitInClause())
		{
			int var = trans_for[ABS(lit)] + 1;

			if (negative(lit))
			{
				cons[var] = -1;
				obj[var]--;
				rhs--;
			}
			else
			{
				cons[var] = 1;
				obj[var]++;
			}
		}

		if (add_constraint(lp, cons, GE, rhs) != TRUE)
		{
			fprintf(stderr, "c ERROR: compute_lp add_constraint\n");
			delete_lp(lp);
			return (lp_return_val = LP_ERROR);
		}

		if (sw.Lap() > LP_TIME_LIMIT)
		{
			res = LP_TIMEOUT;
			break;
		}
	}

	// -- objective
	set_obj_fn(lp, obj);

	set_maxim(lp);

	// == print_lp(lp);

	set_timeout(lp, LP_TIME_LIMIT); // in seconds

	if (sw.Lap() > LP_TIME_LIMIT)
		res = LP_TIMEOUT;
	else
		res = solve(lp);

	if (res == INFEASIBLE)
	{
		//        delete_lp(lp);
		//  printf("s UNSATISFIABLE\n");
		lp_return_val = LP_UNSAT;
	}
	else if (res == TIMEOUT)
	{
		//        delete_lp(lp);
		lp_return_val = LP_TIMEOUT;
	}
	else if (res != OPTIMAL)
	{
		fprintf(stderr, "c ERROR in lp solve!\n");
		//        delete_lp(lp);
		lp_return_val = LP_ERROR;
	}

	// print_objective(lp);
	// print_solution(lp);

	REAL objval;
	if (res == OPTIMAL)
		objval = get_objective(lp) / numActiveClauses; // -- Divide by the number of clauses
	else
		objval = RESERVED_VALUE;
	writeFeature("LP_OBJ", objval);

	if (res == OPTIMAL)
		if (!get_ptr_variables(lp, &solvedInstanceVarPtr))
		{
			fprintf(stderr, "c ERROR in LP comp!\n");
			//        delete_lp(lp);
			lp_return_val = LP_ERROR;
		}

	// == NOTE: double check that indices are correct I.e start with zero

	int nInt = 0;
	for (int i = 0; i < numActiveVars; i++)
	{
		if (res == OPTIMAL)
		{
			cons[i] = (solvedInstanceVarPtr[i] < 1.0 - solvedInstanceVarPtr[i] ? solvedInstanceVarPtr[i] : 1.0 - solvedInstanceVarPtr[i]);
			if (cons[i] == 0)
				nInt++;
		}
		else
			cons[i] = RESERVED_VALUE;
	}

	writeStats(cons, numActiveVars, "LPSLack");

	double lpTime = gSW.TotalLap() - myTime;
	myTime = gSW.TotalLap();

	delete[] obj;
	delete[] cons;
	// if(res!=OPTIMAL)
	// lpTime=-1;

	writeFeature("lpIntRatio", double(nInt) / double(numActiveVars));
	writeFeature("lpTIME", lpTime);

	delete_lp(lp);

	return lp_return_val;
}

// compute diameter
int SATinstance::init_diameter(bool doComp)
{
	if (DEB)
		printf("c start diameter features ...\n");
	if (!doComp)
	{
		double dummy[] = {RESERVED_VALUE, RESERVED_VALUE};
		writeStats(dummy, 2, "DIAMETER");
		writeFeature("DIAMETER-entropy", RESERVED_VALUE);
		writeFeature("DIAMETER-featuretime", RESERVED_VALUE);
		return 0;
	}

	Stopwatch sw;
	sw.Start();
	int i, j, nrofliterals;
	int *var_occ;
	int max_diameter = 0, nr_max = 0, count_diameter = 0;
	int min_diameter = numVars, nr_min = 0;
	double sum_diameter = 0.0;
	int *diameterall;
	CCS = 0;

	dstack = new int[numVars];
	var_depth = new int[numVars + 1];
	var_occ = new int[numVars + 1];
	diameterall = new int[numVars];

	//	numActiveClausesWithVar[ABS(lits[litNum])]++;
	//    clauses[clauseNum][litNum] = lits[litNum];

	nrofliterals = 0;

	for (int i = 1; i <= numVars; i++)
	{
		var_occ[i] = numActiveClausesWithVar[i];
		nrofliterals = nrofliterals + numActiveClausesWithVar[i];
	}

	clause_stamps = new int[numActiveClauses];
	for (i = 0; i < numActiveClauses; i++)
		clause_stamps[i] = 0;

	dia_clause_table = new int[nrofliterals + numVars * 2];
	for (i = 0; i < nrofliterals + numVars * 2; i++)
		dia_clause_table[i] = -1;

	dia_clause_list = (int **)malloc(sizeof(int *) * (numVars + 1));
	int bef_nol = nrofliterals;
	nrofliterals = 0;
	for (i = 1; i <= numVars; i++)
	{
		dia_clause_list[i] = &dia_clause_table[nrofliterals];
		nrofliterals += var_occ[i] + 1;
		var_occ[i] = 0;
	}

	for (i = 0; i < numActiveClauses; i++)
		for (j = 0; j < clauseLengths[i]; j++)
			dia_clause_list[abs(clauses[i][j])][var_occ[abs(clauses[i][j])]++] = i;

	// free(var_occ);
	
	for (i = 1; i <= numVars; i++)
	{
		int diameter = computer_diameter(i, nrofliterals);
		if (sw.Lap() > DIA_TIME_LIMIT)
			break;
		//         printf("c diameter for %d: %d \n", i,diameter);
		if (diameter > max_diameter)
		{
			max_diameter = diameter;
			nr_max = 1;
		}
		else if (diameter == max_diameter)
			nr_max++;

		if (diameter > 0 && diameter < min_diameter)
		{
			min_diameter = diameter;
			nr_min = 1;
		}
		else if (diameter == min_diameter)
			nr_min++;

		if (diameter > 0)
		{
			sum_diameter += diameter;
			diameterall[count_diameter] = diameter;
			count_diameter++;
		}
	}

	//    p("Write diameter related Features...");

	// variable side in the bipartite graph
	writeStats(diameterall, count_diameter, "DIAMETER");
	writeFeature("DIAMETER-entropy", array_entropy(diameterall, count_diameter, max_diameter + 1));

	// printf("c diameter():: MIN: %i (#%i) MAX: %i (#%i) AVG: %.3f\n", min_diameter, nr_min, max_diameter, nr_max, sum_diameter / count_diameter );

	writeFeature("DIAMETER-featuretime", gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();

	delete[] clause_stamps;
	delete[] dia_clause_table;
	delete[] dstack;
	delete[] var_depth;
	delete[] var_occ;
	delete[] diameterall;
	return 1;
}

int SATinstance::computer_diameter(const int varnr, int nrofliterals)
{
	int i, _varnr, current_depth = 0;
	int *_dstackp = dstack;
	int *dstackp = dstack;
	int *myclauses, _myclause;

	if (dia_clause_list[varnr][0] == -1)
		return 0;

	for (i = 1; i < numVars; i++)
		var_depth[i] = 0;

	CCS++;

	DIAMETER_FIX(varnr);
	//	printf ("variable %d has depth of %d %d \n", varnr, *_dstackp, *(dstackp-1));
	for (int i=0; i< nrofliterals + numVars; i++){
		if (dia_clause_table[i] == 1462274){
			printf("asdasdf %d, %d\n", i, dia_clause_table[i]);
		}
	}

	while (_dstackp < dstackp)
	{
		_varnr = *(_dstackp++);
		current_depth = var_depth[_varnr];

		myclauses = dia_clause_list[_varnr];
		while (*myclauses != -1)
		{
			_myclause = *(myclauses++);
			if (clause_stamps[_myclause] == CCS)
				continue;
			clause_stamps[_myclause] = CCS;

			for (i = 0; i < clauseLengths[_myclause]; i++)
			{
				DIAMETER_FIX(abs(clauses[_myclause][i]));
			}
		}
	}
	/*
		for( i = 1; i <= current_depth; i++ )
		{
		int j;
		printf("\ndepth %i :: ", i);
		for( j = 1; j <= nrofvars; j++ )
			if( var_depth[j] == i )
			printf("%i ", j);
		}
	*/
	return current_depth;
}

int *SATinstance::firstClause()
{
	currentClause = 0;
	return nextClause();
}

int *SATinstance::nextClause()
{
	while (currentClause < numClauses && clauseStates[currentClause] != ACTIVE)
		currentClause++;

	return currentClause >= numClauses ? NULL : clauses[currentClause++];
}

int SATinstance::firstLitInClause(int *clause)
{
	currentClauseForLitIter = clause;
	currentLit = 0;
	return nextLitInClause();
}

int SATinstance::nextLitInClause()
{
	while (currentClauseForLitIter[currentLit] != 0 &&
		   varStates[ABS(currentClauseForLitIter[currentLit])] != UNASSIGNED)
		currentLit++;

	if (currentClauseForLitIter[currentLit] == 0)
		return 0;
	else
		return currentClauseForLitIter[currentLit++];
}

int SATinstance::firstLitInClause2(int *clause)
{
	currentClauseForLitIter2 = clause;
	currentLit2 = 0;
	return nextLitInClause2();
}

int SATinstance::nextLitInClause2()
{
	while (currentClauseForLitIter2[currentLit2] != 0 &&
		   varStates[ABS(currentClauseForLitIter2[currentLit2])] != UNASSIGNED)
		currentLit2++;

	if (currentClauseForLitIter2[currentLit2] == 0)
		return 0;
	else
		return currentClauseForLitIter2[currentLit2++];
}

int *SATinstance::firstClauseWithVar(int var, bool ispos)
{
	currentClauseWithVar = 0;
	currentVarForClauseIter = var;
	posClauses = ispos;
	return nextClauseWithVar();
}

int *SATinstance::nextClauseWithVar()
{
	if (posClauses)
	{
		while (currentClauseWithVar < (int)posClausesWithVar[currentVarForClauseIter].size() &&
			   clauseStates[posClausesWithVar[currentVarForClauseIter][currentClauseWithVar]] != ACTIVE)
			currentClauseWithVar++;
		if (currentClauseWithVar == (int)posClausesWithVar[currentVarForClauseIter].size())
			return NULL;
		else
			return clauses[posClausesWithVar[currentVarForClauseIter][currentClauseWithVar++]];
	}
	else
	{
		while (currentClauseWithVar < (int)negClausesWithVar[currentVarForClauseIter].size() &&
			   clauseStates[negClausesWithVar[currentVarForClauseIter][currentClauseWithVar]] != ACTIVE)
			currentClauseWithVar++;
		if (currentClauseWithVar == (int)negClausesWithVar[currentVarForClauseIter].size())
			return NULL;
		else
			return clauses[negClausesWithVar[currentVarForClauseIter][currentClauseWithVar++]];
	}
}

bool SATinstance::isVarInClause(int *clause, int var)
{
	// it may not be worth the trouble, but this could be binary search

	int i = 0;
	while (clause[i] != 0 && ABS(clause[i]) < var)
		i++;

	return (ABS(clause[i]) == var);
}

// mean.  Skips over values of "RESERVED_VALUE"
inline double SATinstance::mean(int *array, int num)
{
	int total = 0, t, reserved_hits = 0;
	for (t = 0; t < num; t++)
	{
		if (array[t] == (int)RESERVED_VALUE)
		{
			reserved_hits++;
			continue;
		}
		total += array[t];
	}

	if (reserved_hits == num)
		return 0;
	return (double)total / (double)(num - reserved_hits);
}
inline double SATinstance::mean(double *array, int num)
{
	double total = 0.0;
	int t, reserved_hits = 0;
	for (t = 0; t < num; t++)
	{

		if (array[t] == (double)RESERVED_VALUE)
		{
			reserved_hits++;
			continue;
		}
		total += array[t];
	}

	if (reserved_hits == num)
		return 0;
	return total / (double)(num - reserved_hits);
}

// standard deviation
inline double SATinstance::stdev(int *array, int num, double mean)
{
	double dtotal = 0.0;
	int reserved_hits = 0;
	for (int t = 0; t < num; t++)
	{
		if (array[t] == (int)RESERVED_VALUE)
		{
			reserved_hits++;
			continue;
		}
		dtotal += square(array[t] - mean);
	}
	if (reserved_hits == num)
		return 0;
	return sqrt(dtotal / (double)(num - reserved_hits));
}
inline double SATinstance::stdev(double *array, int num, double mean)
{
	double dtotal = 0.0;
	int reserved_hits = 0;
	for (int t = 0; t < num; t++)
	{
		if (array[t] == (double)RESERVED_VALUE)
		{
			reserved_hits++;
			continue;
		}
		dtotal += square(array[t] - mean);
	}
	if (reserved_hits == num)
		return 0;
	return sqrt(dtotal / (double)(num - reserved_hits));
}

// min
inline int SATinstance::array_min(int *array, int num)
{
	int m = (1 << 30);
	int reserved_hits = 0;
	for (int t = 0; t < num; t++)
	{
		if (array[t] == (int)RESERVED_VALUE)
		{
			reserved_hits++;
			continue;
		}
		m = (m < array[t] ? m : array[t]);
	}
	if (reserved_hits == num)
		return 0;
	return m;
}

inline double SATinstance::array_min(double *array, int num)
{
	double m = (1 << 30);
	int reserved_hits = 0;
	for (int t = 0; t < num; t++)
	{
		if (array[t] == (double)RESERVED_VALUE)
		{
			reserved_hits++;
			continue;
		}
		m = (m < array[t] ? m : array[t]);
	}
	if (reserved_hits == num)
		return 0;
	return m;
}

// max
inline int SATinstance::array_max(int *array, int num)
{
	int m = 0;
	int reserved_hits = 0;
	for (int t = 0; t < num; t++)
	{
		if (array[t] == (int)RESERVED_VALUE)
		{
			reserved_hits++;
			continue;
		}

		m = (m > array[t] ? m : array[t]);
	}
	if (reserved_hits == num)
		return 0;
	return m;
}

inline double SATinstance::array_max(double *array, int num)
{
	double m = 0;
	int reserved_hits = 0;
	for (int t = 0; t < num; t++)
	{
		if (array[t] == (double)RESERVED_VALUE)
		{
			reserved_hits++;
			continue;
		}
		m = (m > array[t] ? m : array[t]);
	}
	if (reserved_hits == num)
		return 0;
	return m;
}

// entropy
inline double SATinstance::array_entropy(double *array, int num, int vals, int maxval)
{
	int *p = new int[vals + 1];
	double entropy = 0.0, pval;
	int t, res = 0;
	int idx;

	// initialize
	for (t = 0; t <= vals; t++)
		p[t] = 0;

	// make the distribution
	for (t = 0; t < num; t++)
	{
		if (array[t] == (double)RESERVED_VALUE)
		{
			res++;
			continue;
		}
		idx = (int)floor(array[t] / ((double)maxval / (double)vals));
		//      if (idx > maxval) idx = maxval;
		if (idx > vals)
			idx = vals;
		if (idx < 0)
			idx = 0;
		p[idx]++;
	}

	// find the entropy
	for (t = 0; t <= vals; t++)
	{
		if (p[t])
		{
			pval = double(p[t]) / double(num - res);
			entropy += pval * log(pval);
		}
	}

	delete[] p;

	return -1.0 * entropy;
}
inline double SATinstance::array_entropy(int *array, int num, int vals)
{
	int *p = new int[vals];
	double entropy = 0.0, pval;
	int t, res = 0;

	// initialize
	for (t = 0; t < vals; t++)
		p[t] = 0;

	// make the distribution
	for (t = 0; t < num; t++)
	{
		if (array[t] == (int)RESERVED_VALUE)
		{
			res++;
			continue;
		}
#ifdef DEBUG
		if (array[t] < 0 || array[t] >= vals)
		{
			printf("c ERROR: bad array indexing in array_entropy!");
			throw std::runtime_error("Bad array indexing in array_entropy");
		}
#endif
		p[array[t]]++;
	}

	// find the entropy
	for (t = 0; t < vals; t++)
	{

		//      fprintf(stderr, "Bin %d/%d: %d\n", t, vals, p[t]);

		if (p[t])
		{
			pval = double(p[t]) / double(num - res);
			entropy += pval * log(pval);
		}
	}

	delete[] p;
	return -1.0 * entropy;
}

// write out node stats
// could these stats all be computed in one pass?

void SATinstance::writeStats(int *array, int num, const char *name)
{
	double m = mean(array, num);
	char buffer[100];
	sprintf(buffer, "%s-mean", name);
	writeFeature(buffer, m);
	sprintf(buffer, "%s-coeff-variation", name);
	double sd = stdev(array, num, m);
	double cv = (fabs(m) < EPSILON && sd < EPSILON ? 0 : sd / m);
	writeFeature(buffer, cv);
	sprintf(buffer, "%s-min", name);
	writeFeature(buffer, (double)array_min(array, num));
	sprintf(buffer, "%s-max", name);
	writeFeature(buffer, (double)array_max(array, num));
}

void SATinstance::writeStats(double *array, int num, const char *name)
{
	double m = mean(array, num);
	char buffer[100];
	sprintf(buffer, "%s-mean", name);
	writeFeature(buffer, m);
	sprintf(buffer, "%s-coeff-variation", name);
	double sd = stdev(array, num, m);
	double cv = (fabs(m) < EPSILON && sd < EPSILON ? 0 : sd / m);
	writeFeature(buffer, cv);
	sprintf(buffer, "%s-min", name);
	writeFeature(buffer, array_min(array, num));
	sprintf(buffer, "%s-max", name);
	writeFeature(buffer, array_max(array, num));
}

void SATinstance::writeStatsSTDEV(int *array, int num, const char *name)
{
	double m = mean(array, num);
	char buffer[100];
	sprintf(buffer, "%s-mean", name);
	writeFeature(buffer, m);
	sprintf(buffer, "%s-stdev", name);
	writeFeature(buffer, stdev(array, num, m));
	sprintf(buffer, "%s-min", name);
	writeFeature(buffer, (double)array_min(array, num));
	sprintf(buffer, "%s-max", name);
	writeFeature(buffer, (double)array_max(array, num));
}

void SATinstance::writeStatsSTDEV(double *array, int num, const char *name)
{
	double m = mean(array, num);
	char buffer[100];
	sprintf(buffer, "%s-mean", name);
	writeFeature(buffer, m);
	sprintf(buffer, "%s-stdev", name);
	writeFeature(buffer, stdev(array, num, m));
	sprintf(buffer, "%s-min", name);
	writeFeature(buffer, array_min(array, num));
	sprintf(buffer, "%s-max", name);
	writeFeature(buffer, array_max(array, num));
}

void SATinstance::writeStatsQ(double *array, int num, const char *name)
{
	// First sort this array
	vector<double> vc;
	for (int i = 0; i < num; i++)
	{
		vc.push_back(array[i]);
	}
	sort(vc.begin(), vc.end());
	double foonum = num;
	int q90 = (int)floor(foonum * 0.9);
	int q10 = (int)floor(foonum * 0.1);
	int q50 = (int)floor(foonum * 0.5);
	int q75 = (int)floor(foonum * 0.75);
	int q25 = (int)floor(foonum * 0.25);
	char buffer[100];
	sprintf(buffer, "%s-q90", name);
	writeFeature(buffer, vc[q90]);
	sprintf(buffer, "%s-q10", name);
	writeFeature(buffer, vc[q10]);
	sprintf(buffer, "%s-q75", name);
	writeFeature(buffer, vc[q75]);
	sprintf(buffer, "%s-q25", name);
	writeFeature(buffer, vc[q25]);
	sprintf(buffer, "%s-q50", name);
	writeFeature(buffer, vc[q50]);
	/*    sprintf(buffer,"%s-qr9010",name);
	if (vc[q10]>0)
	writeFeature(buffer,vc[q90]/vc[q10]);
	else
	writeFeature(buffer,-1);
	sprintf(buffer,"%s-qr7525",name);
	if (vc[q25]>0){
	writeFeature(buffer,vc[q75]/vc[q25]);
	}
	else
	writeFeature(buffer,-1);
	*/
}

void SATinstance::print()
{
	printf("p cnf %d %d\n", numActiveVars, numActiveClauses);
	for (int clause = 0; clause < numClauses; clause++)
	{
		if (clauseStates[clause] != ACTIVE)
			continue;
		//    printf("%d\t", clauseLengths[clause]);
		for (int lit = 0; clauses[clause][lit] != 0; lit++)
			if (varStates[ABS(clauses[clause][lit])] == UNASSIGNED)
				printf("%d\t", clauses[clause][lit]);
		printf("%d", clauseLengths[clause]);
		printf("\n");
	}
	printf("\n");
}

void SATinstance::testBackTrack()
{

	int dummy1, dummy2;
	unitprop(dummy1, dummy2);

	int *origBinClauseNums = new int[numVars + 1];
	for (int i = 1; i <= numVars; i++)
		origBinClauseNums[i] = numBinClausesWithVar[i];

	int origVars = numActiveVars;
	int origClauses = numActiveClauses;
	retardedSearch();
	for (int i = 0; i < numClauses; i++)
	{

		if (clauseStates[i] != ACTIVE)
			continue;

		int numLits = 0;

		for (int j = 0; clauses[i][j] != 0; j++)
			if (varStates[ABS(clauses[i][j])] == UNASSIGNED)
				numLits++;

		assert(numLits == clauseLengths[i]);
	}
	assert(numActiveVars == origVars);
	assert(numActiveClauses == origClauses);

	for (int i = 1; i <= numVars; i++)
		assert(numBinClausesWithVar[i] == origBinClauseNums[i]);

	if (stupidSearch())
	{
		outputAssignment();
	}
	else
		printf("c No assignment found.\n");

	delete[] origBinClauseNums;

	exit(0);
}

void SATinstance::testAPI()
{
	print();

	printf("\n");
	for (int var = 1; var <= numVars; var++)
	{
		printf("Var: %d true\n", var);
		for (int *clause = firstClauseWithVar(var, true); clause != NULL; clause = nextClauseWithVar())
		{
			for (int lit = firstLitInClause(clause); lit != 0; lit = nextLitInClause())
				printf("%d\t", lit);
			printf("0\n");
		}
		printf("\n");
	}
	for (int var = 1; var <= numVars; var++)
	{
		printf("Var: %d false\n", var);
		for (int *clause = firstClauseWithVar(var, false); clause != NULL; clause = nextClauseWithVar())
		{
			for (int lit = firstLitInClause(clause); lit != 0; lit = nextLitInClause())
				printf("%d\t", lit);
			printf("0\n");
		}
		printf("\n");
	}

	return;
}

// ------------------------------------------------------------------
// Compute a map that maps active variables to their ordinal number
// ------------------------------------------------------------------

void SATinstance::mkVarTranslation(map<int, int> *trans_for, map<int, int> *trans_back)
{
	int var, ord;

	for (ord = 0, var = 1; var <= numVars; var++)
		if (varStates[var] == UNASSIGNED) // -- active
		{
			(*trans_for)[var] = ord;
			(*trans_back)[ord] = var;
			ord++;
		}
}

void SATinstance::buildTranslatedActiveClauses(std::vector<std::vector<int> > &translatedClauses,
											   int &translatedVarCount,
											   std::vector<int> &literalOccurrences,
											   std::vector<int> &clauseLiteralOccurrences)
{
	map<int, int> transFor;
	map<int, int> transBack;
	mkVarTranslation(&transFor, &transBack);

	translatedVarCount = static_cast<int>(transFor.size());
	translatedClauses.clear();
	clauseLiteralOccurrences.clear();
	literalOccurrences.assign(2 * translatedVarCount + 1, 0);

	for (int clauseNum = 0; clauseNum < numClauses; clauseNum++)
	{
		if (clauseStates[clauseNum] != ACTIVE)
			continue;

		std::vector<int> translatedClause;
		for (int litNum = 0; litNum < clauseLengths[clauseNum]; litNum++)
		{
			const int literal = clauses[clauseNum][litNum];
			if (literal == 0 || varStates[ABS(literal)] != UNASSIGNED)
				continue;

			const map<int, int>::const_iterator translatedVar = transFor.find(ABS(literal));
			if (translatedVar == transFor.end())
				continue;

			const int translatedIndex = translatedVar->second + 1;
			const int translatedLiteral = positive(literal) ? translatedIndex : -translatedIndex;
			translatedClause.push_back(translatedLiteral);

			const int literalIndex = translatedLiteral > 0 ? translatedLiteral : translatedVarCount + (-translatedLiteral);
			literalOccurrences[literalIndex]++;
		}

		if (!translatedClause.empty())
		{
			translatedClauses.push_back(translatedClause);
			clauseLiteralOccurrences.push_back(static_cast<int>(translatedClause.size()));
		}
	}
}

int SATinstance::structureFeatures(bool doComp)
{
	if (!doComp)
	{
		writeFeature("vig_modularty", RESERVED_VALUE);
		writeFeature("vig_d_poly", RESERVED_VALUE);
		writeFeature("cvig_db_poly", RESERVED_VALUE);
		writeFeature("variable_alpha", RESERVED_VALUE);
		writeFeature("structure-featuretime", RESERVED_VALUE);
		return 0;
	}

	std::vector<std::vector<int> > clausesVec;
	std::vector<int> literalOccurrences;
	std::vector<int> clauseLiteralOccurrences;
	int translatedVarCount = 0;
	buildTranslatedActiveClauses(clausesVec, translatedVarCount, literalOccurrences, clauseLiteralOccurrences);

	if (translatedVarCount == 0 || clausesVec.empty())
	{
		writeFeature("vig_modularty", 0.0);
		writeFeature("vig_d_poly", 0.0);
		writeFeature("cvig_db_poly", 0.0);
		writeFeature("variable_alpha", 0.0);
		writeFeature("structure-featuretime", gSW.TotalLap() - myTime);
		myTime = gSW.TotalLap();
		return 0;
	}

	std::vector<int> variableCount(translatedVarCount + 1, 0);
	for (size_t clauseIdx = 0; clauseIdx < clausesVec.size(); ++clauseIdx)
		for (size_t litIdx = 0; litIdx < clausesVec[clauseIdx].size(); ++litIdx)
			variableCount[ABS(clausesVec[clauseIdx][litIdx])]++;

	std::vector<int> occurrenceHistogram(clausesVec.size() + 1, 0);
	for (int var = 1; var <= translatedVarCount; ++var)
		occurrenceHistogram[variableCount[var]]++;

	std::vector<std::pair<int, int> > countOccurrences;
	for (size_t count = 0; count < occurrenceHistogram.size(); ++count)
		if (occurrenceHistogram[count] != 0)
			countOccurrences.push_back(std::make_pair(static_cast<int>(count), occurrenceHistogram[count]));

	const double totalOccurrences = std::accumulate(occurrenceHistogram.begin(), occurrenceHistogram.end(), 0.0);
	std::vector<int> X;
	std::vector<double> Y(countOccurrences.size(), 0.0);
	std::vector<double> syLogX(countOccurrences.size(), 0.0);
	std::vector<double> syX(countOccurrences.size(), 0.0);
	X.reserve(countOccurrences.size());
	for (size_t i = 0; i < countOccurrences.size(); ++i)
		X.push_back(countOccurrences[i].first);

	for (int i = static_cast<int>(countOccurrences.size()) - 2; i >= 0; --i)
	{
		Y[i] = Y[i + 1] + (totalOccurrences > 0.0 ? countOccurrences[i].second / totalOccurrences : 0.0);
		if (X[i] > 0)
		{
			syLogX[i] = syLogX[i + 1] + (countOccurrences[i].second / totalOccurrences) * log(static_cast<double>(X[i]));
			syX[i] = syX[i + 1] + (countOccurrences[i].second / totalOccurrences) * static_cast<double>(X[i]);
		}
	}
	double alpha = 0.0;
	if (X.size() >= 4)
		alpha = mostLikelyAlpha(X, Y, syLogX, syX);

	std::vector<std::vector<int> > vigAdj = initUndirectedAdjacency(translatedVarCount);
	std::vector<std::map<int, double> > vigWeighted(translatedVarCount + 1);
	for (size_t clauseIdx = 0; clauseIdx < clausesVec.size(); ++clauseIdx)
	{
		const std::vector<int> &clause = clausesVec[clauseIdx];
		if (clause.size() < 2)
			continue;

		const double weight = 1.0 / static_cast<double>((clause.size() * (clause.size() - 1)) / 2);
		for (size_t i = 0; i < clause.size(); ++i)
		{
			for (size_t j = i + 1; j < clause.size(); ++j)
			{
				const int lhs = ABS(clause[i]);
				const int rhs = ABS(clause[j]);
				vigAdj[lhs].push_back(rhs);
				vigAdj[rhs].push_back(lhs);
				vigWeighted[lhs][rhs] += weight;
				vigWeighted[rhs][lhs] += weight;
			}
		}
	}
	finalizeUndirectedAdjacency(vigAdj);

	std::vector<std::vector<int> > cvigAdj = initUndirectedAdjacency(translatedVarCount + static_cast<int>(clausesVec.size()));
	for (size_t clauseIdx = 0; clauseIdx < clausesVec.size(); ++clauseIdx)
	{
		const int clauseNode = translatedVarCount + static_cast<int>(clauseIdx) + 1;
		for (size_t litIdx = 0; litIdx < clausesVec[clauseIdx].size(); ++litIdx)
		{
			const int varNode = ABS(clausesVec[clauseIdx][litIdx]);
			cvigAdj[varNode].push_back(clauseNode);
			cvigAdj[clauseNode].push_back(varNode);
		}
	}
	finalizeUndirectedAdjacency(cvigAdj);

	writeFeature("vig_modularty", approximateLouvainModularity(vigWeighted, translatedVarCount));
	writeFeature("vig_d_poly", fractalSlope(burningByNodeDegree(vigAdj, translatedVarCount)));
	writeFeature("cvig_db_poly", fractalSlope(burningByNodeDegree(cvigAdj, translatedVarCount + static_cast<int>(clausesVec.size()))));
	writeFeature("variable_alpha", alpha);
	writeFeature("structure-featuretime", gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();
	return FEAT_OK;
}

int SATinstance::newCnfGraphFeatures(bool doComp)
{
	static const char *kGraphPrefixes[] = {
		"v_nd_p_",
		"v_nd_n_",
		"c_nd_p_",
		"c_nd_n_",
		"vg_al_",
		"cg_al_",
		"rg_",
		"big_"};

	if (!doComp)
	{
		for (size_t i = 0; i < sizeof(kGraphPrefixes) / sizeof(kGraphPrefixes[0]); ++i)
			writeReservedMantheyGraphFeatures(this, kGraphPrefixes[i]);
		writeFeature("ncnf-graphs-featuretime", RESERVED_VALUE);
		return 0;
	}

	std::vector<std::vector<int> > clausesVec;
	std::vector<int> literalOccurrences;
	std::vector<int> clauseLiteralOccurrences;
	int translatedVarCount = 0;
	buildTranslatedActiveClauses(clausesVec, translatedVarCount, literalOccurrences, clauseLiteralOccurrences);

	std::vector<double> vPos(translatedVarCount, 0.0), vNeg(translatedVarCount, 0.0);
	std::vector<double> cPos(clausesVec.size(), 0.0), cNeg(clausesVec.size(), 0.0);
	for (size_t clauseIdx = 0; clauseIdx < clausesVec.size(); ++clauseIdx)
	{
		for (size_t litIdx = 0; litIdx < clausesVec[clauseIdx].size(); ++litIdx)
		{
			const int literal = clausesVec[clauseIdx][litIdx];
			if (literal > 0)
			{
				vPos[literal - 1] += 1.0;
				cPos[clauseIdx] += 1.0;
			}
			else
			{
				vNeg[-literal - 1] += 1.0;
				cNeg[clauseIdx] += 1.0;
			}
		}
	}
	const std::vector<double> zeroWeights(1, 0.0);
	writeMantheyGraphFeatures(this, "v_nd_p_", vPos, zeroWeights);
	writeMantheyGraphFeatures(this, "v_nd_n_", vNeg, zeroWeights);
	writeMantheyGraphFeatures(this, "c_nd_p_", cPos, zeroWeights);
	writeMantheyGraphFeatures(this, "c_nd_n_", cNeg, zeroWeights);

	std::vector<std::vector<int> > vgAdj = initUndirectedAdjacency(translatedVarCount);
	std::map<EdgeKey, double> vgWeights;
	for (size_t clauseIdx = 0; clauseIdx < clausesVec.size(); ++clauseIdx)
	{
		const std::vector<int> &clause = clausesVec[clauseIdx];
		for (size_t i = 0; i < clause.size(); ++i)
		{
			int k = 0;
			for (size_t j = i + 1; j < clause.size(); ++j)
			{
				const int lhs = ABS(clause[i]);
				const int rhs = ABS(clause[j]);
				vgAdj[lhs].push_back(rhs);
				vgAdj[rhs].push_back(lhs);
				k++;
				vgWeights[EdgeKey(lhs, rhs)] = pow(2.0, -k);
			}
		}
	}
	finalizeUndirectedAdjacency(vgAdj);
	writeMantheyGraphFeatures(this, "vg_al_", nodeDegreesFromUndirected(vgAdj, translatedVarCount), doubledUndirectedWeights(vgWeights));

	std::vector<std::vector<int> > cgAdj = initUndirectedAdjacency(static_cast<int>(clausesVec.size()));
	std::map<EdgeKey, double> cgWeights;
	for (size_t i = 0; i < clausesVec.size(); ++i)
	{
		const std::set<int> clauseLits(clausesVec[i].begin(), clausesVec[i].end());
		for (size_t j = i + 1; j < clausesVec.size(); ++j)
		{
			int shared = 0;
			for (std::set<int>::const_iterator it = clauseLits.begin(); it != clauseLits.end(); ++it)
				if (std::find(clausesVec[j].begin(), clausesVec[j].end(), *it) != clausesVec[j].end())
					shared++;
			if (shared > 0)
			{
				cgAdj[i + 1].push_back(j + 1);
				cgAdj[j + 1].push_back(i + 1);
				cgWeights[EdgeKey(i + 1, j + 1)] = static_cast<double>(shared);
			}
		}
	}
	finalizeUndirectedAdjacency(cgAdj);
	writeMantheyGraphFeatures(this, "cg_al_", nodeDegreesFromUndirected(cgAdj, static_cast<int>(clausesVec.size())), doubledUndirectedWeights(cgWeights));

	std::vector<std::vector<int> > rgAdj = initUndirectedAdjacency(static_cast<int>(clausesVec.size()));
	std::map<EdgeKey, double> rgWeights;
	for (size_t i = 0; i < clausesVec.size(); ++i)
	{
		const std::set<int> clauseA(clausesVec[i].begin(), clausesVec[i].end());
		for (size_t j = i + 1; j < clausesVec.size(); ++j)
		{
			const std::set<int> clauseB(clausesVec[j].begin(), clausesVec[j].end());
			int pivot = 0;
			int pivotCount = 0;
			for (std::set<int>::const_iterator it = clauseA.begin(); it != clauseA.end(); ++it)
			{
				if (clauseB.count(-(*it)) != 0)
				{
					pivot = *it;
					pivotCount++;
				}
			}
			if (pivotCount != 1)
				continue;

			bool tautologicalResolvent = false;
			std::set<int> resolvent;
			for (std::set<int>::const_iterator it = clauseA.begin(); it != clauseA.end(); ++it)
			{
				if (*it == pivot)
					continue;
				if (resolvent.count(-(*it)) != 0)
				{
					tautologicalResolvent = true;
					break;
				}
				resolvent.insert(*it);
			}
			if (tautologicalResolvent)
				continue;
			for (std::set<int>::const_iterator it = clauseB.begin(); it != clauseB.end(); ++it)
			{
				if (*it == -pivot)
					continue;
				if (resolvent.count(-(*it)) != 0)
				{
					tautologicalResolvent = true;
					break;
				}
				resolvent.insert(*it);
			}
			if (tautologicalResolvent)
				continue;

			rgAdj[i + 1].push_back(j + 1);
			rgAdj[j + 1].push_back(i + 1);
			const size_t unionSize = clauseA.size() + clauseB.size();
			rgWeights[EdgeKey(i + 1, j + 1)] = pow(2.0, -(static_cast<double>(unionSize) - 2.0));
		}
	}
	finalizeUndirectedAdjacency(rgAdj);
	writeMantheyGraphFeatures(this, "rg_", nodeDegreesFromUndirected(rgAdj, static_cast<int>(clausesVec.size())), doubledUndirectedWeights(rgWeights));

	const int literalNodeCount = 2 * translatedVarCount;
	std::vector<std::vector<int> > bigAdj(literalNodeCount + 1);
	std::map<DirectedEdgeKey, double> bigWeights;
	for (size_t clauseIdx = 0; clauseIdx < clausesVec.size(); ++clauseIdx)
	{
		if (clausesVec[clauseIdx].size() != 2)
			continue;
		const int a = clausesVec[clauseIdx][0];
		const int b = clausesVec[clauseIdx][1];
		const DirectedEdgeKey ab(literalNodeIndex(-a, translatedVarCount), literalNodeIndex(b, translatedVarCount));
		const DirectedEdgeKey ba(literalNodeIndex(-b, translatedVarCount), literalNodeIndex(a, translatedVarCount));
		bigAdj[ab.from].push_back(ab.to);
		bigAdj[ba.from].push_back(ba.to);
		bigWeights[ab] = 1.0;
		bigWeights[ba] = 1.0;
	}
	for (int node = 1; node <= literalNodeCount; ++node)
	{
		sort(bigAdj[node].begin(), bigAdj[node].end());
		bigAdj[node].erase(unique(bigAdj[node].begin(), bigAdj[node].end()), bigAdj[node].end());
	}
	writeMantheyGraphFeatures(this, "big_", directedOutDegrees(bigAdj, literalNodeCount), directedWeights(bigWeights));

	writeFeature("ncnf-graphs-featuretime", gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();
	return FEAT_OK;
}

int SATinstance::newCnfConstraintFeatures(bool doComp)
{
	static const char *kConstraintPrefixes[] = {"and_", "band_", "exo_"};
	if (!doComp)
	{
		for (size_t i = 0; i < sizeof(kConstraintPrefixes) / sizeof(kConstraintPrefixes[0]); ++i)
			writeReservedMantheyGraphFeatures(this, kConstraintPrefixes[i]);
		writeFeature("ncnf-constraints-featuretime", RESERVED_VALUE);
		return 0;
	}

	std::vector<std::vector<int> > clausesVec;
	std::vector<int> literalOccurrences;
	std::vector<int> clauseLiteralOccurrences;
	int translatedVarCount = 0;
	buildTranslatedActiveClauses(clausesVec, translatedVarCount, literalOccurrences, clauseLiteralOccurrences);

	const int literalNodeCount = 2 * translatedVarCount;
	std::vector<std::vector<int> > bigAdj(literalNodeCount + 1);
	for (size_t clauseIdx = 0; clauseIdx < clausesVec.size(); ++clauseIdx)
	{
		if (clausesVec[clauseIdx].size() != 2)
			continue;
		const int a = clausesVec[clauseIdx][0];
		const int b = clausesVec[clauseIdx][1];
		bigAdj[literalNodeIndex(-a, translatedVarCount)].push_back(literalNodeIndex(b, translatedVarCount));
		bigAdj[literalNodeIndex(-b, translatedVarCount)].push_back(literalNodeIndex(a, translatedVarCount));
	}
	for (int node = 1; node <= literalNodeCount; ++node)
	{
		sort(bigAdj[node].begin(), bigAdj[node].end());
		bigAdj[node].erase(unique(bigAdj[node].begin(), bigAdj[node].end()), bigAdj[node].end());
	}

	std::vector<std::vector<int> > andAdj = initUndirectedAdjacency(literalNodeCount);
	std::vector<std::vector<int> > bandAdj = initUndirectedAdjacency(literalNodeCount);
	std::vector<std::vector<int> > exoAdj = initUndirectedAdjacency(literalNodeCount);
	std::map<EdgeKey, double> andWeights;
	std::map<EdgeKey, double> bandWeights;
	std::map<EdgeKey, double> exoWeights;

	auto hasImplication = [&](int fromLiteral, int toLiteral) {
		const std::vector<int> &neighbors = bigAdj[literalNodeIndex(fromLiteral, translatedVarCount)];
		return binary_search(neighbors.begin(), neighbors.end(), literalNodeIndex(toLiteral, translatedVarCount));
	};

	for (size_t clauseIdx = 0; clauseIdx < clausesVec.size(); ++clauseIdx)
	{
		const std::vector<int> &clause = clausesVec[clauseIdx];
		if (clause.size() <= 2)
			continue;

		bool isExo = true;
		for (size_t i = 0; i < clause.size(); ++i)
		{
			for (size_t j = 0; j < clause.size(); ++j)
			{
				if (i == j)
					continue;
				if (!hasImplication(clause[i], -clause[j]))
				{
					isExo = false;
					break;
				}
			}
			if (!isExo)
				break;
		}

		const double weight = pow(2.0, -static_cast<double>(clause.size()));
		if (isExo)
		{
			for (size_t i = 0; i < clause.size(); ++i)
			{
				for (size_t j = 0; j < clause.size(); ++j)
				{
					if (i == j)
						continue;
					const int source = literalNodeIndex(clause[j], translatedVarCount);
					const int target = literalNodeIndex(-clause[i], translatedVarCount);
					andAdj[source].push_back(target);
					andAdj[target].push_back(source);
					andWeights[EdgeKey(source, target)] = weight;

					const int exoSource = literalNodeIndex(clause[j], translatedVarCount);
					const int exoTarget = literalNodeIndex(clause[i], translatedVarCount);
					exoAdj[exoSource].push_back(exoTarget);
					exoAdj[exoTarget].push_back(exoSource);
					exoWeights[EdgeKey(exoSource, exoTarget)] = 1.0;
				}
			}
			continue;
		}

		bool blocked = true;
		for (size_t i = 0; i < clause.size(); ++i)
		{
			const int literalIndex = clause[i] > 0 ? clause[i] : translatedVarCount + (-clause[i]);
			if (literalOccurrences[literalIndex] > 3)
			{
				blocked = false;
				break;
			}
		}
		if (!blocked)
			continue;

		for (size_t i = 0; i < clause.size(); ++i)
		{
			for (size_t j = 0; j < clause.size(); ++j)
			{
				if (i == j)
					continue;
				const int source = literalNodeIndex(clause[j], translatedVarCount);
				const int target = literalNodeIndex(-clause[i], translatedVarCount);
				bandAdj[source].push_back(target);
				bandAdj[target].push_back(source);
				bandWeights[EdgeKey(source, target)] = weight;
			}
		}
	}

	finalizeUndirectedAdjacency(andAdj);
	finalizeUndirectedAdjacency(bandAdj);
	finalizeUndirectedAdjacency(exoAdj);

	writeMantheyGraphFeatures(this, "and_", nodeDegreesFromUndirected(andAdj, literalNodeCount), doubledUndirectedWeights(andWeights));
	writeMantheyGraphFeatures(this, "band_", nodeDegreesFromUndirected(bandAdj, literalNodeCount), doubledUndirectedWeights(bandWeights));
	writeMantheyGraphFeatures(this, "exo_", nodeDegreesFromUndirected(exoAdj, literalNodeCount), doubledUndirectedWeights(exoWeights));

	writeFeature("ncnf-constraints-featuretime", gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();
	return FEAT_OK;
}

int SATinstance::newCnfRwhFeatures(bool doComp)
{
	if (!doComp)
	{
		for (int iteration = 0; iteration < 3; ++iteration)
		{
			char featureName[64];
			snprintf(featureName, sizeof(featureName), "rwh_%d_mean", iteration);
			writeFeature(featureName, RESERVED_VALUE);
			snprintf(featureName, sizeof(featureName), "rwh_%d_coeff", iteration);
			writeFeature(featureName, RESERVED_VALUE);
			snprintf(featureName, sizeof(featureName), "rwh_%d_min", iteration);
			writeFeature(featureName, RESERVED_VALUE);
			snprintf(featureName, sizeof(featureName), "rwh_%d_max", iteration);
			writeFeature(featureName, RESERVED_VALUE);
		}
		writeFeature("ncnf-rwh-featuretime", RESERVED_VALUE);
		return 0;
	}

	std::vector<std::vector<int> > clausesVec;
	std::vector<int> literalOccurrences;
	std::vector<int> clauseLiteralOccurrences;
	int translatedVarCount = 0;
	buildTranslatedActiveClauses(clausesVec, translatedVarCount, literalOccurrences, clauseLiteralOccurrences);

	const int maxClauseSize = 10;
	const double gamma = 5.0;
	const double maxDouble = 1e201;
	double mu = 1.0;
	std::vector<double> lastPos(translatedVarCount + 1, 1.0), lastNeg(translatedVarCount + 1, 1.0);
	std::vector<double> thisPos(translatedVarCount + 1, 0.0), thisNeg(translatedVarCount + 1, 0.0);

	for (int iteration = 0; iteration < 3; ++iteration)
	{
		for (size_t clauseIdx = 0; clauseIdx < clausesVec.size(); ++clauseIdx)
		{
			const std::vector<int> &clause = clausesVec[clauseIdx];
			const int clauseLength = static_cast<int>(clause.size());
			if (clauseLength <= 1)
				continue;

			const int exponent = maxClauseSize < clauseLength ? 0 : maxClauseSize - clauseLength;
			const double clauseConstant = pow(gamma, static_cast<double>(exponent)) / pow(mu, static_cast<double>(clauseLength - 1));

			double clauseValue = 1.0;
			bool foundZero = false;
			for (size_t litIdx = 0; litIdx < clause.size(); ++litIdx)
			{
				const int literal = clause[litIdx];
				const double compValue = literal < 0 ? lastPos[-literal] : lastNeg[literal];
				if (fabs(compValue) < EPSILON)
				{
					foundZero = true;
					break;
				}
				clauseValue *= compValue;
			}
			if (foundZero)
				continue;

			clauseValue *= clauseConstant;
			for (size_t litIdx = 0; litIdx < clause.size(); ++litIdx)
			{
				const int literal = clause[litIdx];
				const double compValue = literal < 0 ? lastPos[-literal] : lastNeg[literal];
				const double updateValue = clauseValue / compValue;
				if (literal > 0)
					thisPos[literal] += updateValue;
				else
					thisNeg[-literal] += updateValue;
			}
		}

		std::vector<double> sequence;
		sequence.reserve(2 * translatedVarCount);
		mu = 0.0;
		for (int var = 1; var <= translatedVarCount; ++var)
		{
			const double posValue = MIN(thisPos[var], maxDouble);
			const double negValue = MIN(thisNeg[var], maxDouble);
			sequence.push_back(posValue);
			sequence.push_back(negValue);
			mu += posValue + negValue;
		}
		if (translatedVarCount > 0)
			mu /= (2.0 * static_cast<double>(translatedVarCount));
		if (mu < 1.0)
			mu = 1.0;

		const double rwhMean = sequence.empty() ? 0.0 : std::accumulate(sequence.begin(), sequence.end(), 0.0) / static_cast<double>(sequence.size());
		double variance = 0.0;
		for (size_t i = 0; i < sequence.size(); ++i)
			variance += (sequence[i] - rwhMean) * (sequence[i] - rwhMean);
		const double rwhStd = sequence.empty() ? 0.0 : sqrt(variance / static_cast<double>(sequence.size()));
		const double rwhCoeff = (fabs(rwhMean) < EPSILON && fabs(rwhStd) < EPSILON ? 0.0 : rwhStd / rwhMean);

		char featureName[64];
		snprintf(featureName, sizeof(featureName), "rwh_%d_mean", iteration);
		writeFeature(featureName, rwhMean);
		snprintf(featureName, sizeof(featureName), "rwh_%d_coeff", iteration);
		writeFeature(featureName, rwhCoeff);
		snprintf(featureName, sizeof(featureName), "rwh_%d_min", iteration);
		writeFeature(featureName, sequence.empty() ? 0.0 : *min_element(sequence.begin(), sequence.end()));
		snprintf(featureName, sizeof(featureName), "rwh_%d_max", iteration);
		writeFeature(featureName, sequence.empty() ? 0.0 : *max_element(sequence.begin(), sequence.end()));

		lastPos.swap(thisPos);
		lastNeg.swap(thisNeg);
		std::fill(thisPos.begin(), thisPos.end(), 0.0);
		std::fill(thisNeg.begin(), thisNeg.end(), 0.0);
	}

	writeFeature("ncnf-rwh-featuretime", gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();
	return FEAT_OK;
}

void SATinstance::writeFeature(const char *name, double val)
{
	if (ignoreBadFeats)
		for (int i = 0; i < numBadFeats; i++)
			if (!strcmp(name, badFeatNames[i]))
			{
				// printf("c not including %s feature\n", name);
				return;
			}
	string s = string(string(name)); // featurePrefix)+ string(name);
	nameToIndex[s] = indexCount;
	featureNames[indexCount] = strdup(s.c_str());
	featureVals[indexCount] = val;
	indexCount++;
	if (indexCount >= MAX_FEATURES)
	{
		fprintf(stderr, "c TOO MANY FEATURES COMPUTED!\n");
		throw std::runtime_error("Too many features computed");
	}
	// printf("feature %d: %s; val: %f\n", indexCount, name, val);
}

void SATinstance::writeFeaturesToFile(const char *name)
{
	FILE *f = fopen(name, "a");
	if (f == NULL)
		writeFeaturesToFile(stderr);
	else
	{
		writeFeaturesToFile(f);
		fclose(f);
	}
}

void SATinstance::writeFeaturesToFile(FILE *f)
{
	if (f == NULL)
		f = stderr;
	for (int i = 0; i < indexCount - 1; i++)
	{
		fprintf(f, "%.9lf,", featureVals[i]);
	}
	fprintf(f, "%.9lf\n", featureVals[indexCount - 1]);
}

void SATinstance::writeFeatNamesToFile(const char *name)
{
	FILE *f = fopen(name, "a");
	if (f == NULL)
		writeFeatNamesToFile(stderr);
	else
	{
		writeFeatNamesToFile(f);
		fclose(f);
	}
}

void SATinstance::writeFeatNamesToFile(FILE *f)
{
	if (f == NULL)
		f = stderr;

	for (int i = 0; i < indexCount - 1; i++)
	{
		fprintf(f, "%s,", featureNames[i]);
	}
	fprintf(f, "%s\n", featureNames[indexCount - 1]);
}

void SATinstance::neighborClauses(int *clause, set<int *> *neighbors)
{
	static int lit, *otherClause, var;

	// neighbors->clear();
	// for each literal lit in clause
	for (lit = firstLitInClause(clause); lit != 0; lit = nextLitInClause())
	{
		// get all clauses that contain ~lit
		var = ABS(lit);
		for (otherClause = firstClauseWithVar(var, !positive(lit)); otherClause != NULL; otherClause = nextClauseWithVar())
		{
			neighbors->insert(otherClause);
		}
	}
}

bool SATinstance::conflicted(int *clause1, int *clause2)
{
	int lit1 = firstLitInClause(clause1);
	int lit2 = firstLitInClause2(clause2);

	// fprintf(stderr, "Checking conflicts\n");

	while (lit1 != 0 && lit2 != 0)
	{
		int v1 = ABS(lit1);
		int v2 = ABS(lit2);

		// fprintf(stderr, "Comparing %d vs %d\n", v1, v2);

		if ((v1 == v2) && (positive(lit1) != positive(lit2)))
			return true;

		else if (v1 < v2)
			lit1 = nextLitInClause();

		else
			lit2 = nextLitInClause2();
	}

	// fprintf(stderr, "No Conflict\n");
	return false;
}

#define NUM_VARS_TO_TRY 10
#define NUM_PROBES 5

int SATinstance::unitPropProbe(bool haltOnAssignment, bool doComp)
{

	// testBackTrack();
	if (DEB)
		p("c Unit prop probe...");
	if (!doComp)
	{

		int nextProbeDepth = 1;
		for (int j = 0; j < NUM_PROBES; j++)
		{
			nextProbeDepth *= 4;
			char featNameStr[100];
			sprintf(featNameStr, "vars-reduced-depth-%d", nextProbeDepth / 4);
			writeFeature(featNameStr, RESERVED_VALUE);
		}
		writeFeature("unit-featuretime", RESERVED_VALUE);
		return 0;
	}

	// NOTE: depth is number of vars manually set- not including unitprop.
	int currentDepth = 0;
	int origNumActiveVars = numActiveVars;
	bool reachedBottom = false;

	for (int probeNum = 0; probeNum < NUM_PROBES; probeNum++)
	{
		// this sets depth to 1, 4, 16, 64, 256
		int nextProbeDepth = 1;
		for (int j = 0; j < probeNum; j++)
			nextProbeDepth *= 4;

		while (currentDepth < nextProbeDepth && !reachedBottom)
		{
			int varsInMostBinClauses[NUM_VARS_TO_TRY];
			int numBin[NUM_VARS_TO_TRY];

			int arraySize = 0;
			for (int var = 1; var <= numVars; var++)
			{
				if (varStates[var] != UNASSIGNED)
					continue;
				if (arraySize < NUM_VARS_TO_TRY)
					arraySize++;

				int j = 0;
				while (j < arraySize - 1 && numBinClausesWithVar[var] < numBin[j])
					j++;
				for (int k = arraySize - 1; k > j; k--)
				{
					varsInMostBinClauses[k] = varsInMostBinClauses[k - 1];
					numBin[k] = numBin[k - 1];
				}
				varsInMostBinClauses[j] = var;
				numBin[j] = numBinClausesWithVar[var];
			}

			int maxPropsVar = 0;
			bool maxPropsVal;

			// if there are no binary clauses, just take the first unassigned var
			if (arraySize == 0)
			{
				maxPropsVar = 1;
				while (varStates[maxPropsVar] != UNASSIGNED && maxPropsVar < numVars)
					maxPropsVar++;
				maxPropsVal = true;
			}
			else
			{
				int maxProps = -1;

				for (int varNum = 0; varNum < arraySize; varNum++)
				{
					bool val = true;
					do
					{ // for val = true and val = false

						if (setVarAndProp(varsInMostBinClauses[varNum], val) &&
							numActiveVars <= 0)
						{
							if (haltOnAssignment)
							{
								outputAssignment();
								return DONE;
							}
						}

						int numProps = origNumActiveVars - numActiveVars - currentDepth;

						if (numProps > maxProps)
						{
							maxPropsVar = varsInMostBinClauses[varNum];
							maxPropsVal = val;
						}

						backtrack();

						val = !val;
					} while (val == false);
				}
			}

			assert(maxPropsVar != 0);

			if (!setVarAndProp(maxPropsVar, maxPropsVal))
				reachedBottom = true;

			else if (numActiveClauses == 0)
			{
				if (haltOnAssignment)
				{
					outputAssignment();
					return DONE;
				}
				reachedBottom = true;
			}

			currentDepth++;
		}

		char featNameStr[100];
		sprintf(featNameStr, "vars-reduced-depth-%d", nextProbeDepth);
		writeFeature(featNameStr, (double)(origNumActiveVars - numActiveVars - currentDepth) / numVars);
	}

	while (numActiveVars != origNumActiveVars)
		backtrack();
	writeFeature("unit-featuretime", gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();
	return FEAT_OK;
}

int SATinstance::sp(bool doComp)
{
	if (DEB)
		printf("c do survay propogation for %d second ...\n", SP_TIME_LIMIT);

	if (!doComp)
	{
		double dummy[] = {RESERVED_VALUE, RESERVED_VALUE};
		writeStats(dummy, 2, "SP-bias");
		writeStatsQ(dummy, 2, "SP-bias");
		writeStats(dummy, 2, "SP-unconstraint");
		writeStatsQ(dummy, 2, "SP-unconstraint");
		writeFeature("sp-featuretime", RESERVED_VALUE);
		return 0;
	}
	double **spresult;
	double *uncondfoo, *ratiofoo;
	int spsize;
	spresult = varsat::main(const_cast<char *>(inputFileName), spsize, SP_TIME_LIMIT);
	uncondfoo = new double[spsize];
	ratiofoo = new double[spsize];
	for (int i = 0; i < spsize; i++)
	{
		// cout << "c " << spresult[i][0] <<", "  << spresult[i][1] <<", "  << spresult[i][2] <<", " << endl;
		uncondfoo[i] = spresult[i][2];
		if (spresult[i][0] < 0.0000001 || spresult[i][1] < 0.0000001)
		{
			ratiofoo[i] = 1;
		}
		else
		{
			if (spresult[i][0] > spresult[i][1])
			{
				ratiofoo[i] = 1 - spresult[i][1] / spresult[i][0];
			}
			else
			{
				ratiofoo[i] = 1 - spresult[i][0] / spresult[i][1];
			}
		}
	}

	writeStats(ratiofoo, spsize, "SP-bias");
	writeStatsQ(ratiofoo, spsize, "SP-bias");
	writeStats(uncondfoo, spsize, "SP-unconstraint");
	writeStatsQ(uncondfoo, spsize, "SP-unconstraint");
	writeFeature("sp-featuretime", gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();

	for (int i = 0; i < spsize; i++)
		delete[] spresult[i];
	free(spresult);
	delete[] uncondfoo;
	delete[] ratiofoo;
	return 1;
}

namespace {
void writeReservedLocalSearchFeatures(SATinstance *sat, const char *prefix, const char *timeFeature)
{
	static const char *kFeatureSuffixes[] = {
		"BestSolution_Mean",
		"BestSolution_CoeffVariance",
		"FirstLocalMinStep_Mean",
		"FirstLocalMinStep_CoeffVariance",
		"FirstLocalMinStep_Median",
		"FirstLocalMinStep_Q.10",
		"FirstLocalMinStep_Q.90",
		"BestAvgImprovement_Mean",
		"BestAvgImprovement_CoeffVariance",
		"FirstLocalMinRatio_Mean",
		"FirstLocalMinRatio_CoeffVariance",
	};

	char featureName[128];
	for (const char *suffix : kFeatureSuffixes)
	{
		snprintf(featureName, sizeof(featureName), "%s_%s", prefix, suffix);
		sat->writeFeature(featureName, RESERVED_VALUE);
	}
	sat->writeFeature(timeFeature, RESERVED_VALUE);
}

bool parseLocalSearchOutput(SATinstance *sat, const char *outputFile, const char *prefix, int solvedCode)
{
	ifstream infile(outputFile);
	if (!infile)
		return false;

	char token[1024];
	char ignored[1024];
	char valueToken[1024];
	char featureName[1060];

	while (infile >> token)
	{
		const string key(token);
		const bool isFirst = key.find("First") == 0;
		const bool isBest = key.find("Best") == 0;
		const bool isSuccessfulRuns = key.find("SuccessfulRuns") == 0;
		if (!(isFirst || isBest || isSuccessfulRuns))
			continue;

		if (!(infile >> ignored >> valueToken))
			break;

		if (isFirst || isBest)
		{
			snprintf(featureName, sizeof(featureName), "%s_%s", prefix, token);
			sat->writeFeature(featureName, atof(valueToken));
		}
		else if (atoi(valueToken) == 1)
		{
			sat->solved = solvedCode;
		}
	}

	return true;
}

int runLocalSearchProbe(SATinstance *sat,
						const char *inputFile,
						bool doComp,
						BinSolver *solver,
						const char *prefix,
						const char *timeFeature,
						int solvedCode)
{
	if (!doComp)
	{
		writeReservedLocalSearchFeatures(sat, prefix, timeFeature);
		return 0;
	}

	char outputFile[512];
	if (!createTempPath(outputFile, sizeof(outputFile), "/featuresXXXXXX"))
	{
		fprintf(stderr, "c Error: could not create temporary output file for %s.\n", prefix);
		return 0;
	}

	solver->argv[9] = outputFile;
	solver->execute(inputFile, UBCSAT_TIME_LIMIT_INT);

	if (!parseLocalSearchOutput(sat, outputFile, prefix, solvedCode))
	{
		fprintf(stderr, "c Error: could not read from outputfile %s.\n", outputFile);
		solver->cleanup();
		throw std::runtime_error(std::string("Could not read local-search output file ") + outputFile);
	}

	sat->writeFeature(timeFeature, gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();
	remove(outputFile);
	solver->cleanup();
	return FEAT_OK;
}
}

int SATinstance::localSearchProbeGsat(const char *inputfile, bool doComp)
{
	if (DEB)
		p("c local search probe...");
	return runLocalSearchProbe(this, inputfile, doComp, SolverGsat, "gsat", "ls-gsat-featuretime", 3);
}

int SATinstance::localSearchProbeSaps(const char *inputfile, bool doComp)
{
	if (DEB)
		p("c local search probe...");
	return runLocalSearchProbe(this, inputfile, doComp, SolverSaps, "saps", "ls-saps-featuretime", 4);
}

#undef NUM_LS_PROBE

#define NUM_LOB_PROBE 30000

int SATinstance::lobjoisProbe(bool haltOnAssignment, bool doComp)
{
	if (DEB)
		p("c Lobjois probe...");
	if (!doComp)
	{
		writeFeature("lobjois-mean-depth-over-vars", RESERVED_VALUE);
		writeFeature("lobjois-log-num-nodes-over-vars", RESERVED_VALUE);
		writeFeature("lobjois-featuretime", RESERVED_VALUE);
		return 0;
	}

	Stopwatch sw;
	sw.Start();

	int depths[NUM_LOB_PROBE];

	int origNumActiveVars = numActiveVars;

	int probeNum = 0;

	while (probeNum < NUM_LOB_PROBE && sw.Lap() < LOBJOIS_TIME_LIMIT)
	{

		int var;
		bool val;

		do
		{
			if (numActiveVars == 0)
				if (numActiveClauses == 0 && haltOnAssignment)
				{
					outputAssignment();
					return DONE;
				}
				else
					break;

			var = 0;
			for (int stepsLeft = (rand() % numActiveVars); stepsLeft >= 0; stepsLeft--)
			{
				var++;
				while (varStates[var] != UNASSIGNED)
				{
					var++;
					if (var == numVars)
						var = 0;
				}
			}

			val = (rand() > RAND_MAX / 2);
			//      print();
			//      printf("Setting %d to %s\n", var, val ? "TRUE" : "FALSE");

		} while (setVarAndProp(var, val));

		//    printf("Reached bottom.\n\n");

		depths[probeNum] = origNumActiveVars - numActiveVars;

		while (numActiveVars != origNumActiveVars)
			backtrack();

		probeNum++;
	}

	writeFeature("lobjois-mean-depth-over-vars", mean(depths, probeNum) / numVars);
	int maxDepth = array_max(depths, probeNum);

	double sum = 0;
	for (int i = 0; i < probeNum; i++)
	{
		sum += pow(2.0, depths[i] - maxDepth);
	}

	double lobjois = maxDepth + log(sum / probeNum) / log(2.0);
	if (probeNum == 0)
	{
		lobjois = RESERVED_VALUE;
	}
	writeFeature("lobjois-log-num-nodes-over-vars", lobjois / numVars);
	writeFeature("lobjois-featuretime", gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();
	// printf("c # probes completed: %d\n", probeNum);
	return FEAT_OK;
}

void SATinstance::outputActiveFeat(bool *active)
{
	for (int i = 0; i < indexCount; i++)
		if (active[i])
			printf("%s\n", featureNames[i]);
}

int SATinstance::cl_prob(const char *outfile, bool doComp)
{
	int addClause[10000];
	double addClauseLength[10000];
	int mycl = 0;
	int mylit = 0;
	int returnVal;

	if (!doComp)
	{
		double dummy[] = {RESERVED_VALUE, RESERVED_VALUE};
		writeStats(dummy, 2, "cl-num");
		writeStatsQ(dummy, 2, "cl-num");
		writeStats(dummy, 2, "cl-size");
		writeStatsQ(dummy, 2, "cl-size");
		writeFeature("cl-featuretime", RESERVED_VALUE);

		return 0;
	}

	int mysize = 0;
	if (DEB)
		printf("c start clause learning features ...\n");
	returnVal = SolverZchaff->execute(outfile, CL_TIME_LIMIT);
	if (returnVal == 10 || returnVal == 20)
	{
		if (DEB)
			printf("c This instance is solved by CadiCal with %d!\n", returnVal);
		double dummy[] = {RESERVED_VALUE, RESERVED_VALUE};
		writeStats(dummy, 2, "cl-num");
		writeStatsQ(dummy, 2, "cl-num");
		writeStats(dummy, 2, "cl-size");
		writeStatsQ(dummy, 2, "cl-size");
		writeFeature("cl-featuretime", gSW.TotalLap() - myTime);
		solved = 2;
		return 0;
	}
	ifstream fin(SolverZchaff->outFileName);
	if (!fin)
	{
		fprintf(stderr, "c Error: could not read from %s.\n", outfile);
		throw std::runtime_error(std::string("Could not read clause-learning output from ") + outfile);
	}

	char strbuf[1024];
	while (fin.getline(strbuf, 512))
	{
		// printf("%s", strbuf);
		//   fin >> strbuf;
		string instream(strbuf);
		if (instream.find("c d ") == 0)
		{
			int a;
			int b;
			int c;
			char init1[10];
			char init2[10];
			istringstream ins(instream);
			ins >> init1 >> init2 >> a >> b >> c;
			if (a == 0)
			{
				mycl = c;
				mylit = b;
			}
			else
			{
				addClause[mysize] = c - mycl;
				addClauseLength[mysize] = double(b - mylit) / double(c - mycl);
				// printf("%d, %f\n", addClause[mysize], addClauseLength[mysize]);
				mysize = mysize + 1;
				mycl = c;
				mylit = b;
			}
		}
	}
	if (mysize == 0)
	{
		mysize = 1;
		addClause[0] = 0;
		addClauseLength[0] = 0;
	}
	// printf("c end read %d\n", mysize);
	double *clauseNum = new double[mysize];
	double *claueLen = new double[mysize];
	double previousLen = 0;
	for (int i = 0; i < mysize; i++)
	{
		clauseNum[i] = (double)addClause[i];
		claueLen[i] = addClauseLength[i];
		if (addClause[i] == 0)
		{
			claueLen[i] = previousLen;
		}
		else
		{
			previousLen = claueLen[i];
		}
		// printf("c length %f\n", claueLen[i]);
	}
	SolverZchaff->cleanup();
	writeStats(clauseNum, mysize, "cl-num");
	writeStatsQ(clauseNum, mysize, "cl-num");
	writeStats(claueLen, mysize, "cl-size");
	writeStatsQ(claueLen, mysize, "cl-size");
	writeFeature("cl-featuretime", gSW.TotalLap() - myTime);
	myTime = gSW.TotalLap();
	delete[] clauseNum;
	delete[] claueLen;
	return 0;
}

void SATinstance::finish_computation()
{
	writeFeature("solved", solved);
}
