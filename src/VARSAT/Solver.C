/****************************************************************************************[Solver.C]
VarSat -- Copyright (c) 2008, Eric Hsu.
Built upon MiniSat code by Niklas Een and Niklas Sorensson.
Their original copyright notice is repeated below.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include "Solver.h"
#include "Sort.h"
#include <cmath>

//////////////////////////////////////////////////////////////
// PARAMETERS GOVERNING SOLCOUNT HEURISTICS      EIH 14 JAN 08
//////////////////////////////////////////////////////////////
#define MAXITERATIONS 100
#define EPSILON 0.001

//=================================================================================================
// Constructor/Destructor:

namespace varsat
{
	Solver::Solver() :

					   // Parameters: (formerly in 'SearchParams')
					   var_decay(1 / 0.95), clause_decay(1 / 0.999), random_var_freq(0.02), restart_first(100), restart_inc(1.5), learntsize_factor((double)1 / (double)3), learntsize_inc(1.1)

					   // More parameters:
					   //
					   ,
					   expensive_ccmin(true), polarity_mode(polarity_false), verbosity(0)

					   // Statistics: (formerly in 'SolverStats')
					   //
					   ,
					   starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0), clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0), surveys_attempted(0), surveys_converged(0) // EIH 14 JAN 08
					   ,
					   surveys_time(0.0) // EIH 22 FEB 08

					   ,
					   heuristic_mode(heuristic_default) // EIH 12 JAN 08
					   ,
					   override_threshold(0.6) // EIH 13 FEB 08
					   ,
					   branching_mode(branching_default) // EIH 12 MAR 08
					   ,
					   timeout(43200) // EIH 10 APR 08
					   ,
					   cjm_mode(false) // EIH 17 MAR 08

					   ,
					   ok(true), cla_inc(1), var_inc(1), qhead(0), simpDB_assigns(-1), simpDB_props(0), order_heap(VarOrderLt(activity)), random_seed(1), progress_estimate(0), remove_satisfied(true)
	{
	}

	Solver::~Solver()
	{
		for (int i = 0; i < learnts.size(); i++)
			free(learnts[i]);
		for (int i = 0; i < clauses.size(); i++)
			free(clauses[i]);
	}

	//=================================================================================================
	// Minor methods:

	// Creates a new SAT variable in the solver. If 'decision_var' is cleared, variable will not be
	// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
	//
	Var Solver::newVar(bool sign, bool dvar)
	{
		int v = nVars();
		watches.push(); // (list for positive literal)
		watches.push(); // (list for negative literal)
		reason.push(NULL);
		assigns.push(toInt(l_Undef));
		level.push(-1);
		activity.push(0);
		seen.push(0);

		polarity.push((char)sign);
		decision_var.push((char)dvar);

		insertVarOrder(v);
		return v;
	}

	bool Solver::addClause(vec<Lit> &ps)
	{
		assert(decisionLevel() == 0);

		if (!ok)
			return false;
		else
		{
			// Check if clause is satisfied and remove false/duplicate literals:
			sort(ps);
			Lit p;
			int i, j;
			for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
				if (value(ps[i]) == l_True || ps[i] == ~p)
					return true;
				else if (value(ps[i]) != l_False && ps[i] != p)
					ps[j++] = p = ps[i];
			ps.shrink(i - j);
		}

		if (ps.size() == 0)
			return ok = false;
		else if (ps.size() == 1)
		{
			assert(value(ps[0]) == l_Undef);
			uncheckedEnqueue(ps[0]);
			return ok = (propagate() == NULL);
		}
		else
		{
			Clause *c = Clause_new(ps, false);
			clauses.push(c);
			attachClause(*c);

			// EIH 12 JAN 08 Note: Only original clauses get "added", while learnt (and original) clauses are
			// "attached."  So, by putting registerClause here, we are setting things up so that surveys only
			// consider clauses from the original theory.  We could also put it in "attachClause", in which
			// case future learnt clauses influence subsequent surveys.  In that case, we would also have to
			// manage the process of "detaching" such a learnt clause when we are purging the database.
			registerClause(*c);
		}

		return true;
	}

	// REGISTERCLAUSE: Register a clause as positively/negatively connected to its constituent variables.
	// EIH 12 JAN 08
	void Solver::registerClause(Clause &c)
	{

		// Loop through the literals in the clause.
		Lit p = lit_Undef;
		for (int i = 0; i < c.size(); i++)
		{
			p = c[i];
			// Register the current clause as a positive or negative clause for the literal's variable.
			sign(p) ? neg_clauses[var(p)].push(&c) : pos_clauses[var(p)].push(&c);
		}
	}

	void Solver::attachClause(Clause &c)
	{
		assert(c.size() > 1);
		watches[toInt(~c[0])].push(&c);
		watches[toInt(~c[1])].push(&c);
		if (c.learnt())
			learnts_literals += c.size();
		else
			clauses_literals += c.size();
	}

	void Solver::detachClause(Clause &c)
	{
		assert(c.size() > 1);
		assert(find(watches[toInt(~c[0])], &c));
		assert(find(watches[toInt(~c[1])], &c));
		remove(watches[toInt(~c[0])], &c);
		remove(watches[toInt(~c[1])], &c);

		// On killing clauses from the original theory, need to unregister them from the survey structures.  EIH 13 FEB 08
		if (!(c.learnt()))
		{
			Lit p = lit_Undef;
			// Loop through the clause's literals.
			for (int i = 0; i < c.size(); i++)
			{
				p = c[i];
				// Remove the clause positive or negative clauses list for the literal's variable.
				if (sign(p))
				{
					remove(neg_clauses[var(p)], &c);
				}
				else
				{
					remove(pos_clauses[var(p)], &c);
				}
			}
		}

		if (c.learnt())
			learnts_literals -= c.size();
		else
			clauses_literals -= c.size();
	}

	void Solver::removeClause(Clause &c)
	{
		detachClause(c);
		free(&c);
	}

	bool Solver::satisfied(const Clause &c) const
	{
		for (int i = 0; i < c.size(); i++)
			if (value(c[i]) == l_True)
				return true;
		return false;
	}

	// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
	//
	void Solver::cancelUntil(int level)
	{
		if (decisionLevel() > level)
		{
			for (int c = trail.size() - 1; c >= trail_lim[level]; c--)
			{
				Var x = var(trail[c]);
				assigns[x] = toInt(l_Undef);
				insertVarOrder(x);
			}
			qhead = trail_lim[level];
			trail.shrink(trail.size() - trail_lim[level]);
			trail_lim.shrink(trail_lim.size() - level);
		}
	}

	//=================================================================================================
	// Major methods:

	Lit Solver::pickBranchLit(int polarity_mode, double random_var_freq)
	{
		Lit retlit;

		// Dispatch one of the solution-counting heuristics if selected, unless we are in override mode.  EIH 12 JAN 08
		if ((heuristic_mode != heuristic_default) && (!heuristic_default_override))
		{
			return pickSolCountBranchLit();
		}

		Var next = var_Undef;

		// Random decision:
		if (drand(random_seed) < random_var_freq && !order_heap.empty())
		{
			next = order_heap[irand(random_seed, order_heap.size())];
			if (toLbool(assigns[next]) == l_Undef && decision_var[next])
				rnd_decisions++;
		}

		// Activity based decision:
		while (next == var_Undef || toLbool(assigns[next]) != l_Undef || !decision_var[next])
			if (order_heap.empty())
			{
				next = var_Undef;
				break;
			}
			else
				next = order_heap.removeMin();

		bool sign = false;
		switch (polarity_mode)
		{
		case polarity_true:
			sign = false;
			break;
		case polarity_false:
			sign = true;
			break;
		case polarity_user:
			sign = polarity[next];
			break;
		case polarity_rnd:
			sign = irand(random_seed, 2);
			break;
		default:
			assert(false);
		}

		retlit = (next == var_Undef ? lit_Undef : Lit(next, sign));

		if (verbosity >= 5)
		{
			reportf("Default Branching: ");
			printLit(retlit);
			reportf("\n");
		}

		return retlit;
	}

	/*_________________________________________________________________________________________________
	|
	|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
	|
	|  Description:
	|    Analyze conflict and produce a reason clause.
	|
	|    Pre-conditions:
	|      * 'out_learnt' is assumed to be cleared.
	|      * Current decision level must be greater than root level.
	|
	|    Post-conditions:
	|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
	|
	|  Effect:
	|    Will undo part of the trail, upto but not beyond the assumption of the current decision level.
	|________________________________________________________________________________________________@*/
	void Solver::analyze(Clause *confl, vec<Lit> &out_learnt, int &out_btlevel)
	{
		int pathC = 0;
		Lit p = lit_Undef;

		// Generate conflict clause:
		//
		out_learnt.push(); // (leave room for the asserting literal)
		int index = trail.size() - 1;
		out_btlevel = 0;

		do
		{
			assert(confl != NULL); // (otherwise should be UIP)
			Clause &c = *confl;

			if (c.learnt())
				claBumpActivity(c);

			for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++)
			{
				Lit q = c[j];

				if (!seen[var(q)] && level[var(q)] > 0)
				{
					varBumpActivity(var(q));
					seen[var(q)] = 1;
					if (level[var(q)] >= decisionLevel())
						pathC++;
					else
					{
						out_learnt.push(q);
						if (level[var(q)] > out_btlevel)
							out_btlevel = level[var(q)];
					}
				}
			}

			// Select next clause to look at:
			while (!seen[var(trail[index--])])
				;
			p = trail[index + 1];
			confl = reason[var(p)];
			seen[var(p)] = 0;
			pathC--;

		} while (pathC > 0);
		out_learnt[0] = ~p;

		// Simplify conflict clause:
		//
		int i, j;
		if (expensive_ccmin)
		{
			uint32_t abstract_level = 0;
			for (i = 1; i < out_learnt.size(); i++)
				abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

			out_learnt.copyTo(analyze_toclear);
			for (i = j = 1; i < out_learnt.size(); i++)
				if (reason[var(out_learnt[i])] == NULL || !litRedundant(out_learnt[i], abstract_level))
					out_learnt[j++] = out_learnt[i];
		}
		else
		{
			out_learnt.copyTo(analyze_toclear);
			for (i = j = 1; i < out_learnt.size(); i++)
			{
				Clause &c = *reason[var(out_learnt[i])];
				for (int k = 1; k < c.size(); k++)
					if (!seen[var(c[k])] && level[var(c[k])] > 0)
					{
						out_learnt[j++] = out_learnt[i];
						break;
					}
			}
		}
		max_literals += out_learnt.size();
		out_learnt.shrink(i - j);
		tot_literals += out_learnt.size();

		// Find correct backtrack level:
		//
		if (out_learnt.size() == 1)
			out_btlevel = 0;
		else
		{
			int max_i = 1;
			for (int i = 2; i < out_learnt.size(); i++)
				if (level[var(out_learnt[i])] > level[var(out_learnt[max_i])])
					max_i = i;
			Lit p = out_learnt[max_i];
			out_learnt[max_i] = out_learnt[1];
			out_learnt[1] = p;
			out_btlevel = level[var(p)];
		}

		for (int j = 0; j < analyze_toclear.size(); j++)
			seen[var(analyze_toclear[j])] = 0; // ('seen[]' is now cleared)
	}

	// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
	// visiting literals at levels that cannot be removed later.
	bool Solver::litRedundant(Lit p, uint32_t abstract_levels)
	{
		analyze_stack.clear();
		analyze_stack.push(p);
		int top = analyze_toclear.size();
		while (analyze_stack.size() > 0)
		{
			assert(reason[var(analyze_stack.last())] != NULL);
			Clause &c = *reason[var(analyze_stack.last())];
			analyze_stack.pop();

			for (int i = 1; i < c.size(); i++)
			{
				Lit p = c[i];
				if (!seen[var(p)] && level[var(p)] > 0)
				{
					if (reason[var(p)] != NULL && (abstractLevel(var(p)) & abstract_levels) != 0)
					{
						seen[var(p)] = 1;
						analyze_stack.push(p);
						analyze_toclear.push(p);
					}
					else
					{
						for (int j = top; j < analyze_toclear.size(); j++)
							seen[var(analyze_toclear[j])] = 0;
						analyze_toclear.shrink(analyze_toclear.size() - top);
						return false;
					}
				}
			}
		}

		return true;
	}

	/*_________________________________________________________________________________________________
	|
	|  analyzeFinal : (p : Lit)  ->  [void]
	|
	|  Description:
	|    Specialized analysis procedure to express the final conflict in terms of assumptions.
	|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
	|    stores the result in 'out_conflict'.
	|________________________________________________________________________________________________@*/
	void Solver::analyzeFinal(Lit p, vec<Lit> &out_conflict)
	{
		out_conflict.clear();
		out_conflict.push(p);

		if (decisionLevel() == 0)
			return;

		seen[var(p)] = 1;

		for (int i = trail.size() - 1; i >= trail_lim[0]; i--)
		{
			Var x = var(trail[i]);
			if (seen[x])
			{
				if (reason[x] == NULL)
				{
					assert(level[x] > 0);
					out_conflict.push(~trail[i]);
				}
				else
				{
					Clause &c = *reason[x];
					for (int j = 1; j < c.size(); j++)
						if (level[var(c[j])] > 0)
							seen[var(c[j])] = 1;
				}
				seen[x] = 0;
			}
		}

		seen[var(p)] = 0;
	}

	void Solver::uncheckedEnqueue(Lit p, Clause *from)
	{
		assert(value(p) == l_Undef);
		assigns[var(p)] = toInt(lbool(!sign(p))); // <<== abstract but not uttermost effecient
		level[var(p)] = decisionLevel();
		reason[var(p)] = from;
		trail.push(p);
	}

	/*_________________________________________________________________________________________________
	|
	|  propagate : [void]  ->  [Clause*]
	|
	|  Description:
	|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
	|    otherwise NULL.
	|
	|    Post-conditions:
	|      * the propagation queue is empty, even if there was a conflict.
	|________________________________________________________________________________________________@*/
	Clause *Solver::propagate()
	{
		Clause *confl = NULL;
		int num_props = 0;

		while (qhead < trail.size())
		{
			Lit p = trail[qhead++]; // 'p' is enqueued fact to propagate.
			vec<Clause *> &ws = watches[toInt(p)];
			Clause **i, **j, **end;
			num_props++;

			for (i = j = (Clause **)ws, end = i + ws.size(); i != end;)
			{
				Clause &c = **i++;

				// Make sure the false literal is data[1]:
				Lit false_lit = ~p;
				if (c[0] == false_lit)
					c[0] = c[1], c[1] = false_lit;

				assert(c[1] == false_lit);

				// If 0th watch is true, then clause is already satisfied.
				Lit first = c[0];
				if (value(first) == l_True)
				{
					*j++ = &c;
				}
				else
				{
					// Look for new watch:
					for (int k = 2; k < c.size(); k++)
						if (value(c[k]) != l_False)
						{
							c[1] = c[k];
							c[k] = false_lit;
							watches[toInt(~c[1])].push(&c);
							goto FoundWatch;
						}

					// Did not find watch -- clause is unit under assignment:
					*j++ = &c;
					if (value(first) == l_False)
					{
						confl = &c;
						qhead = trail.size();
						// Copy the remaining watches:
						while (i < end)
							*j++ = *i++;
					}
					else
						uncheckedEnqueue(first, &c);
				}
			FoundWatch:;
			}
			ws.shrink(i - j);
		}
		propagations += num_props;
		simpDB_props -= num_props;

		return confl;
	}

	/*_________________________________________________________________________________________________
	|
	|  reduceDB : ()  ->  [void]
	|
	|  Description:
	|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
	|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
	|________________________________________________________________________________________________@*/
	// struct reduceDB_lt { bool operator () (Clause* x, Clause* y) { return x->size() > 2 && (y->size() == 2 || x->activity() < y->activity()); } };

	void Solver::reduceDB()
	{
		int i, j;
		double extra_lim = cla_inc / learnts.size(); // Remove any clause below this activity

		//    sort(learnts, reduceDB_lt());
		for (i = j = 0; i < learnts.size() / 2; i++)
		{
			if (learnts[i]->size() > 2 && !locked(*learnts[i]))
				removeClause(*learnts[i]);
			else
				learnts[j++] = learnts[i];
		}
		for (; i < learnts.size(); i++)
		{
			if (learnts[i]->size() > 2 && !locked(*learnts[i]) && learnts[i]->activity() < extra_lim)
				removeClause(*learnts[i]);
			else
				learnts[j++] = learnts[i];
		}
		learnts.shrink(i - j);
	}

	void Solver::removeSatisfied(vec<Clause *> &cs)
	{
		int i, j;
		for (i = j = 0; i < cs.size(); i++)
		{
			if (satisfied(*cs[i]))
				removeClause(*cs[i]);
			else
				cs[j++] = cs[i];
		}
		cs.shrink(i - j);
	}

	/*_________________________________________________________________________________________________
	|
	|  simplify : [void]  ->  [bool]
	|
	|  Description:
	|    Simplify the clause database according to the current top-level assigment. Currently, the only
	|    thing done here is the removal of satisfied clauses, but more things can be put here.
	|________________________________________________________________________________________________@*/
	bool Solver::simplify()
	{
		assert(decisionLevel() == 0);

		if (!ok || propagate() != NULL)
			return ok = false;

		if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
			return true;

		// Remove satisfied clauses:
		removeSatisfied(learnts);
		if (remove_satisfied) // Can be turned off.
			removeSatisfied(clauses);

		// Remove fixed variables from the variable heap:
		order_heap.filter(VarFilter(*this));

		simpDB_assigns = nAssigns();
		simpDB_props = clauses_literals + learnts_literals; // (shouldn't depend on stats really, but it will do for now)

		return true;
	}

	/*_________________________________________________________________________________________________
	|
	|  search : (nof_conflicts : int) (nof_learnts : int) (params : const SearchParams&)  ->  [lbool]
	|
	|  Description:
	|    Search for a model the specified number of conflicts, keeping the number of learnt clauses
	|    below the provided limit. NOTE! Use negative value for 'nof_conflicts' or 'nof_learnts' to
	|    indicate infinity.
	|
	|  Output:
	|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
	|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
	|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
	|________________________________________________________________________________________________@*/
	lbool Solver::search(int nof_conflicts, int nof_learnts)
	{
		assert(ok);
		int backtrack_level;
		int conflictC = 0;
		vec<Lit> learnt_clause;

		starts++;

		// For solution counting heuristics.  EIH 14 JAN 08
		initializeBiases();
		// if (verbosity >= 4) reportBiases();

		// If we are restarting, we need to make sure to allow surveys again, even if we had turned them
		// off beacuse the biases started to fall below the threshold.   EIH 13 FEB 08
		heuristic_default_override = false;

		bool first = true;

		for (;;)
		{
			Clause *confl = propagate();
			if (confl != NULL)
			{
				// CONFLICT
				conflicts++;
				conflictC++;
				if (decisionLevel() == 0)
					return l_False;

				first = false;

				learnt_clause.clear();
				analyze(confl, learnt_clause, backtrack_level);
				cancelUntil(backtrack_level);
				assert(value(learnt_clause[0]) == l_Undef);

				if (learnt_clause.size() == 1)
				{
					uncheckedEnqueue(learnt_clause[0]);
				}
				else
				{
					Clause *c = Clause_new(learnt_clause, true);
					learnts.push(c);
					attachClause(*c);
					claBumpActivity(*c);
					uncheckedEnqueue(learnt_clause[0], c);
				}

				varDecayActivity();
				claDecayActivity();
			}
			else
			{
				// NO CONFLICT

				if (nof_conflicts >= 0 && conflictC >= nof_conflicts)
				{
					// Reached bound on number of conflicts:
					progress_estimate = progressEstimate();
					cancelUntil(0);
					return l_Undef;
				}

				// Simplify the set of problem clauses:
				if (decisionLevel() == 0 && !simplify())
					return l_False;

				if (nof_learnts >= 0 && learnts.size() - nAssigns() >= nof_learnts)
					// Reduce the set of learnt clauses:
					reduceDB();

				Lit next = lit_Undef;
				while (decisionLevel() < assumptions.size())
				{
					// Perform user provided assumption:
					Lit p = assumptions[decisionLevel()];
					if (value(p) == l_True)
					{
						// Dummy decision level:
						newDecisionLevel();
					}
					else if (value(p) == l_False)
					{
						analyzeFinal(~p, conflict);
						return l_False;
					}
					else
					{
						next = p;
						break;
					}
				}

				if (next == lit_Undef)
				{
					// New variable decision:
					decisions++;
					next = pickBranchLit(polarity_mode, random_var_freq);
					if (next == lit_Undef)
						// Model found:
						return l_True;
				}

				// Increase decision level and enqueue 'next'
				assert(value(next) == l_Undef);
				newDecisionLevel();
				uncheckedEnqueue(next);
			}
		}
	}

	double Solver::progressEstimate() const
	{
		double progress = 0;
		double F = 1.0 / nVars();

		for (int i = 0; i <= decisionLevel(); i++)
		{
			int beg = i == 0 ? 0 : trail_lim[i - 1];
			int end = i == decisionLevel() ? trail.size() : trail_lim[i];
			progress += pow(F, i) * (end - beg);
		}

		return progress / nVars();
	}

	bool Solver::solve(const vec<Lit> &assumps)
	{
		model.clear();
		conflict.clear();

		if (!ok)
			return false;

		assumps.copyTo(assumptions);

		double nof_conflicts = restart_first;
		double nof_learnts = nClauses() * learntsize_factor;
		lbool status = l_Undef;

		if (verbosity >= 1)
		{
			reportf("============================[ Search Statistics ]==============================\n");
			reportf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
			reportf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
			reportf("===============================================================================\n");
		}

		// Search:
		while (status == l_Undef)
		{
			if (verbosity >= 1)
				reportf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n", (int)conflicts, order_heap.size(), nClauses(), (int)clauses_literals, (int)nof_learnts, nLearnts(), (double)learnts_literals / nLearnts(), progress_estimate * 100), fflush(stdout);
			status = search((int)nof_conflicts, (int)nof_learnts);
			nof_conflicts *= restart_inc;
			nof_learnts *= learntsize_inc;

			// EIH 27 MAR 08 TIMEOUTS
			if (getCpuTime() > timeout)
			{
				reportf("\nTIMEOUT after %d seconds.\n", timeout);
				exit(-1);
			}
		}

		if (verbosity >= 1)
			reportf("===============================================================================\n");

		if (status == l_True)
		{
			// Extend & copy model:
			model.growTo(nVars());
			for (int i = 0; i < nVars(); i++)
				model[i] = value(i);
#ifndef NDEBUG
			verifyModel();
#endif
		}
		else
		{
			assert(status == l_False);
			if (conflict.size() == 0)
				ok = false;
		}

		cancelUntil(0);
		return status == l_True;
	}

	//=================================================================================================
	// Debug methods:

	void Solver::verifyModel()
	{
		bool failed = false;
		for (int i = 0; i < clauses.size(); i++)
		{
			assert(clauses[i]->mark() == 0);
			Clause &c = *clauses[i];
			for (int j = 0; j < c.size(); j++)
				if (modelValue(c[j]) == l_True)
					goto next;

			reportf("unsatisfied clause: ");
			printClause(*clauses[i]);
			reportf("\n");
			failed = true;
		next:;
		}

		assert(!failed);

		reportf("Verified %d original clauses.\n", clauses.size());
	}

	void Solver::checkLiteralCount()
	{
		// Check that sizes are calculated correctly:
		int cnt = 0;
		for (int i = 0; i < clauses.size(); i++)
			if (clauses[i]->mark() == 0)
				cnt += clauses[i]->size();

		if ((int)clauses_literals != cnt)
		{
			fprintf(stderr, "literal count: %d, real value = %d\n", (int)clauses_literals, cnt);
			assert((int)clauses_literals == cnt);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////// MAIN AREA FOR IMPLEMENTING VARSAT STUFF ///////////////////////////////////
	///////////////////////////////////////// 14 JAN 08 //////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	/*_________________________________________________________________________________________________
	|
	|  pickSolCountBranchLit : () -> Lit
	|
	|  Description:
	|    As with regular pickBranchLit, chooses a literal to branch on next, but using one of our
	|    solution counting techniques instead of activity/random.  Will do so according to the
	|	 heuristic_mode for the solver (BP, EMBP, SP, or EMSP).  Also, we can use different scoring
	|    methods to choose the next literal upon consulting the survey computed by these methods.
	|
	|  Helpers:
	|	 void computeSurvey()
	|
	|  EIH 12 JAN 08
	|________________________________________________________________________________________________@*/
	Lit Solver::pickSolCountBranchLit()
	{
		Var v = 0;
		Lit retlit = lit_Undef;
		Var next = var_Undef;
		bool sign = false;

		double maxscore = 0.0;

		// First compute a survey.
		//
		bool survey_converged;
		if (heuristic_mode == heuristic_LC)
		{
			survey_converged = computeLCSurvey();
		}
		else
		{
			survey_converged = computeSurvey();
		}

		// Handle non-convergence.  Here we just use whatever values the survey left off with.
		// Other options would be to use the regular pickSolCountBranchLit, or make to pick randomly.
		//
		if (!survey_converged)
		{
			if (verbosity >= 4)
			{
				reportf("Survey did not converge, proceeding with intermediate values.\n");
			}
		}

		// Go through the unfixed variables and find the one with the highest score according to chosen metric.
		// FUTURE: integrate other scoring functions besides scoreVarSkew.
		//
		for (v = 0; v < nVars(); v++)
		{
			// Skip any variables that are already fixed.
			if (!(toLbool(assigns[v]) == l_Undef && decision_var[v]))
				continue;

			if (scoreVarSkew(v) >= maxscore)
			{
				next = v;
				maxscore = scoreVarSkew(v);
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////
		// Set polarity of the chosen variable according to the chosen rule.
		// NOTE: setting sign to true means we want the literal to be false, i.e. sign means "bar"!
		/////////////////////////////////////////////////////////////////////////////////////////////

		// SUCCEED FIRST
		if ((branching_mode == branching_solution) ||
			(branching_mode == branching_default))
		{
			sign = (pos_biases[next] > neg_biases[next]) ? false : true;

			// FAIL FIRST
		}
		else if (branching_mode == branching_conflict)
		{
			sign = (pos_biases[next] > neg_biases[next]) ? true : false;

			// STREAMLINE (fail first, then succeed)
		}
		else if (branching_mode == branching_streamline)
		{
			if (trail.size() / nVars() < 0.1)
			{																 // if we are fixing the first 10% of variables
				sign = (pos_biases[next] > neg_biases[next]) ? true : false; // fail first
			}
			else
			{
				sign = (pos_biases[next] > neg_biases[next]) ? false : true; // succeed first
			}

			// MIXED (random 50/50 between fail and succeed)
		}
		else if (branching_mode == branching_mixed)
		{
			sign = irand(random_seed, 2);

			// undefined branching rule
		}
		else
		{
			assert(false);
		}

		if (pos_biases[next] >= 0.95)
			sign = false; // don't fail first on highly likely vars
		if (neg_biases[next] >= 0.95)
			sign = true;

		retlit = (next == var_Undef ? lit_Undef : Lit(next, sign));

		if (verbosity >= 3)
		{
			reportf("SolCount Branching: ");
			printLit(retlit);
			reportf("\n");
		}

		// Check whether we should turn off the surveys.
		// FUTURE: try other rules, like log-likelihood, fixed percentage, etc.
		//
		// if (verbosity >=3) {reportf("Skew: %f %f %f\n", scoreVarSkew(next), pos_biases[next], neg_biases[next]);}
		if (maxscore < override_threshold)
		{
			heuristic_default_override = true;
		}

		return retlit;
	}

	/*_________________________________________________________________________________________________
	|
	|  computeSurvey : () -> bool
	|
	|  Description:
	|    Starting with the current setting of the variables' bias arrays, computes a survey by
	|    repeatedly going through all the variables, updating them according to the rule specified
	|    in Solver::heuristic_mode.  Each pass through all the variables is called an "iteration"
	|    and each calculation for a single variable is called an "update".
	|
	|    Before using this method for the first time, the bias arrays should be initialized with
	|    Solver::initializeBiases.  From then on, you can reinitialize every time you start to
	|    compute a new survey, or you can leave the arrays as they are in order to start off from
	|    same point where the last survey left off.  This latter approach might have advantages in
	|    terms of speed and in terms of making sure we are "working within the same cluster" across
	|    surveys.
	|
	|  Effects:
	|    Updates a solver's pos_biases, neg_biases, and star_biases.
	|
	|  Returns:
	|    true for successful convergence
	|    false if process times out without converging
	|
	|  Helpers:
	|    double updateBP(Var v)
	|    double updateEMBPL(Var v)
	|    double updateEMBPG(Var v)
	|    double updateEMBPGV2(Var v)
	|    double updateSP(Var v)
	|    double updateEMSPL(Var v)
	|    double updateEMSPG(Var v)
	|    double updateEMSPGV2(Var v)
	|    double updateCC(Var v)
	|    double calculateFreedom(Var v, vec<Clause*>& clist, int& count)        i.e. \prod_clist (1 - \sigma(v,c))
	|    double calculateRestriction(Var v, vec<Clause*>& clist, int& count)    i.e. \sum_clist \sigma(v,c)
	|    double calculateNeed(Clause& c, Var v)                                 i.e. \sigma(v,c)
	|	 double calculateRatioBP(Var v, Agg aggregator)
	|	 double calculateRatioSP(Var v, Agg aggregator)
	|
	|  EIH 12 JAN 08
	|________________________________________________________________________________________________@*/
	bool Solver::computeSurvey()
	{
		double max_change = 0.0;
		double change = 0.0;
		int iter = 0;
		Var v = 0;
		double timefoo = 0;

		double start_time = getCpuTime();

		// Iterate until convergence or a timeout.
		//
		do
		{
			if (verbosity >= 3)
				reportf(".");
			max_change = 0.0;

			// Go through all the variables.
			//
			for (v = 0; v < nVars(); v++)
			{

				// Skip any variables that are already fixed.
				if (!(toLbool(assigns[v]) == l_Undef && decision_var[v]))
					continue;

				// Update each unfixed variable according to the appropriate rule.
				//
				switch (heuristic_mode)
				{
				case heuristic_BP:
					change = updateBP(v);
					break;
				case heuristic_EMBPL:
					change = updateEMBPL(v);
					break;
				case heuristic_EMBPG:
					change = updateEMBPG(v);
					break;
				case heuristic_EMBPGV2:
					change = updateEMBPGV2(v);
					break;
				case heuristic_SP:
					change = updateSP(v);
					break;
				case heuristic_EMSPL:
					change = updateEMSPL(v);
					break;
				case heuristic_EMSPG:
					change = updateEMSPG(v);
					break;
				case heuristic_EMSPGV2:
					change = updateEMSPGV2(v);
					break;
				case heuristic_CC:
					change = updateCC(v);
					break;
				default:
					assert(false);
				}

				// Keep track of maximal change in order to test for convergence.
				if (change > max_change)
					max_change = change;
			}
			timefoo = (getCpuTime() - start_time);
			//		printf("time out is %f, now it is %f \n", (double)timeout, timefoo);

		} while (max_change > EPSILON && iter++ < MAXITERATIONS && timefoo < (double)timeout);

		// Update stats, print log messages, and return flag according to convergence or timeout.
		//
		surveys_time += (getCpuTime() - start_time);
		surveys_attempted++;
		if (max_change <= EPSILON)
		{
			surveys_converged++;
			if (verbosity >= 3)
				reportf(":-)\n");
			if (verbosity >= 4)
				reportBiases();
			return true;
		}
		else
		{
			if (verbosity >= 3)
				reportf("[%f]:-(\n", max_change);
			if (verbosity >= 4)
				reportBiases();
			return false;
		}
	}

	// UPDATEBP: Update a variable according to BP.
	double Solver::updateBP(Var v)
	{
		int nneg, npos; // number of active negative and positive clauses for v.
		double alpha = calculateFreedom(v, neg_clauses[v], nneg);
		double beta = calculateFreedom(v, pos_clauses[v], npos);

		return determineRatioBP(v, nneg, npos, alpha, beta);
	}

	// UPDATEEMBPL: Update a variable accoring to EMBP-L.
	double Solver::updateEMBPL(Var v)
	{
		int nneg, npos; // number of active negative and positive clauses for v.
		double sum_alpha = calculateRestriction(v, neg_clauses[v], nneg);
		double sum_beta = calculateRestriction(v, pos_clauses[v], npos);

		double c = double(nneg + npos);
		return determineRatioBP(v, nneg, npos, c - sum_alpha, c - sum_beta);
	}

	// UPDATEEMBPG: Update a variable accoring to EMBP-G.
	double Solver::updateEMBPG(Var v)
	{
		int nneg, npos; // number of active negative and positive clauses for v.
		double alpha = calculateFreedom(v, neg_clauses[v], nneg);
		double beta = calculateFreedom(v, pos_clauses[v], npos);

		return determineRatioBP(v, nneg, npos, (double)nneg * alpha + (double)npos, (double)npos * beta + (double)nneg);
	}

	// UPDATEEMBPGV2: Update a variable accoring to EMBP-G-V2.
	double Solver::updateEMBPGV2(Var v)
	{
		int nneg, npos; // number of active negative and positive clauses for v.
		double alpha = calculateFreedom(v, neg_clauses[v], nneg);
		double beta = calculateFreedom(v, pos_clauses[v], npos);

		double c = double(nneg + npos);
		return determineRatioBP(v, nneg, npos, c * alpha, c * beta);
	}

	// UPDATESP: Update a variable according to SP.
	double Solver::updateSP(Var v)
	{
		int nneg, npos;	   // number of active negative and positive clauses for v.
		double rho = 0.95; // important to prevent the trivial core, as in for instance, Maneva '06

		double alpha = calculateFreedom(v, neg_clauses[v], nneg);
		double beta = calculateFreedom(v, pos_clauses[v], npos);

		double prod = alpha * beta;
		return determineRatioSP(v, nneg, npos, alpha - rho * prod, beta - rho * prod, prod);
	}

	// UPDATEEMSPL: Update a variable accoring to EMSP-L.
	double Solver::updateEMSPL(Var v)
	{
		int nneg, npos; // number of active negative and positive clauses for v.
		double sum_alpha = calculateRestriction(v, neg_clauses[v], nneg);
		double sum_beta = calculateRestriction(v, pos_clauses[v], npos);

		double c = double(nneg + npos);
		return determineRatioSP(v, nneg, npos, c - sum_alpha, c - sum_beta, c - sum_alpha - sum_beta);
	}

	// UPDATEEMSPG: Update a variable accoring to EMSP-G.
	double Solver::updateEMSPG(Var v)
	{
		int nneg, npos; // number of active negative and positive clauses for v.
		double alpha = calculateFreedom(v, neg_clauses[v], nneg);
		double beta = calculateFreedom(v, pos_clauses[v], npos);

		return determineRatioSP(v, nneg, npos, (double)nneg * alpha + (double)npos * (1.0 - beta),
								(double)npos * beta + (double)nneg * (1.0 - alpha),
								(double)(npos + nneg) * alpha * beta);
	}

	// UPDATEEMSPGV2: Update a variable accoring to EMSP-G-V2.
	double Solver::updateEMSPGV2(Var v)
	{
		int nneg, npos; // number of active negative and positive clauses for v.

		double alpha = calculateFreedom(v, neg_clauses[v], nneg);
		double beta = calculateFreedom(v, pos_clauses[v], npos);

		double prod = alpha * beta;
		double c = double(nneg + npos);
		return determineRatioSP(v, nneg, npos, c * (alpha - prod), c * (beta - prod), c * prod);
	}

	// UPDATECC: Update a variable accoring to clause count, this is a simple control heuristic.
	double Solver::updateCC(Var v)
	{
		double junk;
		int nneg, npos; // number of active negative and positive clauses for v.

		// (We don't care about the actual freedom, we just want to get clause counts.)
		junk = calculateFreedom(v, pos_clauses[v], npos);
		junk = calculateFreedom(v, neg_clauses[v], nneg);

		// (Again we use existing infrastructure, to handle boundary conditions.  Weights are just #'s of clauses.)
		return determineRatioBP(v, nneg, npos, npos, nneg);
	}

	// CALCULATEFREEDOM: Calculate probability (by _product_) that no clause in the given list needs variable v;
	// as a SIDE EFFECT, store the number of active clauses in the list, in given address.  This corresponds to the
	// product over the clauses of 1 - sigma.
	double Solver::calculateFreedom(Var v, vec<Clause *> &clist, int &count)
	{
		double retval = 1.0;
		int retcount = 0;

		int i;
		Clause *c = NULL;

		// Loop through the clauses.
		for (i = 0; i < clist.size(); i++)
		{
			c = clist[i];

			// Skip any clauses that were killed. (Sometimes clauses are killed at the top level by simplify().) (Actually this can't happen anymore.  EIH 15 FEB 08)
			// if (c == NULL) continue;

			// Skip any clauses that are empty or already satisfied. (Apparently Minisat doesn't ever create empty clauses during its operation.)
			// FUTURE may be faster to turn this off--no need to keep checking, but decrease in accuracy.
			if (satisfied(*c))
				continue;

			// So here we have an active clause; update counter.
			retcount++;

			// Update product of probabilities.  Remember that we want the prob that the clause
			// _doesn't_ need the var's help, so we subtract the running product from one.
			retval *= (1.0 - calculateNeed(*c, v));
		}

		count = retcount; // SIDE EFFECT
		return retval;
	}

	// CALCULATERESTRICTION: Calculate _sum_ of probabilities that the clauses in the given list need variable v;
	// as a SIDE EFFECT, store the number of active clauses in the list, in given address.  This corresponds to the
	// sum of sigmas over all the clauses.
	double Solver::calculateRestriction(Var v, vec<Clause *> &clist, int &count)
	{
		double retval = 0.0;
		int retcount = 0;

		int i;
		Clause *c = NULL;

		// Loop through the clauses.
		for (i = 0; i < clist.size(); i++)
		{
			c = clist[i];

			// Skip any clauses that were killed. (Sometimes clauses are killed at the top level by simplify().) (Actually this can't happen anymore.  EIH 15 FEB 08)
			// if (c == NULL) continue;

			// Skip any clauses that are empty or already satisfied. (Apparently Minisat doesn't ever create empty clauses during its operation.)
			// FUTURE may be faster to turn this off--no need to keep checking, but decrease in accuracy.
			if (satisfied(*c))
				continue;

			// So here we have an active clause; update counter.
			retcount++;

			// Update sum of probabilities.
			retval += calculateNeed(*c, v);
		}

		count = retcount; // SIDE EFFECT
		return retval;
	}

	// CALCULATENEED: Calculate probability that the clause c needs variable v for support.  This corresponds to
	// \sigma in the writeup.  It's the product of wrong-way biases for the active literals in the clause (i.e.
	// the event that all the positive variables turned out negative, and all the negative literals turned out
	// positive.)  In other words, returns probability that "v is the sole-support of c."
	double Solver::calculateNeed(Clause &c, Var v)
	{
		int i;
		double prob_need = 1.0;
		Lit p = lit_Undef;

		// Loop through the clause's literals.
		//
		for (i = 0; i < c.size(); i++)
		{
			p = c[i];

			// Skip the literal that corresponds to the given variable, as
			// well as any that represents a var that has already been fixed.
			//
			if (var(p) == v)
				continue;
			// if (value(p) != l_Undef) continue;
			if (assigns[var(p)] != toInt(l_Undef))
				continue;

			// Update running probability that clause needs var's help.
			//
			if (sign(p))
			{
				// negative literals all turned out positive
				prob_need *= pos_biases[var(p)];
			}
			else
			{
				// positive lits all turned out negative
				prob_need *= neg_biases[var(p)];
			}
		}

		return prob_need;
	}

	// DETERMINERATIOBP: Update a two-state variable's bias by given weights; return delta to test for convergence.
	// (Check boundary conditions via nneg and npos.)
	double Solver::determineRatioBP(Var v, int nneg, int npos, double wplus, double wminus)
	{
		double old_pos; // for figuring out how much the bias changed
		double new_pos; // caches a temporary value to avoid recalculation

		// Actually with the new code, I think boundary conditions can't happen.
		new_pos = wplus / (wplus + wminus);

		bool lacks_pos_clause = (npos == 0) ? true : false; // for handling boundary conditions
		bool lacks_neg_clause = (nneg == 0) ? true : false;

		// Handle boundary conditions, or else actually determine the bias by ratio in their absence.
		//
		if (wplus + wminus == 0.0)
		{
			// Prevents divide-by-zero's; this usually would also invoke lacks_pos_clause and lacks_neg_clause,
			// but might also result from pure neighbor literals that themselves lack certain polarities of clauses.
			wplus = 1.0;
			wminus = 1.0;
			// reportf("Hey there was a divide by zero problem with EM.\nExiting.\n"); exit(-1);
		}

		if (lacks_pos_clause && lacks_neg_clause)
		{
			// Technically this case is not handled by BP, but we know that half the solutions have v set
			// positively and the other half have it set negatively.
			new_pos = 0.5;
		}
		else if (lacks_pos_clause)
		{
			// By BP's assumption that every var is a sole support, we must conclude v is constrained to be neg.
			new_pos = 0.0;
		}
		else if (lacks_neg_clause)
		{
			// By BP's assumption that every var is a sole support, we must conclude v is constrained to be pos.
			new_pos = 1.0;
		}
		else
		{
			// In the absence of a boundary condition, the variable must balance its responsibilities to both
			// positive and negative clauses.  So, to do a regular update, the variable's new positive bias is:
			// wplus / (wplus + wminus).  That is, the probability that it must be positive is estimated as the
			// probability that it doesn't have to be negative, normalized.  And bcs we are just doing the
			// two-state, BP version of things, the chance of being negative is just 1 minus that.
			new_pos = wplus / (wplus + wminus);
		}

		// Store the resulting values.
		old_pos = pos_biases[v];
		pos_biases[v] = new_pos;
		neg_biases[v] = 1.0 - new_pos;

		// Finally, return a float measuring how much this variable actually moved.
		return fabs(old_pos - new_pos);
	}

	// DETERMINERATIOSP: Update a three-state variable's bias by given weights; return delta to test for convergence.
	// (Check boundary conditions via nneg and npos.)
	double Solver::determineRatioSP(Var v, int nneg, int npos, double wplus, double wminus, double wstar)
	{
		double old_pos; // for figuring out how much the bias changed
		double old_neg;
		double new_pos; // caches a temporary value to avoid recalculation
		double new_neg;

		// reportf("DETERMINERATIOSP(%d): %f %f %f \n", v + 1, wplus, wminus, wstar);

		// Actually with the new code, I think boundary conditions can't happen.
		new_pos = wplus / (wplus + wminus + wstar);
		new_neg = wminus / (wplus + wminus + wstar);

		bool lacks_pos_clause = (npos == 0) ? true : false; // for handling boundary conditions
		bool lacks_neg_clause = (nneg == 0) ? true : false;

		// Handle boundary conditions, or else actually determine the bias by ratio in their absence.
		//
		if (wplus + wminus + wstar == 0.0)
		{
			// Prevents divide-by-zero's; this usually would also invoke lacks_pos_clause and lacks_neg_clause,
			// but might also result from pure neighbor literals that themselves lack certain polarities of clauses.
			wplus = 1.0;
			wminus = 1.0;
			wstar = 1.0;
			// reportf("Hey there was a divide by zero problem with EM.\nExiting.\n"); exit(-1);
		}

		if (lacks_pos_clause && lacks_neg_clause)
		{
			// reportf("case1");
			//  Variable is unconstrained, put all its weight on star.
			new_pos = 0.0;
			new_neg = 0.0;
		}
		else if (lacks_pos_clause)
		{
			// reportf("case2");
			//  Variable is either constrained negative, or it's star.  FIXTHIS: can divide-by-zero happen?
			new_pos = 0.0;
			new_neg = wminus / (wminus + wstar);
		}
		else if (lacks_neg_clause)
		{
			// Variable is either constrained positive, or it's star.  FIXTHIS: can divide-by-zero happen?
			// reportf("case3");
			new_pos = wplus / (wplus + wstar);
			new_neg = 0.0;
		}
		else
		{
			// reportf("case4");
			//  In the absence of a boundary condition, the variable must balance its three responsibilities
			//  by normalizing over its weights.  Because we are doing three states for SP, we can just
			//  determine the chances of being positive and negative; the chances of being star are just 1 minus
			//  these two quantities.
			new_pos = wplus / (wplus + wminus + wstar);
			new_neg = wminus / (wplus + wminus + wstar);
		}

		// reportf(":%f %f %f\n", new_pos, new_neg, 1.0 - new_pos - new_neg);
		//  Store the resulting values.
		old_pos = pos_biases[v];
		old_neg = neg_biases[v];
		pos_biases[v] = new_pos;
		neg_biases[v] = new_neg;
		star_biases[v] = 1.0 - new_pos - new_neg;

		// Finally, return a float measuring how much this variable actually moved.
		return fabs(old_pos - new_pos) + fabs(old_neg - new_neg);
	}

	///////////////////////////////
	///////////////////////////////
	/////  SCORING FUNCTIONS
	///////////////////////////////
	///////////////////////////////

	// SCOREVARSKEW: Score a survey variable by the difference between its positive and negative readings.
	double Solver::scoreVarSkew(Var v)
	{
		return fabs(pos_biases[v] - neg_biases[v]);
	}

	///////////////////////////////
	///////////////////////////////
	/////  MISCELLANEOUS
	///////////////////////////////
	///////////////////////////////

	// COMPUTELCSURVEY: Compute a survey for the 'LC' heuristic based on an extremely simplified version of Dubois and Dequen's method.  Only works for 3-SAT.
	bool Solver::computeLCSurvey()
	{
		double start_time = getCpuTime();

		int v;
		int i, j;
		Lit l;
		int p1;
		int prod;
		int countneg, countpos;

		// Loop through the variables.
		for (v = 0; v < nVars(); v++)
		{

			// Loop through the current variable's unsatisfied positive clauses.
			countpos = 0;
			for (i = 0; i < pos_clauses[v].size(); i++)
			{
				Clause &c = *pos_clauses[v][i];
				if (satisfied(c))
					continue;

				// We will only process three-clauses.
				if (c.size() == 3)
				{
					// Go through the clause, skipping the literal for the current variable.
					prod = 1;
					for (j = 0; j < c.size(); j++)
					{
						l = c[j];
						if (var(l) == v)
							continue;
						p1 = (sign(l) ? pos_clauses[var(l)].size() : neg_clauses[var(l)].size());
						prod *= p1;
					}
					countpos += prod;
				}
			}

			// Loop through the current variable's unsatisfied negative clauses.
			countneg = 0;
			for (i = 0; i < neg_clauses[v].size(); i++)
			{
				Clause &c = *neg_clauses[v][i];
				if (satisfied(c))
					continue;

				// We will only process three-clauses.
				if (c.size() == 3)
				{
					// Go through the clause, skipping the literal for the current variable.
					prod = 1;
					for (j = 0; j < c.size(); j++)
					{
						l = c[j];
						if (var(l) == v)
							continue;
						p1 = (sign(l) ? pos_clauses[var(l)].size() : neg_clauses[var(l)].size());
						prod *= p1;
					}
					countneg += prod;
				}
			}

			// Update biases for current variable.
			if (countpos + countneg == 0)
			{
				pos_biases[v] = 0.5;
				neg_biases[v] = 0.5;
				star_biases[v] = 0.0;
			}
			else
			{
				pos_biases[v] = (double)countpos / ((double)countpos + (double)countneg);
				neg_biases[v] = 1.0 - pos_biases[v];
				star_biases[v] = 0.0;
			}
		}

		// Compile statistics and output desired diagnostics.
		surveys_time += (getCpuTime() - start_time);
		surveys_attempted++;
		surveys_converged++;
		if (verbosity >= 3)
			reportf(":-)\n");
		if (verbosity >= 4)
			reportBiases();

		return true;
	}

	// DOCJM: Alternate mode where we just run one survey and give output for Christian's code.
	void Solver::doCJM()
	{
		bool survey_converged = computeSurvey();

		// Even survay is not converge, still need report something

		/*	if (!survey_converged) {
				reportf("Error: survey did not converge.  Exiting.\n");
				exit(-1);
			}
		*/
		Var v;

		//	for (v = 0; v < nVars(); v++) {
		// if (assigns[v] == toInt(l_Undef)) {
		//			reportf("%d,%1.5f,%1.5f,%1.5f\n",
		//				v + 1, pos_biases[v], neg_biases[v], star_biases[v]);
		//} else {
		//	reportf(" %8d Fix: %2d (%1.3f)   (%1.3f)   (%1.3f)\n",
		//		v + 1, assigns[v], pos_biases[v], neg_biases[v], star_biases[v]);
		//}
		//	}
	}

	// INITIALIZEBIASES: Seed biases with random values.
	void Solver::initializeBiases()
	{
		int i;
		double pos, neg, star, sum;

		// We don't have to do anything if we're not using a solution-counting heuristic.
		if (heuristic_mode == heuristic_default)
			return;

		// Go through all the variables.
		for (i = 0; i < nVars(); i++)
		{

			// Initialize pos and neg biases if we are using a BP-type heuristic.
			if ((heuristic_mode == heuristic_BP) ||
				(heuristic_mode == heuristic_EMBPL) ||
				(heuristic_mode == heuristic_EMBPG))
			{
				pos = drand(random_seed);
				neg = drand(random_seed);
				sum = pos + neg;
				pos_biases[i] = pos / sum;
				neg_biases[i] = 1.0 - pos_biases[i];
				// pos_biases[i] = 0.5;
				// neg_biases[i] = 0.5;
				star_biases[i] = 0.0; //(not actually necessary)
			}

			// Initialize pos, neg, and star biases if we are using a SP-type heuristic.
			if ((heuristic_mode == heuristic_SP) ||
				(heuristic_mode == heuristic_EMSPL) ||
				(heuristic_mode == heuristic_EMSPG))
			{
				pos = drand(random_seed);
				neg = drand(random_seed);
				star = drand(random_seed);
				sum = pos + neg + star;
				pos_biases[i] = pos / sum;
				neg_biases[i] = neg / sum;
				star_biases[i] = 1.0 - pos_biases[i] - neg_biases[i];
			}
		}
	}

	// REPORTBIASES: Report biases for all the variables.
	void Solver::reportBiases()
	{
		Var v;
		reportf("==============================[ Survey Report ]================================\n");
		reportf("     Variable       Pos       Neg       Star\n");
		for (v = 0; v < nVars(); v++)
		{
			if (assigns[v] == toInt(l_Undef))
			{
				reportf(" %8d          %1.3f     %1.3f     %1.3f\n",
						v + 1, pos_biases[v], neg_biases[v], star_biases[v]);
			}
			else
			{
				reportf(" %8d Fix: %2d (%1.3f)   (%1.3f)   (%1.3f)\n",
						v + 1, assigns[v], pos_biases[v], neg_biases[v], star_biases[v]);
			}
		}
		// reportf("===============================================================================\n");
	}

}
