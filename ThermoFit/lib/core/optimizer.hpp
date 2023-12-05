
#ifndef OPTIMIZER_HPP

#define OPTIMIZER_HPP


#pragma once
#include <iostream>
#include "options.hpp"

class Optimizer {
public:

	// Pointer to relative options
	options* opts_ptr;

	// Solves problem
	virtual void solve() { std::cerr << "solve() not overridden!" << std::endl; }

	// Returns outflag value
	virtual int outflag() { return outflag_; }

	// Returns solution
	//std::vector<double> solution() {};

	// Sets new Problem dimensions and constraints
	virtual void setDimensions(std::vector<std::vector<double>> constr, int dim) { };

	// Returns zero_function value of solution
	virtual double solution_function() { std::cerr << "solution_function() not overridden!" << std::endl; return DBL_MAX; };

	// Returns number of function evaluations
	virtual int function_calls() { std::cerr << "function_calls() not overridden!" << std::endl;  return 0; };

	// Groups of point
	int groups = 1;

	friend class RL;

	// Inherited classes must be befriended to get accsess private members and functions 
	friend class PSO;
	friend class GA;
private:

	// 0 .... solution found
	// 1 .... max_Iter reached
	// 2 .... max_FunCalls reached
	// 3 .... function returned nan
	// 4 .... function returned inf
	// x .... others
	// 10 ... function stability reached
	// 61 ... still training
	// 62 ... training finished succsessful
	// 65 ... still running
	// 66 ... running finished succsessful
	// Solves problem with RL enhancement
	virtual std::vector<double> solve_rl(int iter) { std::cerr << "solve_rl() not overridden!" << std::endl; return {}; }

	// Indicator if solution was found
	bool done = 0;

	// Outflag for optimizer status
	int outflag_ = 99;
};
#endif // !OPTIMIZER_HPP
