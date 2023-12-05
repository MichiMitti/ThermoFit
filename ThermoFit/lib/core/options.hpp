
#ifndef OPTIONS_HPP

#define OPTIONS_HPP

#pragma once
#include <vector>
#include "global.hpp"

// Base Class for all optimizer Options
struct options {

	double tolerance = 1e-05;	// tolerance for function falue and stability
	int rl_on = 0;				// Turns on Reinforcement Learning (if rl_on > 0 this network will be loaded)
	int max_Iter = 1001;		// Maximal iterations
	int maxFunEvals = 10000;	// Maximal function calls
	int break_tol = 100;		// Number of allowed successive violations of the termination criteria |fun_val_new-fun_val_old|<tolerance
	int cnt_calfun = 0;			// counter for function evaluations
	const unsigned int seed = SEED;

	std::vector<std::vector<double>> constraints = { {0.0, 1.0} }; // Constraints for boundary conditions

	std::vector<std::vector<double>> action;						// Holds action choosen by the RL agent
};

#endif // !OPTIONS_HPP


