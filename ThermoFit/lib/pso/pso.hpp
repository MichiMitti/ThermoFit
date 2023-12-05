/*
//Particle Swarm Optimizer
//     ____  _____ ____ 
//    / __ \/ ___// __ \
//   / /_/ /\__ \/ / / /
//  / ____/___/ / /_/ / 
// /_/    /____/\____/                 
//
@author: Michael Mitterlindner, TU Graz
@email: mitterlindner@gmx.at
@github: MichiMitti
@year: 2023
*/

#ifndef PSO_HPP

#define PSO_HPP

#define _CRT_SECURE_NO_WARNINGS		// for ctime in createCSV

#pragma once
#include <iostream>		// write to console
#include <fstream>		// write to file
#include <chrono>		// high precision time
#include <ctime>		// current time and date
#include <filesystem>	// manage files
#include <vector>		// vectors
#include <random>		// random number engine and distribution
#include <numeric>		// for accumulate function
#include "../core/global.hpp"

#include "../core/zero_function.hpp"
#include "../core/options.hpp"
#include "../core/optimizer.hpp"
#include "networks_psorl.hpp"

//===== To Do =====\\ 
//* test getVals function
//* convergence flag
//* check if random first velocity is better or random with constraints or as it is (with velo calc with initialized values)
//* get output dimensions of zero_function to size fvals right
//* test outflags for nan and inf (add breaks)

struct options_pso : options {
	std::string options_typ = "pso";								// stores options_powell name for checking if right options_powell are used in CallClass

	double tolerance = 1e-08;										// tolerance for solution
	int max_Iter = 1001;											// maximal Iterations

	int maxFunEvals = 20000;										// maximal function evaluations
	//int cnt_calfun = 0;												// counter for function evaluations

	int swarm_size = 30;											// number of Particles
	std::vector<std::vector<double>> constraints = { {0.0, 20.0} };	// constraints variable space (upper and lower bound) applied to all dimensions

	double omega_start = 0.7;										// starting inertia weight linearly decreast from _start to _end typically [1.2, 0.9]
	double omega_end = 0.7;											// end inertia weight after max_Iter
	double c1 = 1.5;												// cognitive weitht (particle)
	double c2 = 1.5;												// social weight (swarm)

	double velocity_max = 0.5;										// clips maximal particle velocity (multiplied with range)

	bool debug = 0;													// select debug mode

	bool write_CSV = 0;												// writes steps to a csv file
	std::string filename = "log.csv";								// filename for saveing
	std::string problem_description = " ";							// string for short problem description

	std::vector<std::vector<double>> starting_points = { {} };		// Starting point matrix
};

class Particle {
public:
	friend class PSO;						// PSO can acess all private members

private:
	std::vector<double> position;			// current position of the particle
	std::vector<double> velocity;			// current velocity of the particle
	double fvalue;							// stores current value of particle

	std::vector<double> best_position;		// best particle position in particle
	double best_fvalue;						// best function value of particle
};


/* Particle Swarm Optimizer
with dampening wall constraints
 swarm_size ... number of particles
*/
class PSO : public Optimizer {

public:
	// Ctor
	PSO();
	// Ctor
	PSO(Zero_Function* fun, std::vector<double>* initGuess, options_pso* opts);

	// Reads in PSO values so it has not to be created each Call
	void getVals(Zero_Function* fun, std::vector<double>* initGuess, options_pso* opts);

	// Solves problem
	void solve();

	// Returns solution
	std::vector<double> solution();

	// Returns zero_function value of solution
	double solution_function();

	// Returns number of function evaluations
	int function_calls();

	// Calls zero function
	void calfun(std::vector<double>* vec, std::vector<double>* fvals, options_pso* opts);

	friend class RL;

private:

	Zero_Function* fun_;								// Zero function

	void solve_std();
	void solve_rl();

	std::vector<double> initGuess_;						// Initial guess

	options_pso opts_;									// Stores option struct with hyper-parameters

	int N_;												// Number of dimensions

	std::vector<std::vector<double>> constraints_;		// Vector for constraints

	int swarm_size_;									// Number of particles in swarm

	void initializeSwarm();								// Initializes swarm
	void setDimensions(std::vector<std::vector<double>> constr, int dim);
	void setConstraints();								// Sets constraints to Optimizer

	void updateParticle(int p);							// Updates partile value, position, velocity,...
	void calcPosition(int p);							// Calculates particle velocity

	double omega_;										// Inerta weight
	double c1_;											// Cognitive weitht (particle)
	double c2_;											// Social weight (swarm)
	double velocity_max_;								// Clips maximal particle velocity (multiplied with range)
	std::vector<double> velocity_max_vec;				// Holds maximal velocities of each dimension (constraint_[1]-constraint_[0])*velocity_max_

	std::vector<Particle> Swarm;						// Particle vector *swarm*

	double best_swarm_value;							// Stores best particle value for swarm
	std::vector<double> best_swarm_position;			// Stores best position for swarm
	double old_best_swarm_value;						// Stores old best particle value for swarm

	std::vector<double> fvals;							// Stores zero function values
	double fval;										// Stores zero function sum

	std::default_random_engine engine;					// Random engine std::mt19937
	std::uniform_real_distribution<double> rand_dbl;	// Double random number between 0 and 1

	int iter;											// Iteration counter
	int i_break;										// Function breaks counter

	int getVal_funcCtr_;								// Gets Values (so no ctor must be called) NOT YET IMPLEMENTED

	void createCSV();									// Creates csv file
	void writeParticleToCSV();							// Writes particle positions to a csv file
	void closeCSV();									// Writes end message to csv file

	std::string agent_file;								// path to csv file
	std::fstream fout;									// Contains file stream (fileout)
	char delimiter;										// Delimiter for writing to file

	std::chrono::high_resolution_clock::time_point start_time;	// Start time for stopwatch
	std::chrono::high_resolution_clock::time_point end_time;	// End time for stopwatch

	void solve_rl_step(int step);		// Solves for RL method (gives back state)
	void calcState();									// Calculates state
	void transformState();
	std::vector<double> state;							// state vector
	std::vector<double> transformed_state;				// state vector

	int groups;											// Particle groups for rl
	std::vector<std::vector<double>> hyper_params;		// Hyperparameter vector from rl

	double diversity;									// Diversity (average distance between all particles)
	std::vector<double> swarm_mean;						// Mean swarm position
	double distance_to_mean;							// Distance between each particle and the mean swarm postion
	double iteration;									// Normalized iteration = iter / max_Iter
	int iter_last_improvement;							// Iteration where the best_swarm_value was improved
	double no_improvement;								// Normalize number of iteration when no improvement = (iter - iter_last_improvement)/max_Iter
	Network* Net;										// Agent network
	Network_1 Net_1;									// Agent network1
	Network_2 Net_2 ;									// Agent network2
	Network_3 Net_3;									// Agent network3

	void selectNet();

	void getAction();
	void convertAction();

	std::vector<double> action_vector;
	std::vector<std::vector<double>> action;
};




#endif // !PSO_HPP
