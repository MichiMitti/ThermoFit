/*
Genetic Algorithm
//      ________
//    / ____/   |
//   / / __/ /| |
//  / /_/ / ___ |
//  \____/_/  |_|
//
@author: Michael Mitterlindner, TU Graz
@email: mitterlindner@gmx.at
@github: MichiMitti
@year: 2023
*/

#ifndef GA_HPP

#define GA_HPP

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
#include <algorithm>	// sorting algorithm
#include "../core/global.hpp"

#include "../core/zero_function.hpp"
#include "../core/options.hpp"
#include "../core/optimizer.hpp"

//===== To Do =====\\ 
//* include mutation in crossover function (could be faster)
//* crossover -> not checking if gene is already in in offspring from parent 2
//* work around return vectors for getOffspting() i think this is slow
//* convergence flag

struct options_ga : options {
	std::string options_typ = "ga";										// stores options_powell name for checking if right options_powell are used in Call_Class

	double tolerance = 1e-02;											// requested tolerance for solution
	int max_Iter = 1001;												// maximal Iterations (generations)

	int maxFunEvals = 20000;											// maximal function evaluations
	int cnt_calfun = 0;													// counter for function evaluations

	int population_size = 30;											// number of Individuals			

	std::vector<std::vector<double>> constraints = { {0.0, 10.0} };		// constraints variable space (upper and lower bound) applied to all dimensions

	double mutation_rate = 0.05;										// how likely one genome is mutated

	bool debug = 0;														// select debug mode

	bool write_CSV = 0;													// writes steps to a csv file
	std::string filename = "log.csv";									// filename for saveing
	std::string problem_description = " ";								// string for short problem description
};

class Individual {
public:
	std::vector<double> chromosome;			// holds chromosomes (Individual position)
	double fitness;							// holds fitness (zero_funcion valued)
};

class GA : public Optimizer {
public:
	// Ctor
	GA(Zero_Function* fun, std::vector<double>* initGuess, options_ga* opts);

	// Solves problem
	void solve();

	// Returns solution
	std::vector<double> solution();

	// Retruns zero_function value of solution
	double solution_function();

	// Returns number of function evaluations
	int function_calls();

	// Calls zero_function
	void calfun(std::vector<double>* vec, std::vector<double>* fval, options_ga* opts);

	friend class RL;

private:

	Zero_Function* fun_;													// Zero function

	std::vector<double> initGuess_;											// Initial guess

	options_ga opts_;														// Stores option struct with hyper-parameters

	int N_;																	// Number of dimensions
	std::vector<std::vector<double>> constraints_;							// Vector for constraints

	double tolerance_;														// Requested tolerance for solution

	int population_size_;													// Number of Individuals
	int chromosome_lenght_;													// Lengh of one chromosome (dimensions)
	int max_generations_;													// Maximal Generation (iterations)
	double mutation_rate_;													// How likely one genome is mutated

	std::vector<Individual> Population;										// Vector of Individuals
	std::vector<Individual> Next_Gen;										// Vector of next generation Individuals
	std::vector<Individual> Combined_Pop;									// Vector of Individuals (includes Population and Next_Gen)

	double best_pop_fitness;													// Stores best individual value for population
	std::vector<double> best_pop_chromosome;									// Stores best position for population
	double old_best_pop_fitness;												// Stores old best individual value for population

	Individual* Parent1;													// Holds ptr to parent
	Individual* Parent2;													// Holds ptr to parent

	Individual Offspring1;													// Holds new offspring
	Individual Offspring2;													// Holds new offspring

	std::default_random_engine engine;										// Random engine std::mt19937
	std::vector<std::uniform_real_distribution<double>> random_engine;		// Vector with engines for problem (because boundaries must be set for each dimensions)

	std::uniform_real_distribution<double> rand_dbl;						// Double random number between 0 and 1
	std::uniform_int_distribution<int> rand_int_dim;						// Integer random number from 0 to number of dimensions
	std::uniform_int_distribution<int> rand_int_pop;						// Integer random number from 0 to number of valid parents

	std::vector<double> fvals;												// Stores zero function values
	double fval;															// Stores zero function sum

	void initializePopulation();											// Initializes Population

	void evaluatePopulation();												// Calculates Populations individual fittness
	void evaluateIndividual(Individual* Ind);								// Calculates individual fittness

	void sortPopulation(std::vector<Individual>* population);				// Sorts Population (best Individual first)

	void createNextGeneration();											// Creates next Generation

	Individual* getParent();												// Selects Individual for mateing
	Individual* tournamentSelection();										// Selects parent tournament like
	Individual* biasedRandomSelection();									// Uses a biasedRandomSelection with a cummulative probability vector
	void createCumPropVec();												// Creates Cummulative probability vector [0,1] for biasedRandomSelection
	std::vector<double> cum_prop_of_pop;									// Cummulative probability vector [0,1] for biasedRandomSelection

	void getOffspring(Individual* A, Individual* B);						// Creates offspring
	Individual doCrossover(Individual* Ind_1, Individual* Ind_2);			// Selects crossover method
	Individual randomCrossover(Individual* Ind_1, Individual* Ind_2);		// Performs random crossover for each dimension r*Parent1 + (1-r)*Parent2
	Individual uniformCrossover(Individual* Ind_1, Individual* Ind_2);		// Performs crossover for each dimension individual 
	Individual singleCrossover(Individual* Ind_1, Individual* Ind_2);		// Performs crossover with single crossover point

	void mutate(Individual* individual);									// Performs mutation on individual

	int iter;																// Iteration counter
	int i_break;

	void createCSV();									// Creates csv file
	void writeIndividualToCSV();						// Writes individual positions to a csv file
	void closeCSV();									// Writes end message to csv file

	std::string agent_file;								// path to csv file
	std::fstream fout;									// Contains file stream (fileout)
	char delimiter = ',';								// Delimiter for writing to file

	std::chrono::high_resolution_clock::time_point start_time;	// Start time for stopwatch
	std::chrono::high_resolution_clock::time_point end_time;	// End time for stopwatch
};


#endif // !GA_HPP
