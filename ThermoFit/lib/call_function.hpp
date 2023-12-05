/*
Legendary CallClass

@author: Michael Mitterlindner, TU Graz
@email: mitterlindner@gmx.at
@github: MichiMitti
@year: 2023
*/

#ifndef CALL_FUNCTION_HPP

#define CALL_FUNCTION_HPP
#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <chrono>
#include "core/global.hpp"

// Available Methods
#include "pso/pso.hpp"
#include "powell/powell.hpp"
#include "simplex/simplex.hpp"
#include "ga/ga.hpp"

// Includes Zero_Function class and transfer_struct
#include "core/zero_function.hpp"
#include "core/optimizer.hpp"
#include "core/options.hpp"


//===== To Do =====\\ 
//* write results to log file
//* use getVAls function so object must not be created each time
//* Use options ptr instead of void ptr
//* 


// Available optimizers:
// pso, pso_rl, powell and simplex 
// I posses two Ctor if you have a single problem please pass one pointer to your problem in my class call
// if you have two problems please pass two pointers first the Inner pointer second the Outer.
// You can pass a long boundary vector right in the ctor and it will then be seperated for the Inner and Outer boundaries.
// After this specify how you want your problems to be solved therefore use for a single problem the function:
// .solve_Inner_problem_with and for two problems please first also call 
// .solve_Inner_problem_with and then .solve_Outer_problem_with
// it is important that you call your Inner problem in the Outer problem with this line of code:
// call_class_ptr->call_inner_problem(call_class_ptr);
class CallClass {
public:

	//CallClass(); //Ctor for useing inherited functions (currently not necessary)
	
			// Ctor for single problem
	CallClass(Zero_Function* Fun_Inner, std::vector<std::vector<double>>* boundaries = nullptr);
		
			// Ctor for two problems
	CallClass(Zero_Function* Fun_Inner, Zero_Function* Fun_Outer, std::vector<std::vector<double>>* boundaries = nullptr);

			// Dtor 
	~CallClass();

			// Declare method as string to solve the INNER problem options_powell are: "pso", "powell", "simplex", "pso_rl"
			// If no boundaries are given the boundaries will be used from the class call or default values will be used
			// If initGuess_Inner is { } the USER is asked to specify the number of dimensions and the initial guess will be searched by randomly
			// but other methods can be choosen with the set_initalGuess_options() function
			// If no options_powell are specified the standart options_powell are used
	void solve_Inner_problem_with(std::string method_Inner, std::vector<std::vector<double>>* boundaries = nullptr, std::vector<double>* initGuess_Inner = nullptr, void* opts_Inner = nullptr);

			// Declare method as string to solve the OUTER problem options_powell are: "powell", "simplex", "pso", "pso_rl"
			// If no boundaries are given the boundaries will be used from the class call or default values will be used
			// If initGuess_Inner is { } the USER is asked to specify the number of dimensions and the initial guess will be searched randomly
			// but other methods can be choosen with the set_initalGuess_options() function (if for Inner and Outer problem no inital guess provided
			// it will be searched randomly for it)
			// If no options_powell are specified the standart options_powell are used
	void solve_Outer_problem_with(std::string method_Outer, std::vector<std::vector<double>>* boundaries = nullptr, std::vector<double>* initGuess_Outer = nullptr, void* opts_Outer = nullptr);

			// Solves the declared problem(s)
	void solve();

	std::vector<double> result();
	double function_value();
	int function_evaluations();
	
			// Function should be called from Outer Function to solve Inner Function returns a vector
			// The syntax is: "call_class_ptr->call_inner_problem();" the result is set in the transfer_struct
			// for this the following two header files must be included: 
			// "lib/call_function.hpp" and "lib/core/zero_function.hpp"
	void call_inner_problem();

			// This method searches an inital guess randomly
			// For finding Inner guess please set all dependencies most important the inital guess for the outer problem 
			// set_Outer_initial_guess_for_single_problem(&initGuess_Outer); (the guess will be passed to the transfer_struct_.outer)
			// For finding Outer guess the dependencies must be set with the following methods:
			// set_initalGuess_options("pso");	
			// set_Boundary_conditions(&boundaries_inner, "inner");
			// set_Inner_initial_guess_for_single_problem(&initGuess_Inner);
	void random_guess_search(std::vector<std::vector<double>>* boundaries, std::string problem = "inner");

			// This method searches an inital guess randomly for two specified problems within boudaries randomly
			// (the Inner zero_function will not be solved!)
	void random_guess_search(std::vector<std::vector<double>>* boundaries_inner, std::vector<std::vector<double>>* boundaries_outer);

			// Sets base for how many guesses for random searched will be used guesses = base^dimensions
	void set_base_for_random_guess_search(int base);

			// Sets Outer initial guesses for a single problem if Inner problem is dependent on them (single problem)
	void set_Outer_initial_guess(std::vector<double>* initGuess_Outer);
			// Sets Inner initial guesses for a single problem if Outer problem is dependent on them
	void set_Inner_initial_guess(std::vector<double>* initGuess_Inner);

			// Sets Boundary conditions if no boundaries are given std boundaries are used but applied only for one dimension
			// Specifiy with "inner" or "outer" wich boundaries you want to set. Default is "inner".
	void set_Boundary_conditions(std::vector<std::vector<double>>* boundaries = nullptr, std::string problem = "inner");


			// Sets the Options for Inital Guess calculation (e.g. for solveing Inner problem)
			// method..."random" sets random method, "pso" sets pso method, "inner" sets method for solveing nested Inner problem!
			// options_Int_Guess ... solver name to which options_powell should be set
			// options_initGuess ... options_powell
	void set_inital_guess_options(std::string method = "random", std::string options_Init_Guess = "pso", void* options_initGuess = nullptr);

			// Sets special_variables in transfer_struct
	void set_transfer_struct_special_variables(std::vector<double> special_variables);

			//  forward transformation with the sigmoid function whithin the Inner boundaries
	std::vector<double> sigmoid_forward_Inner(std::vector<double> vec);
			//  forward transformation with the sigmoid function whithin the Outer boundaries
	std::vector<double> sigmoid_forward_Outer(std::vector<double> vec);
			//  backward transformation with the sigmoid function whithin the Inner boundaries
	std::vector<double> sigmoid_backward_Inner(std::vector<double> vec);
			//  backward transformation with the sigmoid function whithin the Outer boundaries
	std::vector<double> sigmoid_backward_Outer(std::vector<double> vec);

			// Squares and summs up a vector
	double solution_conversion_n_to_1(std::vector<double>* vec);

			// Displays results for Inner and Outer problem
	void display_results();	

			// Displays results for inital guess calculation
	void display_inital_guess_results();

			// Returns outflag
	int outflag(std::string problem = "inner");

	void display_outflag(std::string problem = "inner");

			// Displays banner
	void displayBanner(std::string banner = "ThermoFit");

	// Saves matrix to csv file
	void printVectorToFile(std::vector<std::vector<double>>* matrix, std::string file_name = "file", std::vector<std::vector<std::string>>* header = nullptr, std::string extension = ".csv", char delimiter = ';');
	void printVectorToFile(std::vector<double>* matrix, std::string file_name = "file", std::vector<std::vector<std::string>>* header = nullptr, std::string extension = ".csv", char delimiter = ';');

	double shift_for_std_boundary;					// Shifts std boundary conditions std={0,1}
	double multiplier_for_std_boundary;				// Scales std boundary conditions std={0,1}

	std::vector<double> vec_Inner;					// vector of Inner solution
	std::vector<double> vec_Outer;					// vector of Outer solution

	double fval_Inner;								// Function values of Inner solution (of zero_function)
	double fval_Outer;								// Function values of Outer solution (of zero_function)

	int fun_eval_Inner;								// Number of function evaluations
	int fun_eval_Outer;								// Number of function evaluations

	std::string error_massage;						// Error string to find were problem is

	double solveing_time;							// Time needed to solve specifed problem in milliseconds

private:

	Zero_Function* Fun_Inner_;								// Pointer to the Inner zero function
	Zero_Function* Fun_Outer_;								// Pointer to the Outer zero function

	void solve_Inner();										// Solves Inner problem
	void solve_Outer();										// Solves Outer problem

	transfer_struct transfer_struct_;						// Holds transfer variables used in the Inner and Outer problem

	int multiple_problems;									// | 0...single problem | 1...multiple problems |
	int problems_Specified;									// Determines if solve_..._Problem_with() was called (0...not specified, 1...only Inner, 2...Inner and Outer)

	std::vector<std::vector<double>>* boundary_;			// Passed Boundary conditions first entries are for Inner problem second are for Outer
	std::vector<std::vector<double>> boundary_Inner;		// Inner boundary conditions
	std::vector<std::vector<double>> boundary_Outer;		// Outer boundary conditions

	std::vector<double> initGuess_Inner_calc;				// Holds Inner inital guess if it is necessary to calculate
	std::vector<double> initGuess_Outer_calc;				// Holds Outer inital guess if it is necessary to calculate

	std::vector<double>* initGuess_Inner_;					// (vector)Pointer to the Inner inital guess
	std::vector<double>* initGuess_Outer_;					// (vector)Pointer to the Outer inital guess

	int dimensions_Inner;									// Hold Inner dimensions
	int dimensions_Outer;									// Hold Outer dimensions

	int method_Inner_save_;									// | 0...pso | 1...powell | 2...simplex |  10...pso_rl(from solve_..._problem_with())
	int method_Outer_save_;									// | 0...pso | 1...powell | 2...simplex |  10...pso_rl(from solve_..._problem_with())

	int method_Inner_;										// | 0...pso | 1...powell | 2...simplex | 3...ga | 10...pso_rl
	int method_Outer_;										// | 0...pso | 1...powell | 2...simplex | 3...ga | 10...pso_rl

	PSO* Problem_Inner_PSO;									// Heap pointer for Inner PSO (Particle Swarm Optimizer)
	PSO* Problem_Outer_PSO;									// Heap pointer for Outer PSO (Particle Swarm Optimizer)
	powell* Problem_Inner_Powell;							// Heap pointer for Inner Powell
	powell* Problem_Outer_Powell;							// Heap pointer for Outer Powell
	DHS* Problem_Inner_DHS;									// Heap pointer for Inner DHS (Downhill Simplex)
	DHS* Problem_Outer_DHS;									// Heap pointer for Outer DHS (Downhill Simplex)
	GA* Problem_Inner_GA;									// Heap pointer for Inner GA (Genetic Algorithm)
	GA* Problem_Outer_GA;									// Heap pointer for Outer GA (Genetic Algorithm)

	int outflag_Inner;										// Ouflag for Inner
	int outflag_Outer;										// Ouflag for Outer

	void calculate_initial_guess();							// Calls inital guess methods
												
	int base_random_guesses;								// Base of how many guesses will be used guesses = base^dimensions
	void calc_init_guess_pso();								// Calculates Inner inital Guess with PSO
	void calc_init_guess_random();							// Calculates Inner inital Guess randomly
	std::vector<double> best_guess_vec;						// Holds best guess from random guess search
	std::vector<std::vector<double>> best_guess_vec_2;		// Holds best guess from random guess search {{inner},{outer}}
	
	void set_transfer_struct(std::string problem);			// Method to set transfer_struct_ for inital guess calculation

	std::string method_Init_Guess_;							// Method for Init_Gues_Calculation
	void* options_initalGuess_ptr_;							// Pointer to options for calculating the inital guess

	void check_dimensions(std::string problem = "inner");			// Checks if dimensions match for inner and outer problem (dimensions_..., boundary_..., initGuess_..._, transfer_struct_.(boundary_, and passing vectors))
	void check_options(std::string method, std::string options_);	// Checks if selected method and options_powell are equal

	void* opts_Inner_;										// (void)Pointer to the Inner options struct (gets dereferenced when decided if passed options_powell or std options_powell are used)
	void* opts_Outer_;										// (void)Pointer to the Outer options struct (gets dereferenced when decided if passed options_powell or std options_powell are used)

	options_pso options_pso_std_Inner;						// Standard options for PSO optimizer (Inner because boundaries could be diffrent)
	options_pso options_pso_std_Outer;						// Standard options_powell for PSO optimizer (Outer because boundaries could be diffrent)
	options_pso* options_pso_Inner_ptr;						// Pointer to Inner PSO options
	options_pso* options_pso_Outer_ptr;						// Pointer to Outer PSO options
		
	options_powell options_powell_std;						// Standard options_powell for Powell solver
	options_powell* options_powell_Inner_ptr;				// Pointer to Inner Powell options options_powell
	options_powell* options_powell_Outer_ptr;				// Pointer to Outer Powell options options_powell

	options_dhs options_dhs_std;							// Standard options for Simplex optimizer
	options_dhs* options_dhs_Inner_ptr;						// Pointer to Inner Simplex options
	options_dhs* options_dhs_Outer_ptr;						// Pointer to Outer Simplex options

	options_ga options_ga_std_Inner;						// Standard options for GA optimizer (Inner because boundaries could be diffrent)
	options_ga options_ga_std_Outer;						// Standard options for GA optimizer (Outer because boundaries could be diffrent)
	options_ga* options_ga_Inner_ptr;						// Pointer to Inner GA options
	options_ga* options_ga_Outer_ptr;						// Pointer to Outer GA options

	int initGuess_case;										// Option for displaying results (ugly but I did not know a better way^^) (2...display inner method when calculating outer guess)

	std::string method_to_string(int method);				// Converts the method (int) into a string

	std::chrono::steady_clock::time_point time_start;		// Start time for stopwatch
	std::chrono::steady_clock::time_point time_end;			// End time for stopwatch
	std::chrono::microseconds duration;						// Time between time_start and time_end in microseconds
};

#endif // !CALL_FUNCTION_HPP