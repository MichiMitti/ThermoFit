/*
Test classes for CallClass

@author: Michael Mitterlindner, TU Graz
@email: mitterlindner@gmx.at
@github: MichiMitti
@year: 2023
*/

#pragma once
#include "parameters.hpp"
#include <vector>
#include <string>

#include "../../lib/core/zero_function.hpp"
#include "../../lib/call_function.hpp"

// Problem class for solveing a two-dimensional optimization problem
// this class must be inherit from the class Zero_Function
class Call_Test_Inner : public Zero_Function {
public:
	std::string problem_description = "=== Inner probelem (single) ===";

	void zero_function(std::vector<double>* vec, std::vector<double>* fval);

	parameters* params_ = nullptr;
};

// Problem class for solveing a two-dimensional optimization problem 
// with an two-dimensional Inner problem
// which will be called in this function
// this class must be inherit from the class Zero_Function
class Call_Test_Outer : public Zero_Function {
public:

	std::string problem_description = "=== Outer problem (two problems) ===";

	void zero_function(std::vector<double>* vec, std::vector<double>* fval);

	parameters* params_ = nullptr;
};

// Benchmark Functions:
// "sphere"
// "rastrigin"
// "rosenbrock"
class Benchmark_Functions : public Zero_Function {
public:
		// std ctor no function is set please set it with set_function(...)!
	Benchmark_Functions();	

		// Ctor also set test function with string available: "sphere", "rastrigin", "rosenbrock"
	Benchmark_Functions(std::string function_str);		

		// "sphere", "rastrigin", "rosenbrock"
	void set_function(std::string function_str);

		// Set test function with int (so all can be easly used through iteration) 1...Sphere, 2...Rastrigin, 3...Rosenbrock
	void set_function(int function_int);										

	std::string problem_description;

	parameters* params_;														// Parameters are not used

private:
	
	int function_;																// 1...Sphere, 2...Rastrigin, 3...Rosenbrock
	int dimensions_;															// Number of dimensions

	void zero_function(std::vector<double>* vec, std::vector<double>* fval);	// Inehrited from Zero_Function

		// Sphere test function		(zero at (0,0) without offset,	boundaries ]-inf, inf[)
	void sphere(std::vector<double>* vec, std::vector<double>* fval);

		// Rastrigin test function	(zero at (0,0) without offset,	boundaries [-5.12, 5.12])
	void rastrigin(std::vector<double>* vec, std::vector<double>* fval);
		
		// Rosenbrock test function	(zero at (b,b^2) (for 2d),		boundaries ]-inf, inf[)
	void rosenbrock(std::vector<double>* vec, std::vector<double>* fval);		

	const double pi = 3.14159265358979323846;									// Cirular number ;)
};
