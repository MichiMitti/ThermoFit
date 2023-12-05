/*
Training class for RL

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
//#include "../../lib/call_function.hpp"

// Benchmark Functions:
// "sphere"
// "rastrigin"
// "rosenbrock"
// "beale"
// "camel"
class Training_Function : public Zero_Function {
public:
	// std ctor no function is set please set it with set_function(...)!
	Training_Function();

	// "sphere", "rastrigin", "rosenbrock"
	void set_function(std::string function_str);

	void resetTrainingFunction(int number);

	// Set test function with int (so all can be easly used through iteration)
	void set_function(int function_int);
	void set_offset(double offset_dbl);

	int available_functions;													// Number of Zero_Functions implemented

	std::vector<std::vector<std::vector<double>>> constraints_;								// Boundaries for Functions

	std::string problem_description;

	parameters* params_;														// Parameters are not used

	friend class RL;

private:

	int function_;																// 0...Sphere, 1...Rastrigin, 2...Rosenbrock, 3...Beale, 4...Three-Hump Camel
	int dimensions_;															// Number of dimensions
	double offset;																// Offset for training

	void zero_function(std::vector<double>* vec, std::vector<double>* fval);	// Inehrited from Zero_Function

	// 0  Sphere test function		(zero at (0,0) without offset, boundaries ]-inf, inf[)
	void sphere(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> sphere_domain;

	// 1 Rastrigin test function	(zero at (0,0) without offset, boundaries [-5.12, 5.12])
	void rastrigin(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> rastrigin_domain;

	// 2 Rosenbrock test function	(zero at (b,b^2) (for 2d), boundaries ]-inf, inf[)
	void rosenbrock(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> rosenbrock_domain;

	// 3 Beale test function	(zero at (3, 0.5) (only 2d), boundaries [-4.5, 4.5])
	void beale(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> beale_domain;

	// 4 Three-hump camel test function	(zero at (0, 0) (only 2d), boundaries [-5, 5])
	void three_camel(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> three_camel_domain;

	// 5 Ackley test function (zero at (0, 0) (only 2d), boundaries [-5, 5]
	void ackley(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> ackley_domain;

	// 6 Schwefel test function (zero at (0, 0), boundaries [-500, 500]
	void schwefel(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> schwefel_domain;

	// 7 Levy test function (zero at (1, 1), boundaries [-10, 10]
	void levy(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> levy_domain;

	// 8 Griewank test function (zero at (0, 0), boundaries [-600, 600]
	void griewank(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> griewank_domain;

	// 9 Six-hump camel test function (zero at (0.0898, -0.7126) and (-0.0898, 0.7126), boundaries {[-3, 3], [-2,2]}
	void six_camel(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> six_camel_domain;

	// 10 Alpine test function (zero at (0.0, 0.0), boundaries [-10, 10]
	void alpine(std::vector<double>* vec, std::vector<double>* fval);
	std::vector<std::vector<double>> alpine_domain;

	std::vector<std::vector<int>> dimensions;

	void setDescription();

	void displayDescription();

	const double pi;								// Cirular number ;)
	const double euler;
};