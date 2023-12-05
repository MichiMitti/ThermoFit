#include "test_functions.hpp"

void Call_Test_Inner::zero_function(std::vector<double>* vec, std::vector<double>* fval) {
	// Inner zero function

	double x;
	double y;

	x = (*vec)[0];
	y = (*vec)[1];

	//Test function to find simple minima

	double a = transfer_->outer[0]; //0;
	double b = transfer_->outer[1]; //0;

	//Sphere function with offset
	(*fval)[0] = pow(x - a - 2, 2);
	(*fval)[1] = pow(y - b - 2, 2);

	//Himmelblau function (difficult because function has four roots)
	//(*fval)[0] = pow(pow(x,2) + y - (*params_).a_soll, 2.0);
	//(*fval)[1] = pow((x + pow(y, 2.0)- (*params_).b_soll), 2.0);

	// Rosenbrock function (not for multiple problem because depends only on a root(a,a^2))
	//(*fval)[0] = pow(((*params_).a_soll - x), 2.0);
	//(*fval)[1] = (*params_).b_soll * pow((y - pow(x, 2.0)), 2.0);

	if ((isnan(x)) || (isnan(y))) {
		(*fval)[0] = 1e10;
		(*fval)[1] = 1e10;
	}

	transfer_->inner[0] = x;
	transfer_->inner[1] = y;

	transfer_->boundary_inner[0]; // test if boundaries are set

	//std::cout << "=== inner input values ===" << std::endl;
	//std::cout << "(*vec)[0] " << (*vec)[0] << " | (*vec)[1] " << (*vec)[1] << std::endl;
	//std::cout << "=== inner output values ===" << std::endl;
	//std::cout << "(*fval)[0] " << (*fval)[0] << " | (*fval)[1] " << (*fval)[1] << std::endl;
	//std::cin.get();
	return;
}

void Call_Test_Outer::zero_function(std::vector<double>* vec, std::vector<double>* fval) {
	// Outer zero function

	std::vector<double> result(2, 0);

	transfer_->outer[0] = (*vec)[0];
	transfer_->outer[1] = (*vec)[1];

	transfer_->boundary_inner[0]; // test if boundaries are set
	transfer_->boundary_outer[0]; // test if boundaries are set

	//============= Inner Function Call Section =============\\ 

	//===== First variant =====\\ 
	//CallClass::call_inner_problem(ptr_);				// for this the CallClass must be inherited and a default ctor is necessary

	//===== Second variant =====\\  
	//(*ptr_).call_inner_problem(ptr_);

	//===== Third variant =====\\  
	call_class_ptr->call_inner_problem();				// calls inner problem


	//===== Old variant =====\\
	//options_pso opts_pso;
	//opts_pso.swarm_size = 10;
	//opts_pso.constraints = { 0, 10 };
	//opts_pso.max_Iter = 100;
	//opts_pso.omega_end = opts_pso.omega_start = 0.7;
	//opts_pso.tolerance = 1e-08;
	//opts_pso.debug = 0;
	//opts_pso.seed = 100;
	//test_fun_powell_pso_3 Test_pso;  // Initalizes test_fun Klass as Test!
	//Test_pso.params_ = params_;
	//Test_pso.transfer_ = transfer_;
	//initGuess = { 4, 4 };
	//PSO Inner_Problem(&Test_pso, &initGuess, &opts_pso);
	//Inner_Problem.solve();
	// 
	//result = Inner_Problem.solution();

	//=======================================================\\

	result = transfer_->inner;

	(*fval)[0] = pow(((*params_).a_soll - (*transfer_).inner[0]), 2.0);
	(*fval)[1] = pow(((*params_).b_soll - (*transfer_).inner[1]), 2.0);

	//std::cout << "======== outer input values ========" << std::endl;
	//std::cout << "(*vec)[0] " << (*vec)[0] << " | (*vec)[1] " << (*vec)[1] << std::endl;
	//std::cout << "======== outer output values ========" << std::endl;
	//std::cout << "(*fval)[0] " << (*fval)[0] << " | (*fval)[1] " << (*fval)[1] << std::endl;
	//std::cin.get();
	return;
}

Benchmark_Functions::Benchmark_Functions() :
	params_{ nullptr },
	function_{ 0 },
	dimensions_{ 0 }
{
	problem_description = "=== Benchmark Problems ===";
}

Benchmark_Functions::Benchmark_Functions(std::string function_str) :
	params_{ nullptr },
	function_{ 0 },
	dimensions_{ 0 }
{
	if (function_str == "sphere") {
		problem_description = "=== Benchmark Sphere ===";
		function_ = 1;
	}
	else if (function_str == "rastrigin") {
		problem_description = "=== Benchmark Rastrigin ===";
		function_ = 2;
	}
	else if (function_str == "rosenbrock") {
		problem_description = "=== Benchmark Rosenbrock ===";
		function_ = 3;
	}
	else {
		std::cout << std::endl << "Choosen benchmark function not available!" << std::endl << std::endl;
	}
}

void Benchmark_Functions::zero_function(std::vector<double>* vec, std::vector<double>* fval){

	dimensions_ = (*vec).size();

	if (function_ == 1) {
		sphere(vec, fval);
	}
	else if (function_ == 2) {
		rastrigin(vec, fval);
	}
	else if (function_ == 3) {
		rosenbrock(vec, fval);
	}

	transfer_->inner = (*vec);
}

void Benchmark_Functions::set_function(std::string function_str) {

	if (function_str == "sphere") {
		problem_description = "=== Benchmark Sphere ===";
		function_ = 1;
	}
	else if (function_str == "rastrigin") {
		problem_description = "=== Benchmark Rastrigin ===";
		function_ = 2;
	}
	else if (function_str == "rosenbrock") {
		problem_description = "=== Benchmark Rosenbrock ===";
		function_ = 3;
	}
	else {
		std::cout << std::endl << "Choosen benchmark function not available!" << std::endl << std::endl;
	}
}

void Benchmark_Functions::set_function(int function_int) {

	if (1 <= function_int && function_int <= 3) {
		function_ = function_int;

		std::ostringstream oss;
		oss << "=== Benchmark Function Nr. " << function_ << " ===";
		problem_description = oss.str();

	}
	else {
		std::cout << std::endl << "Choosen benchmark function not available!" << std::endl << std::endl;
	}
}

void Benchmark_Functions::sphere(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at (0,0) without offset, boundaries ]-inf, inf[

	// Testcase with only one returned zero_function value
	double offset = 0.0;
	(*fval)[0] =  0.0 ;
	for (int d = 0; d < dimensions_; d++) {
		(*fval)[0] += pow((*vec)[d] - offset, 2);
	}
}

void Benchmark_Functions::rastrigin(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at (0,0) without offset, boundaries [-5.12, 5.12]

	double A = 10;
	double offset = 0.0;
	
	for (int d = 0; d < dimensions_; d++) {
		(*fval)[d] = A + (((*vec)[d]-offset) * ((*vec)[d] - offset) - A * std::cos(2 * pi * ((*vec)[d] - offset)));
	}
}

void Benchmark_Functions::rosenbrock(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at b^2, boundaries ]-inf, inf[

	double a = 100;
	double b = 1; 
	
	for (int d = 0; d < dimensions_-1; d++) {
		(*fval)[d] = a * (std::pow((*vec)[d + 1] * (*vec)[d + 1] -(*vec)[d], 2) + std::pow((*vec)[d]-b,2));
	}
	(*fval)[dimensions_ - 1] = std::pow(b - (*vec)[dimensions_ - 1], 2);
}
