#include "training_functions.hpp"

Training_Function::Training_Function() :
	params_{ nullptr },
	function_{ 0 },
	dimensions_{ 0 },
	offset{ 10.0 },
	available_functions{ },
	sphere_domain{ {-1000.0, 1000.0} },
	rastrigin_domain{ { -5.12, 5.12} },
	rosenbrock_domain{ {-1000.0, 1000.0} },
	beale_domain{ {-4.5, 4.5} },
	three_camel_domain{ {-5.0, 5.0} },
	ackley_domain{ {-32.768, 32.768} },
	schwefel_domain{ {-500.0, 500.0} },
	levy_domain{ {-10.0, 10.0} },
	griewank_domain{ {-600, 600} },
	six_camel_domain{ {{-3.0, 3.0}, {-2.0, 2.0}} },
	alpine_domain{ {-10.0, 10.0} },
	dimensions{ },
	pi{ 3.14159265358979323846 },
	euler{ 2.7182818284590452 }
{
	problem_description = "=== Benchmark Problems ===";

	// constraints for 
	constraints_ = {
		sphere_domain,
		rastrigin_domain,
		rosenbrock_domain,
		beale_domain,
		three_camel_domain,
		ackley_domain,
		schwefel_domain,
		levy_domain,
		griewank_domain,
		six_camel_domain,
		alpine_domain
	};

	dimensions = {		// dimensions of function (0 means n dimensional!)
		{0, 0, 1},		// sphere		| rand dims			| use in training
		{0, 0, 1},		// rastrigin	| rand dims			| use in training
		{0, 0, 1},		// rosenbrock	| rand dims			| use in training
		{2, 1, 1},		// beale		| constr for all	| use in training
		{2, 1, 1},		// three_camel	| constr for all	| use in training
		{2, 1, 1},		// ackley		| constr for all	| use in training
		{0, 0, 1},		// schwefel 	| rand dims			| use in training
		{0, 0, 1},		// levy 		| rand dims			| use in training
		{0, 0, 1},		// griewank 	| rand dims			| use in training
		{2, 1, 1},		// six_camel 	| fixed constr		| use in training
		{0, 1, 0}		// alpine 		| constr for all	| not used
		//{x, 1, 1},	// nrtl			| fixed constr		| model
		//{x, 1, 2},	// pcsaft		| fixed constr		| model
	};

	if (constraints_.size() != dimensions.size()) {
		std::cout << "\nError in Training_Function! \n" <<
			"Dimensions of constraionts_(" << constraints_.size() << ") does not equal " <<
			"the dimensions of the dimensions(" << dimensions.size() << "). \n" << 
			"When implementing new functions this must also be updated" << std::endl;
	}
	else {
		available_functions = constraints_.size(); 
	}
}

void Training_Function::resetTrainingFunction(int number) {
	set_function(number);
	set_offset(0.0);
	dimensions_ = 2;
}


void Training_Function::zero_function(std::vector<double>* vec, std::vector<double>* fval) {

	dimensions_ = (*vec).size();

	if (function_ == 0) {
		sphere(vec, fval);
	}
	else if (function_ == 1) {
		rastrigin(vec, fval);
	}
	else if (function_ == 2) {
		rosenbrock(vec, fval);
	}
	else if (function_ == 3) {
		dimensions_ = 2;
		beale(vec, fval);
	}
	else if (function_ == 4) {
		dimensions_ = 2;
		three_camel(vec, fval);
	}
	else if (function_ == 5) {
		dimensions_ = 2;
		ackley(vec, fval);
	}
	else if (function_ == 6) {
		schwefel(vec, fval);
	}
	else if (function_ == 7) {
		levy(vec, fval);
	}
	else if (function_ == 8) {
		griewank(vec, fval);
	}
	else if (function_ == 9) {
		dimensions_ = 2;
		six_camel(vec, fval);
	}
	else if (function_ == 10) {
		alpine(vec, fval);
	}
	else {
		std::cerr << "Zero function not available!" << std::endl;
		std::cin.get();
		std::cin.get();
	}

	transfer_->inner = (*vec);
}

void Training_Function::set_function(std::string function_str) {

	if (function_str == "sphere") {
		problem_description = "=== Benchmark Sphere ===";
		function_ = 0;
	}
	else if (function_str == "rastrigin") {
		problem_description = "=== Benchmark Rastrigin ===";
		function_ = 1;
	}
	else if (function_str == "rosenbrock") {
		problem_description = "=== Benchmark Rosenbrock ===";
		function_ = 2;
	}
	else if (function_str == "beale") {
		problem_description = "=== Benchmark Beale ===";
		function_ = 3;
	}
	else if (function_str == "three_camel") {
		problem_description = "=== Benchmark Three-Hump Camel ===";
		function_ = 4;
	}
	else if (function_str == "ackley") {
		problem_description = "=== Benchmark Ackley ===";
		function_ = 5;
	}
	else if (function_str == "schwefel") {
		problem_description = "=== Benchmark Schwefel ===";
		function_ = 6;
	}
	else if (function_str == "levy") {
		problem_description = "=== Benchmark Levy ===";
		function_ = 7;
	}
	else if (function_str == "griewank") {
		problem_description = "=== Benchmark Griewank ===";
		function_ = 8;
	}
	else if (function_str == "six_camel") {
		problem_description = "=== Benchmark Six-Hump Camel ===";
		function_ = 9;
	}
	else if (function_str == "alpine") {
		problem_description = "=== Benchmark Alpine ===";
		function_ = 10;
	}
	else {
		std::cout << std::endl << "Choosen benchmark function not available!" << std::endl << std::endl;
	}
}

void Training_Function::set_function(int function_int) {

	if (0 <= function_int && function_int <= available_functions-1) {
		function_ = function_int;

		setDescription();
	}
	else {
		std::cout << std::endl << "Choosen benchmark function not available!" << std::endl << std::endl;
	}
}

void Training_Function::set_offset(double offset_dbl) {
	offset = offset_dbl;
}


void Training_Function::sphere(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at (0,0) without offset, boundaries ]-inf, inf[

	for (int d = 0; d < dimensions_; d++) {
		(*fval)[d] = pow((*vec)[d], 2) + offset/dimensions_;
	}
}

void Training_Function::rastrigin(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at (0,0) without offset, boundaries [-5.12, 5.12]

	for (int d = 0; d < dimensions_; d++) {
		(*fval)[d] = 10 + (((*vec)[d]) * ((*vec)[d]) - 10 * std::cos(2 * pi * ((*vec)[d]))) + offset/dimensions_;
	}
}

void Training_Function::rosenbrock(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at [1,...,1], boundaries ]-inf, inf[

	for (int d = 0; d < dimensions_ - 1; d++) {
		(*fval)[d] = 100 * (std::pow((*vec)[d + 1] - (*vec)[d] * (*vec)[d], 2) + std::pow(1-(*vec)[d], 2)) + offset/dimensions_;
	}
	(*fval)[dimensions_-1] = 100 * (std::pow((*vec)[0] - (*vec)[dimensions_ - 1] * (*vec)[dimensions_ - 1], 2) + std::pow(1 - (*vec)[dimensions_ - 1], 2)) + offset / dimensions_;
}

void Training_Function::beale(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at [3,0.5]

	std::fill(fval->begin(), fval->end(), 0.0);

	(*fval)[0] =	std::pow(1.5 - (*vec)[0] + (*vec)[0] * (*vec)[1],2) + 
					std::pow(2.25 - (*vec)[0] + (*vec)[0] * (*vec)[1] * (*vec)[1], 2) +
					std::pow(2.625 - (*vec)[0] + (*vec)[0] * (*vec)[1] * (*vec)[1] * (*vec)[1], 2) +
					offset;
}

void Training_Function::three_camel(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at [1,...,1], boundaries ]-inf, inf[

	std::fill(fval->begin(), fval->end(), 0.0);

	(*fval)[0] =	2 * std::pow((*vec)[0], 2) - 
					1.05 * std::pow((*vec)[0], 4) + 
					std::pow((*vec)[0], 6) / 6.0 +
					(*vec)[0] * (*vec)[1] + std::pow((*vec)[1], 2) + 
					offset;
}

void Training_Function::ackley(std::vector<double>* vec, std::vector<double>* fval) {

	std::fill(fval->begin(), fval->end(), 0.0);

	(*fval)[0] =	-20.0 * std::exp(-0.2 * std::sqrt(0.5 * (std::pow((*vec)[0], 2) + std::pow((*vec)[1], 2)))) -
					std::exp(0.5 * (std::cos(2 * pi * (*vec)[0]) + std::cos(2 * pi * (*vec)[1]))) + euler + 20.0 + 
					offset;
}

void Training_Function::schwefel(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at (0,0) without offset, boundaries [-5.12, 5.12]

	for (int d = 0; d < dimensions_; d++) {
		(*fval)[d] = 418.9829 - (*vec)[d] * std::sin(std::sqrt(std::abs((*vec)[d]))) + offset / dimensions_;
	}
}

void Training_Function::levy(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at (1, ... ,1) without offset, boundaries [-10, 10]

	double wi = 1.0 + ((*vec)[0] - 1.0) / 4.0;

	(*fval)[0] = std::pow(std::sin(pi * wi),2) + offset / dimensions_;

	double wd = 1.0 + ((*vec)[dimensions_ - 1.0] - 1.0) / 4.0;
	wd = std::pow(wd - 1.0,2) * (1.0 + std::pow(std::sin(2.0 * pi * wd),2)); 	// is always constant

	for (int d = 1; d < dimensions_; d++) { 				// fist dimension is already calculated
		wi = 1.0 + ((*vec)[d] - 1) / 4.0;
		
		(*fval)[d] = std::pow(wi - 1.0,2) * (1.0 + 10.0 * std::pow(std::sin(pi * wi + 1.0),2)) + wd + offset / dimensions_;
	}
}

void Training_Function::griewank(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at (0, ... ,0) without offset, boundaries [-600, 600]

	double product = 0.0;
	for (int d = 0; d < dimensions_; d++) {
		product *= std::cos((*vec)[d] / std::sqrt((double)d+1)); 			// hope the sqrt of an int works as aspected
	}

	product = product / dimensions_;			// plus changes to minus

	for (int d = 0; d < dimensions_; d++) {
		(*fval)[d] = std::pow((*vec)[d],2) / 4000.0 - product + offset / dimensions_;
	}
}

void Training_Function::six_camel(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at (0.0898, -0.7126) and (-0.0898, 0.7126) without offset, boundaries {[-3, 3], [-2,2]}


	std::fill(fval->begin(), fval->end(), 0.0);

	(*fval)[0] = (4.0 - 2.1 * (*vec)[0] * (*vec)[0] + std::pow((*vec)[0], 4) / 3.0) * (*vec)[0]*(*vec)[0] +
		(*vec)[0] * (*vec)[1] +
		(-4.0 + 4.0 * (*vec)[1]*(*vec)[1]) * (*vec)[1]*(*vec)[1] +
		1.0316 + offset;
}

void Training_Function::alpine(std::vector<double>* vec, std::vector<double>* fval) {
	// zero at (0,0) without offset, boundaries ]-10, 10[

	for (int d = 0; d < dimensions_; d++) {
		(*fval)[d] = std::abs((*vec)[d] * std::sin((*vec)[d]) + 0.1 * (*vec)[d]);
	}
}

void Training_Function::setDescription() {

	if (function_ == 0) {
		problem_description = "=== Benchmark Sphere ===";
	}
	else if (function_ == 1) {
		problem_description = "=== Benchmark Rastrigin ===";
	}
	else if (function_ == 2) {
		problem_description = "=== Benchmark Rosenbrock ===";
	}
	else if (function_ == 3) {
		problem_description = "=== Benchmark Beale ===";
	}
	else if (function_ == 4) {
		problem_description = "=== Benchmark Three-Hump Camel ===";
	}
	else if (function_ == 5) {
		problem_description = "=== Benchmark Ackley ===";
	}
	else if (function_ == 6) {
		problem_description = "=== Benchmark Schwefel ===";
	}
	else if (function_ == 7) {
		problem_description = "=== Benchmark Levy ===";
	}
	else if (function_ == 8) {
		problem_description = "=== Benchmark Griewank ===";
	}
	else if (function_ == 9) {
		problem_description = "=== Benchmark Six-Hump Camel ===";
	}
	else if (function_ == 10) {
		problem_description = "=== Benchmark Alpine ===";
	}
	else {
		std::cout << std::endl << "Choosen benchmark function not available!" << std::endl << std::endl;
	}
}

void Training_Function::displayDescription() {

	setDescription();

	std::cout << problem_description << std::endl;
}