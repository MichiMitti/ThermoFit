// Call function to call solver and optimizers
#pragma once
#include "call_function.hpp"

// Ctor for single Problem
CallClass::CallClass(Zero_Function* Fun_Inner, std::vector<std::vector<double>>* boundaries) :
	Fun_Inner_{ Fun_Inner },
	Fun_Outer_{ nullptr },
	transfer_struct_{ },
	multiple_problems{ 0 },
	problems_Specified{ 0 },
	boundary_{ boundaries },
	boundary_Inner{ },
	boundary_Outer{ },
	initGuess_Inner_calc{ },
	initGuess_Outer_calc{ },
	initGuess_Inner_{ &initGuess_Inner_calc },
	initGuess_Outer_{ &initGuess_Outer_calc },
	dimensions_Inner{ 0 },
	dimensions_Outer{ 0 },
	method_Inner_save_{ 0 },
	method_Outer_save_{ 0 },
	method_Inner_{ 0 },
	method_Outer_{ 0 },
	Problem_Inner_PSO{ nullptr },
	Problem_Outer_PSO{ nullptr },
	Problem_Inner_Powell{ nullptr },
	Problem_Outer_Powell{ nullptr },
	Problem_Inner_DHS{ nullptr },
	Problem_Outer_DHS{ nullptr },
	Problem_Inner_GA{ nullptr },
	Problem_Outer_GA{ nullptr },
	outflag_Inner{ 99 },
	outflag_Outer{ 99 },
	base_random_guesses{ 5 },
	best_guess_vec{ },
	best_guess_vec_2{ },
	method_Init_Guess_{ "random" },
	options_initalGuess_ptr_{ nullptr },
	opts_Inner_{ nullptr },
	opts_Outer_{ nullptr },
	options_pso_std_Inner{ },
	options_pso_std_Outer{ },
	options_powell_std{ },
	options_dhs_std{ },
	options_ga_std_Inner{ },
	options_ga_std_Outer{ },
	initGuess_case{ 0 },
	shift_for_std_boundary{ 0.0 },
	multiplier_for_std_boundary{ 1.0 },
	fval_Inner{ DBL_MAX },
	fval_Outer{ DBL_MAX },
	fun_eval_Inner{ 0 },
	fun_eval_Outer{ 0 },
	solveing_time{ 0.0 }
	{
	displayBanner();

	options_pso_Inner_ptr = &options_pso_std_Inner;
	options_pso_Outer_ptr = &options_pso_std_Inner;

	options_powell_Inner_ptr = &options_powell_std;
	options_powell_Outer_ptr = &options_powell_std;

	options_dhs_Inner_ptr = &options_dhs_std;
	options_dhs_Outer_ptr = &options_dhs_std;

	options_ga_Inner_ptr = &options_ga_std_Inner;
	options_ga_Outer_ptr = &options_ga_std_Inner;
}	

// Ctor for multiple Problem
CallClass::CallClass(Zero_Function* Fun_Inner, Zero_Function* Fun_Outer, std::vector<std::vector<double>>* boundaries) :
	Fun_Inner_{ Fun_Inner },
	Fun_Outer_{ Fun_Outer },
	transfer_struct_{ },
	multiple_problems{ 1 },
	problems_Specified{ 0 },
	boundary_{ boundaries },
	boundary_Inner{ },
	boundary_Outer{ },
	initGuess_Inner_calc{ },
	initGuess_Outer_calc{ },
	initGuess_Inner_{ &initGuess_Inner_calc },
	initGuess_Outer_{ &initGuess_Outer_calc },
	dimensions_Inner{ 0 },
	dimensions_Outer{ 0 },
	method_Inner_save_{ 0 },
	method_Outer_save_{ 0 },
	method_Inner_{ 0 },
	method_Outer_{ 0 },
	Problem_Inner_PSO{ nullptr },
	Problem_Outer_PSO{ nullptr },
	Problem_Inner_Powell{ nullptr },
	Problem_Outer_Powell{ nullptr },
	Problem_Inner_DHS{ nullptr },
	Problem_Outer_DHS{ nullptr },
	Problem_Inner_GA{ nullptr },
	Problem_Outer_GA{ nullptr },
	outflag_Inner{ 99 },
	outflag_Outer{ 99 },
	base_random_guesses{ 5 },
	best_guess_vec{ },
	best_guess_vec_2{ },
	method_Init_Guess_{ "random" },
	options_initalGuess_ptr_{ nullptr },
	opts_Inner_{ nullptr },
	opts_Outer_{ nullptr },
	options_pso_std_Inner{ },
	options_pso_std_Outer{ },
	options_powell_std{ },
	options_dhs_std{ },
	options_ga_std_Inner{ },
	options_ga_std_Outer{ },
	initGuess_case{ 0 },
	shift_for_std_boundary{ 0.0 },
	multiplier_for_std_boundary{ 1.0 },
	fval_Inner{ DBL_MAX },
	fval_Outer{ DBL_MAX },
	fun_eval_Inner{ 0 },
	fun_eval_Outer{ 0 },
	solveing_time{ 0.0 }
	{
	displayBanner();

	options_pso_Inner_ptr = &options_pso_std_Inner;
	options_pso_Outer_ptr = &options_pso_std_Inner;

	options_powell_Inner_ptr = &options_powell_std;
	options_powell_Outer_ptr = &options_powell_std;

	options_dhs_Inner_ptr = &options_dhs_std;
	options_dhs_Outer_ptr = &options_dhs_std;

	options_ga_Inner_ptr = &options_ga_std_Inner;
	options_ga_Outer_ptr = &options_ga_std_Inner;
}

CallClass::~CallClass() {
	if (Problem_Outer_PSO != nullptr) {
		delete Problem_Outer_PSO;
		Problem_Outer_PSO = nullptr;
	}	
	if (Problem_Inner_Powell != nullptr) {
		delete Problem_Inner_Powell;
		Problem_Inner_Powell = nullptr;
	}	
	if (Problem_Outer_Powell != nullptr) {
		delete Problem_Outer_Powell;
		Problem_Outer_Powell = nullptr;
	}	
	if (Problem_Inner_DHS != nullptr) {
		delete Problem_Inner_DHS;
		Problem_Inner_DHS = nullptr;
	}	
	if (Problem_Outer_DHS != nullptr) {
		delete Problem_Outer_DHS;
		Problem_Outer_DHS = nullptr;
	}
	if (Problem_Inner_GA != nullptr) {
		delete Problem_Inner_GA;
		Problem_Inner_GA = nullptr;
	}
	if (Problem_Outer_GA != nullptr) {
		delete Problem_Outer_GA;
		Problem_Outer_GA = nullptr;
	}

	if (Problem_Inner_PSO != nullptr) {
		delete Problem_Inner_PSO;
		Problem_Inner_PSO = nullptr;
	}
}

void CallClass::solve_Inner_problem_with(std::string method_Inner, std::vector<std::vector<double>>* boundaries, std::vector<double>* initGuess_Inner, void* opts_Inner) {
	// method: specifies solving method, boundary conditions, initial Guess, method options_powell
	// if no inital guess is provided the dimensions of the problem must be choosen

	// sets variables
	opts_Inner_ = opts_Inner;
	initGuess_Inner_ = initGuess_Inner;

	// Method selection section (options_powell are set)
	if (method_Inner == "pso") {
		method_Inner_save_ = 0;
		method_Inner_ = 0;
		if (opts_Inner_ == nullptr) {
			options_pso_Inner_ptr = &options_pso_std_Inner;
		}
		else {
			options_pso_Inner_ptr = (options_pso*)opts_Inner;
		}

		// checks if option type is correct
		check_options(method_Inner, options_pso_Inner_ptr->options_typ);

	}
	else if (method_Inner == "pso_rl") {
		method_Inner_save_ = 10;
		method_Inner_ = 10;
		if (opts_Inner_ == nullptr) {
			options_pso_Inner_ptr = &options_pso_std_Inner;
			options_pso_Inner_ptr->rl_on = 1;
		}
		else {
			options_pso_Inner_ptr = (options_pso*)opts_Inner;
			if (options_pso_Inner_ptr->rl_on < 1) { options_pso_Inner_ptr->rl_on = 1; }
		}
		// checks if option type is correct
		check_options(method_Inner, options_pso_Inner_ptr->options_typ);

	}
	else if (method_Inner == "powell") {
		method_Inner_save_ = 1;
		method_Inner_ = 1;
		if (opts_Inner_ == nullptr) {
			options_powell_Inner_ptr = &options_powell_std;
		}
		else {
			options_powell_Inner_ptr = (options_powell*)opts_Inner;
		}

		// checks if option type is correct
		//check_options(method_Inner, options_powell_Inner_ptr->options_typ);
	}
	else if (method_Inner == "simplex") {
		method_Inner_save_ = 2;
		method_Inner_ = 2;
		if (opts_Inner_ == nullptr) {
			options_dhs_Inner_ptr = &options_dhs_std;
		}
		else {
			options_dhs_Inner_ptr = (options_dhs*)opts_Inner;
		}
		// checks if option type is correct
		check_options(method_Inner, options_dhs_Inner_ptr->options_typ);
	}
	else if (method_Inner == "ga") {
		method_Inner_save_ = 3;
		method_Inner_ = 3;
		if (opts_Inner_ == nullptr) {
			options_ga_Inner_ptr = &options_ga_std_Inner;
		}
		else {
			options_ga_Inner_ptr = (options_ga*)opts_Inner;
		}

		// checks if option type is correct
		check_options(method_Inner, options_ga_Inner_ptr->options_typ);

	}
	else {
		std::cerr << "Please select right method in init_call! \nThe options are: \npowell \npso \nsimplex" << std::endl;
		std::cin.get();
	}

	// defines how many variabels must be solved (dimensions)
	if (initGuess_Inner_ == nullptr || (*initGuess_Inner_).empty()) { 
		
		initGuess_Inner_ = nullptr;															// if inital Guesses are not provided all are set to nullptr

		//if (boundaries != nullptr and (*boundaries).size() > 1) {							// if the size of the boundary condition are > 1 then they are used for problem dimensions
		//	std::cout << "Dimensions are taken from the Inner Boundaries!" << std::endl;
		//	dimensions_Inner = (*boundaries).size();
		//}
																				// if not the user is asked to provide a number
		std::cout << "You selected no initial guess for the Inner problem, this will start a Generator for finding it.\n"
			"Please choose the numbers of unknown variables: " << std::endl;
		std::cin >> dimensions_Inner;

	}
	else {
		dimensions_Inner = (*initGuess_Inner_).size();
		transfer_struct_.inner = (*initGuess_Inner_);
	}

	// sets boundaries
	set_Boundary_conditions(boundaries, "inner");
	std::cout << std::endl;

	problems_Specified++; // checks if Inner is called before Outer

	return;
}

void CallClass::solve_Outer_problem_with(std::string method_Outer, std::vector<std::vector<double>>* boundaries, std::vector<double>* initGuess_Outer, void* opts_Outer) {
	// method specifies method, boundary conditions, initial Guess, method options_powell
	// if no inital guess is provided the dimensions of the problem must be choosen

	// sets outer variabels
	initGuess_Outer_ = initGuess_Outer;
	opts_Outer_ = opts_Outer;

	// Method selection section (options_powell are set)
	if (method_Outer == "pso") {
		method_Outer_save_ = 0;
		method_Outer_ = 0;
		if (opts_Outer_ == nullptr) {
			options_pso_Outer_ptr = &options_pso_std_Outer;
		}
		else {
			options_pso_Outer_ptr = (options_pso*)opts_Outer;
		}

		// checks if option type is correct
		check_options(method_Outer, options_pso_Outer_ptr->options_typ);
	}
	else if (method_Outer == "pso_rl") {
		method_Outer_save_ = 10;
		method_Outer_ = 10;
		if (opts_Outer_ == nullptr) {
			options_pso_Outer_ptr = &options_pso_std_Outer;
			options_pso_Outer_ptr->rl_on = 1;
		}
		else {
			options_pso_Outer_ptr = (options_pso*)opts_Outer;
			if (options_pso_Outer_ptr->rl_on < 1) { options_pso_Outer_ptr->rl_on = 1; }
		}

		// checks if option type is correct
		check_options(method_Outer, options_pso_Outer_ptr->options_typ);
	}
	else if (method_Outer == "powell") {
		method_Outer_save_ = 1;
		method_Outer_ = 1;
	}

	// checks if option type is correct
	//if (options_powell_Inner_ptr->options_typ != method_Inner) {
	//	std::cerr << "Error! Selected Inner solver (" << method_Inner << ") and passed solver options_powell (" << options_powell_Inner_ptr->options_typ
	//		<< ") does not match!" << std::endl;
	//}
	else if (method_Outer == "simplex") {
		method_Outer_save_ = 2;
		method_Outer_ = 2;
		if (opts_Outer_ == nullptr) {
			options_dhs_Outer_ptr = &options_dhs_std;
		}
		else {
			options_dhs_Outer_ptr = (options_dhs*)opts_Outer;
		}

		// checks if option type is correct
		check_options(method_Outer, options_dhs_Outer_ptr->options_typ);

	}
	else if (method_Outer == "ga") {
		method_Outer_save_ = 3;
		method_Outer_ = 3;
		if (opts_Outer_ == nullptr) {
			options_ga_Outer_ptr = &options_ga_std_Outer;
		}
		else {
			options_ga_Outer_ptr = (options_ga*)opts_Outer;
		}

		// checks if option type is correct
		check_options(method_Outer, options_ga_Outer_ptr->options_typ);
	}
	else {
		std::cerr << "Please select right method in init_call! \nThe options are: \npowell \npso \nsimplex" << std::endl;
		std::cin.get();
		return;
	}

	// if no initial guess how many variabels must be solved? 
	if (initGuess_Outer_ == nullptr || (*initGuess_Outer_).empty()) {

		initGuess_Outer_ = nullptr;									// if inital Guesses are not provided all ptr are set to nullptr

		//if (boundaries != nullptr and (*boundaries).size() > 1) {	// if the size of the boundary condition are > 1 then they are used for problem dimensions
		//	std::cout << "Dimensions are taken from the Outer Boundaries!" << std::endl;
		//	dimensions_Outer = (*boundaries).size();
		//}
														// if not the user is asked to provide a number
		std::cout << "You selected no initial guess for the Outer problem, this will start a Generator for finding it.\n"
			"Please choose the numbers of unknown variables: " << std::endl;
		std::cin >> dimensions_Outer;

	}
	else {
		dimensions_Outer = initGuess_Outer_->size();
		transfer_struct_.outer = (*initGuess_Outer_);
	}

	// sets boundaries
	set_Boundary_conditions(boundaries, "outer");
	std::cout << std::endl;

	if (problems_Specified == 0) {
		std::cerr << "Please specifiy Inner problem first!" << std::endl;
		std::cin.get();
		std::cin.get();
	}
	else if (problems_Specified == 1) {
		problems_Specified++;
	}

	return;
}

void CallClass::solve() {
	// solves specified problem(s)

	// sets transfer struct and pointer to THIS class in problem
	(*Fun_Inner_).transfer_ = &transfer_struct_;
	(*Fun_Inner_).call_class_ptr = this;

	if (multiple_problems == 1) {
		// sets transfer_struct and pointer for Outer problem
		(*Fun_Outer_).transfer_ = &transfer_struct_;
		(*Fun_Outer_).call_class_ptr = this;
	}

	// calculates inital guess
	
	if (initGuess_Inner_ == nullptr || initGuess_Outer_ == nullptr) {

		// calculates inital guess for a sinlge problem and sets it to transfer struct
		if (multiple_problems == 0) {

			calculate_initial_guess();

			// sets Inner initial guess in the transfer_struct_
			transfer_struct_.inner = *(std::vector<double>*)initGuess_Inner_; // (warning dereferencing a nullptr ptr)
		}
		// calculates inital guess for multiple problem and sets it to transfer struct
		else if (multiple_problems == 1) {

			calculate_initial_guess();

			// sets Inner initial guess(es) in the transfer_struct_
			transfer_struct_.inner = *(std::vector<double>*)initGuess_Inner_; // (warning dereferencing a nullptr ptr)
			transfer_struct_.outer = *(std::vector<double>*)initGuess_Outer_; // (warning dereferencing a nullptr ptr)
		}
	}

	// restores options_powell from solve_..._problem_with() if they were changed by the random guess generator
	method_Inner_ = method_Inner_save_;
	if (multiple_problems == 1) {
		method_Outer_ = method_Outer_save_;
	}

// solver section for single problem
	if (multiple_problems == 0 && problems_Specified == 1) {		// if only one problem is specified this section will be executed

		std::cout << "Solving Inner problem with: " << method_to_string(method_Inner_) << "\n" << std::endl;

		time_start = std::chrono::high_resolution_clock::now();		// start stopwatch

		check_dimensions("inner");	// checks if all everything is set (also in transfer_struct_)

		// solves Outer problem
		solve_Inner();

		time_end = std::chrono::high_resolution_clock::now();		// end stopwatch and calc time passed
		duration = std::chrono::duration_cast<std::chrono::microseconds> (time_end - time_start);
		// Calculate solving_time in milliseconds
		solveing_time = static_cast<double>(duration.count()) / 1000.0;
	}

// solver section for multiple problems
	else if (multiple_problems == 1 && problems_Specified == 2) {		// if multiple problems are specified this section will be executed

		std::cout << "Solving Inner problem with: " << method_to_string(method_Inner_) << "\n";
		std::cout << "Solving Outer problem with: " << method_to_string(method_Outer_) << "\n" << std::endl;

		check_dimensions("inner");		// checks if all everything is set (also in transfer_struct_)
		check_dimensions("outer");		// checks if all everything is set (also in transfer_struct_)

		time_start = std::chrono::high_resolution_clock::now();		// start stopwatch

		// solves Outer problem
		solve_Outer();

		time_end = std::chrono::high_resolution_clock::now();		// end stopwatch and calc time passed
		duration = std::chrono::duration_cast<std::chrono::microseconds> (time_end - time_start);
		// Calculate solving_time in milliseconds
		solveing_time = static_cast<double>(duration.count()) / 1000.0;
	}
	else {
		std::cerr << "Please specify the problems with \"solve_Inner_problem_with(..)\" and \"solve_Outer_problem_with(...)\" in this order!" << std::endl;
		std::cin.get();
	}

	return;
}

void CallClass::solve_Inner() {

	// solves with PSO
	if (method_Inner_ == 0) {	// PSO

		// sets boundaries to inner PSO options
		options_pso_Inner_ptr->constraints = boundary_Inner;

		if (Problem_Inner_PSO != nullptr) {
			delete Problem_Inner_PSO;
		}

		// solves problem
		Problem_Inner_PSO = new PSO(Fun_Inner_, initGuess_Inner_, options_pso_Inner_ptr);

		Problem_Inner_PSO->solve();
		outflag_Inner = Problem_Inner_PSO->outflag();

		vec_Inner = Problem_Inner_PSO->solution();
		fval_Inner = Problem_Inner_PSO->solution_function();
		fun_eval_Inner = Problem_Inner_PSO->function_calls();

	}
	// solves with RL inforced PSO
	else if (method_Inner_ == 10) {	// PSO_RL

		// sets boundaries to inner PSO options
		options_pso_Inner_ptr->constraints = boundary_Inner;

		if (Problem_Inner_PSO != nullptr) {
			delete Problem_Inner_PSO;
		}

		if (options_pso_Inner_ptr->rl_on < 1) {
			options_pso_Inner_ptr->rl_on = 3;
		}

		// solves problem
		Problem_Inner_PSO = new PSO(Fun_Inner_, initGuess_Inner_, options_pso_Inner_ptr);

		Problem_Inner_PSO->solve();

		vec_Inner = Problem_Inner_PSO->solution();
		fval_Inner = Problem_Inner_PSO->solution_function();
		fun_eval_Inner = Problem_Inner_PSO->function_calls();
		outflag_Inner = Problem_Inner_PSO->outflag();

	}
	// solves with Powell
	else if (method_Inner_ == 1) {	// Powell

		if (Problem_Inner_Powell != nullptr) {
			delete Problem_Inner_Powell;
		}

		//powell.GetVals!!!;
		Problem_Inner_Powell = new powell(Fun_Inner_, initGuess_Inner_, options_powell_Inner_ptr);
		Problem_Inner_Powell->solve(); // needs probably initial Guess

		vec_Inner = Problem_Inner_Powell->solution();
		fval_Inner = Problem_Inner_Powell->solution_function();

		outflag_Inner = Problem_Inner_Powell->outflg();

		//if (Problem_Inner_Powell->outflg() == 1)
		//	std::cout << "Powell terminated: maximum number of Iterations reached." << std::endl;
		//if (Problem_Inner_Powell->outflg() == 2)
		//	std::cout << "Powell terminated: maximum number of function evaluations reached." << std::endl;
		//if (Problem_Inner_Powell->outflg() == 3)
		//	std::cout << "Powell terminated: 5 calls of calfun ineffective." << std::endl;
		//if (Problem_Inner_Powell->outflg() == 4)
		//	std::cout << "Powell terminated: failed to decrease Fsq_ using a new Jacobian." << std::endl;
		//if (Problem_Inner_Powell->outflg() == 5)
		//	std::cout << "Powell terminated: nearby stationary point predicted." << std::endl;

		//powell_ptr = (Zero_Function*)Fun_Inner_;
		//powell.GetVals!!!;
		//powell Problem_Inner_Powell(powell_ptr, &initGuess_Inner_, (options_powell*)opts_Inner_);
		//Problem_Inner_Powell.solve(); // needs probably initial Guess
		//result_Inner = Problem_Inner_Powell.solution();
		//if (Problem_Inner_Powell.outflg() == 1)
		//	std::cout << "Powell terminated: maximum number of Iterations reached." << std::endl;
		//if (Problem_Inner_Powell.outflg() == 2)
		//	std::cout << "Powell terminated: maximum number of function evaluations reached." << std::endl;
		//if (Problem_Inner_Powell.outflg() == 3)
		//	std::cout << "Powell terminated: 5 calls of calfun ineffective." << std::endl;
		//if (Problem_Inner_Powell.outflg() == 4)
		//	std::cout << "Powell terminated: failed to decrease Fsq_ using a new Jacobian." << std::endl;
		//if (Problem_Inner_Powell.outflg() == 5)
		//	std::cout << "Powell terminated: nearby stationary point predicted." << std::endl;

	}
	// solves with Simplex
	else if (method_Inner_ == 2) {	// Simplex

		if (Problem_Inner_DHS != nullptr) {
			delete Problem_Inner_DHS;
		}

		// solves Inner problem with Simplex (NO boundaries are set in the options_powell)
		Problem_Inner_DHS = new DHS(Fun_Inner_, options_dhs_Inner_ptr);
		vec_Inner = Problem_Inner_DHS->Simplex_Solve(initGuess_Inner_);
		fval_Inner = Problem_Inner_DHS->fval;
		
		fun_eval_Inner = Problem_Inner_DHS->function_calls();

		//outflag_Inner = Problem_Inner_DHS->outflag();
	}
	else if (method_Inner_ == 3) { // GA

		// sets boundaries to inner GA options_powell
		options_ga_Inner_ptr->constraints = boundary_Inner;

		if (Problem_Inner_GA != nullptr) {
			delete Problem_Inner_GA;
		}

		// solves problem
		Problem_Inner_GA = new GA(Fun_Inner_, initGuess_Inner_, options_ga_Inner_ptr);
		Problem_Inner_GA->solve();
		vec_Inner = Problem_Inner_GA->solution();
		fval_Inner = Problem_Inner_GA->solution_function();
		outflag_Inner = Problem_Inner_GA->outflag();
		fun_eval_Inner = Problem_Inner_GA->function_calls();
	}
	else {
		std::cerr << "Inner solver call was not sucsessful! \n See \"call_inner_problem()\" function" << std::endl;
		std::cin.get();
	}
}

void CallClass::solve_Outer() {

	// solves with PSO
	if (method_Outer_ == 0) {	// PSO

		// sets boundaries to Outer PSO options_powell
		options_pso_Outer_ptr->constraints = boundary_Outer;

		if (Problem_Inner_PSO != nullptr) {
			delete Problem_Outer_PSO;
		}

		// solves problem
		Problem_Outer_PSO = new PSO(Fun_Outer_, initGuess_Outer_, options_pso_Outer_ptr);
		Problem_Outer_PSO->solve();
		vec_Outer = Problem_Outer_PSO->solution();
		fval_Outer = Problem_Outer_PSO->solution_function();
		fun_eval_Outer = Problem_Outer_PSO->function_calls();
		outflag_Outer = Problem_Outer_PSO->outflag();


	}
	else if (method_Outer_ == 10) {

		std::cout << "RL PSO ON!" << std::endl; // Solves with RL enforcement

		// sets boundaries to outer PSO options
		options_pso_Outer_ptr->constraints = boundary_Outer;

		if (Problem_Inner_PSO != nullptr) {
			delete Problem_Outer_PSO;
		}

		if (options_pso_Outer_ptr->rl_on < 1) {
			options_pso_Outer_ptr->rl_on = 3;
		}

		// solves problem
		Problem_Outer_PSO = new PSO(Fun_Outer_, initGuess_Outer_, options_pso_Outer_ptr);

		Problem_Outer_PSO->solve();

		vec_Outer = Problem_Outer_PSO->solution();
		fval_Outer = Problem_Outer_PSO->solution_function();
		fun_eval_Outer = Problem_Outer_PSO->function_calls();
		outflag_Outer = Problem_Outer_PSO->outflag();
	}
	// solves with Powell
	else if (method_Outer_ == 1) {	// Powell

		if (Problem_Outer_Powell != nullptr) {
			delete Problem_Outer_Powell;
		}

		//powell.GetVals!!!;
		Problem_Outer_Powell = new powell(Fun_Outer_, initGuess_Outer_, options_powell_Outer_ptr);
		Problem_Outer_Powell->solve(); // needs probably initial Guess

		vec_Outer = Problem_Outer_Powell->solution();
		fval_Outer = Problem_Outer_Powell->solution_function();

		outflag_Outer = Problem_Outer_Powell->outflg();

		//powell_ptr = (Zero_Function*)Fun_Outer_;
		//powell.GetVals!!!;
		//powell Problem_Outer_Powell(powell_ptr, &initGuess_Outer_, (options_powell*)opts_Inner_);
		//Problem_Inner_Powell.solve(); // needs probably initial Guess
		//result_Inner = Problem_Outer_Powell.solution();
		//if (Problem_Outer_Powell.outflg() == 1)
		//	std::cout << "Powell terminated: maximum number of Iterations reached." << std::endl;
		//if (Problem_Outer_Powell.outflg() == 2)
		//	std::cout << "Powell terminated: maximum number of function evaluations reached." << std::endl;
		//if (Problem_Outer_Powell.outflg() == 3)
		//	std::cout << "Powell terminated: 5 calls of calfun ineffective." << std::endl;
		//if (Problem_Outer_Powell.outflg() == 4)
		//	std::cout << "Powell terminated: failed to decrease Fsq_ using a new Jacobian." << std::endl;
		//if (Problem_Outer_Powell.outflg() == 5)
		//	std::cout << "Powell terminated: nearby stationary point predicted." << std::endl;
		//std::cout << "x = " << result_Outer[0] << std::endl;
		//std::cout << "y = " << result_Outer[1] << std::endl;
	}
	// solves with Simplex
	else if (method_Outer_ == 2) {

		if (Problem_Outer_DHS != nullptr) {
			delete Problem_Outer_DHS;
		}

		// solves Inner problem with Simplex (NO boundaries are set in the options_powell)
		Problem_Outer_DHS = new DHS(Fun_Outer_, options_dhs_Outer_ptr);
		vec_Outer = Problem_Outer_DHS->Simplex_Solve(initGuess_Outer_);
		fval_Outer = Problem_Outer_DHS->fval;

		fun_eval_Outer = Problem_Outer_DHS->function_calls();
		//fun_eval_Outer = Problem_Outer_DHS->
		//outflag_Outer = Problem_Inner_Simplex->outflag();
	}
	// solves with GA
	else if (method_Outer_ == 3) {
		// sets boundaries to Outer GA options_powell
		options_ga_Outer_ptr->constraints = boundary_Outer;

		if (Problem_Outer_GA != nullptr) {
			delete Problem_Outer_GA;
		}

		// solves problem
		Problem_Outer_GA = new GA(Fun_Outer_, initGuess_Outer_, options_ga_Outer_ptr);
		Problem_Outer_GA->solve();
		vec_Outer = Problem_Outer_GA->solution();
		fval_Outer = Problem_Outer_GA->solution_function();
		fun_eval_Outer = Problem_Outer_GA->function_calls();
		outflag_Outer = Problem_Inner_GA->outflag();
	}
	else {
		std::cerr << "Someting went wrong while selecting the method for multiple problems!" << std::endl;
		std::cin.get();
	}

}


void CallClass::call_inner_problem() {
	// calls inner problem

	// monitor the inner and outer pointers
	//std::cout << " ptr Inner: " << ptr_to_class->Fun_Inner_ << std::endl;
	//std::cout << " ptr Outer: " << ptr_to_class->Fun_Outer_ << std::endl;

	// necessary when useing random inital guess search when both initial guesses are not provided
	if (method_Init_Guess_ == "random_method_no_inner_solution_necessary!") {
		return;
	}
	
	// solves Inner Problem
	solve_Inner();
}



void CallClass::calculate_initial_guess() {
	// selects how inital guess will be calculated

	if (method_Init_Guess_ == "pso") {	
		calc_init_guess_pso();				// PSO
	}
	else {
		calc_init_guess_random();			// random
	}

	display_inital_guess_results();

	return;
}

void CallClass::calc_init_guess_pso() {

	// calculation for inner problem
	if (initGuess_Inner_ == nullptr) {

		// set everything to transfer_struct and class and checks it
		set_transfer_struct("inner");
		check_dimensions("inner");
		if (multiple_problems == 1) {
			check_dimensions("outer");
		}

		initGuess_case = 1;

		// sets boundaries to PSO options_powell
		(*(options_pso*)options_initalGuess_ptr_).constraints = boundary_Inner;

		// sets PSO problem 
		PSO Find_Init_Guess(Fun_Inner_, initGuess_Inner_, (options_pso*)options_initalGuess_ptr_);
		// solves with PSO
		Find_Init_Guess.solve();

		fval_Inner = Find_Init_Guess.solution_function();

		// sets inital guess to solution found
		*initGuess_Inner_ = Find_Init_Guess.solution();

	}
	// calculation for outer probelm
	else if (initGuess_Outer_ == nullptr) {

		initGuess_case = 2;

		// set everything to transfer_struct and class and checks it
		set_transfer_struct("outer");
		check_dimensions("inner");
		check_dimensions("outer");
		
		// if no options_powell are selected std options_powell will be used
		if (options_initalGuess_ptr_ == nullptr) {
			std::cout << "PSO standart options will be used!" << std::endl << std::endl;
			options_initalGuess_ptr_ = &options_pso_std_Outer;
		}
		else {
			std::cout << "Set PSO options will be used!" << std::endl << std::endl;
		}

		// sets transfer_struct for inialal guess calcualtion
		//set_transfer_struct_Outer();

		// sets boundaries to PSO options_powell
		(*(options_pso*)options_initalGuess_ptr_).constraints = boundary_Outer;

		// sets PSO problem 
		PSO Find_Init_Guess(Fun_Outer_, initGuess_Outer_, (options_pso*)options_initalGuess_ptr_);
		// solves with PSO
		Find_Init_Guess.solve();

		fval_Outer = Find_Init_Guess.solution_function();

		// sets inital guess to solution found
		*initGuess_Outer_ = Find_Init_Guess.solution();

	}
	else {
		std::cerr << "No inital guess calculated! (problem in \"calc_init_guess_pso\")";
		std::cin.get();
		std::cin.get();
	}
}

void CallClass::calc_init_guess_random() {

	// if both guesses are not set
	if (initGuess_Inner_ == nullptr && initGuess_Outer_ == nullptr) {

		random_guess_search(&boundary_Inner, &boundary_Outer);	// calculates inital guess random for both
	}
	// calculation for inner problem
	else if (initGuess_Inner_ == nullptr) {

		random_guess_search(&boundary_Inner, "inner");			// calculates inital guess random for Inner
	}
	// calculation for outer probelm
	else if (initGuess_Outer_ == nullptr) {

		random_guess_search(&boundary_Outer, "outer");			// calculates inital guess random for Outer
	}

	else {
		std::cerr << "No inital guess calculated! (problem in \"calc_init_guess_random\")";
		std::cin.get();
		std::cin.get();
	}

}

void CallClass::random_guess_search(std::vector<std::vector<double>>* boundaries, std::string problem) {

	int dimensions = (*boundaries)[0].size();		// uses size of boundaries so it can be used for Inner and Outer problem

	if (problem == "outer") {
		dimensions_Outer = dimensions;
		initGuess_Outer_ = &initGuess_Outer_calc;
		boundary_Outer = *boundaries;

		initGuess_case = 2;								//just to show which Inner method was used (not very good solution XD)
	}
	else {
		dimensions_Inner = dimensions;
		initGuess_Inner_ = &initGuess_Inner_calc;
		boundary_Inner = *boundaries;

		initGuess_case = 1;
	}

	(*Fun_Inner_).transfer_ = &transfer_struct_;
	(*Fun_Inner_).call_class_ptr = this;

	if (multiple_problems == 1) {
		(*Fun_Outer_).transfer_ = &transfer_struct_;
		(*Fun_Outer_).call_class_ptr = this;
	}

	// prints out warning because a big calculation is performed
	std::cerr << std::endl;
	std::cerr << "Warning! Initial guess will be searched with random brute-force! \nTherefore !!! "
		<< pow(base_random_guesses, dimensions)
		<< " !!! cases will be tried out!" << std::endl;
	std::cerr << std::endl;

	std::cin.get();
	std::cin.get();

	// sets everthing to transfer_struct_
	if (multiple_problems == 0) {
		set_transfer_struct("inner");
	}
	else if (multiple_problems == 1) {
		if (problem == "outer") {
			set_transfer_struct("outer");

			if (problems_Specified == 0) {
				if (boundary_Inner.empty()) { std::cerr << "Error! Please specify Inner boundary conditions!" << std::endl << std::endl; }
				else { set_Boundary_conditions(&boundary_Inner, "inner"); }
			}
			if (options_initalGuess_ptr_ == nullptr) { options_pso_Inner_ptr = &options_pso_std_Inner; }
			else									 { options_pso_Inner_ptr = (options_pso*)options_initalGuess_ptr_; }
		}
		else {
			set_transfer_struct("inner");	// also does stuff for outer
		}
	}

	// checks if everthing is set correctly
	check_dimensions("inner");
	if (multiple_problems == 1) {
		check_dimensions("outer");
	}

	int number_of_guesses = pow(base_random_guesses, dimensions);			// number of guesses (possible combinations)

	std::vector<std::vector<double>> guesses(number_of_guesses, std::vector<double>(dimensions, 0.0)); // vector for guesses

	std::vector<double> results(dimensions, 0.0);							// results of zero_function
	double result;															// holds squared sum of results

	std::vector<double> best_vec(dimensions, 0.0);							// vector which produced best result
	double best_fval = DBL_MAX;												// best result so far


	std::default_random_engine engine{ SEED };										// random engine std::mt19937
	std::vector<std::uniform_real_distribution<double>> random_engine;		// vector with engines for problem (because boundaries must be set for each dimensions)

	// fills random_engine vector
	for (int d = 0; d < dimensions; d++) {
		random_engine.push_back(std::uniform_real_distribution<double>{ (*boundaries)[d][0], (*boundaries)[d][1] });
	}

	int iter = 0; // just couts iterations vor writing out some results

	// do for all guesses
	for (int swarm_size = 0; swarm_size < number_of_guesses; swarm_size++) {
		for (int d = 0; d < dimensions; d++) {
			guesses[swarm_size][d] = random_engine[d](engine);
		} // finished constructing one guess

		// calls zero_function
		if (problem == "outer") { // pointer would possibly be faster but this is clearer
			(*Fun_Outer_).zero_function(&guesses[swarm_size], &results);
		}
		else {
			(*Fun_Inner_).zero_function(&guesses[swarm_size], &results);
		}

		result = solution_conversion_n_to_1(&results);	// squares and summs up results

		// sets best function values and vectors
		if (result < best_fval) {
			best_vec = guesses[swarm_size];
			best_fval = result;

			// debugging best values
			/*std::cout << "New best values: ";
			for (int d = 0; d < dimensions; d++) {
				std::cout << best_vec[d] << "   \t";
			}
			std::cout << " || " << best_fval << std::endl;*/
		}

		// debugging guesses
		/*if (iter % 5 == 0) {
			for (int d = 0; d < dimensions; d++) {
				std::cout << guesses[swarm_size][d] << "   \t";
			}
			std::cout << " || " << result << std::endl;
		}*/
		iter++;

	}

	if (problem == "outer") {
		fval_Outer = best_fval;
		*initGuess_Outer_ = best_vec;
	}
	else {
		fval_Inner = best_fval;
		*initGuess_Inner_ = best_vec;
	}

	method_Init_Guess_ = "random";

	// sets inital Guess to best guess (best_vec)
	best_guess_vec = best_vec;

}

void CallClass::random_guess_search(std::vector<std::vector<double>>* boundaries_inner, std::vector<std::vector<double>>* boundaries_outer) {

	initGuess_case = 3;

	int dimensions_inner = (*boundaries_inner)[0].size();	// number of variables from boundaries
	int dimensions_outer = (*boundaries_outer)[0].size();	// number of variables from boundaries

	dimensions_Inner = dimensions_inner;					// set dimension they will be used later in set_boundary_conditions
	dimensions_Outer = dimensions_outer;					// set dimension they will be used later in set_boundary_conditions

	set_Boundary_conditions(boundaries_inner,"inner");		// sets boundaries also to transfer_struct_
	set_Boundary_conditions(boundaries_outer,"outer");		// sets boundaries also to transfer_struct_

	(*initGuess_Inner_).resize(dimensions_Inner);			// sets initial guess to right dimensions
	(*initGuess_Outer_).resize(dimensions_Outer);			// sets initial guess to right dimensions
	
	transfer_struct_.inner.resize(dimensions_Inner);		// sets transfer_struct_.inner to right dimensions
	transfer_struct_.outer.resize(dimensions_Outer);		// sets transfer_struct_.outer to right dimensions

	// prints out warning because a big calculation is performed
	std::cerr << std::endl;
	std::cerr << "Warning! Initial guess will be searched with random brute-force! \nTherefore !!! "
		<< pow(base_random_guesses, dimensions_inner) * pow(base_random_guesses, dimensions_outer)
		<< " !!! cases will be tried out!" << std::endl;
	std::cerr << std::endl;

	std::cin.get();
	std::cin.get();

	std::cout << "random searching for the best Inner and Outer inital guess!" << std::endl;

	// sets this string because Inner function should not be solved!
	set_inital_guess_options("random_method_no_inner_solution_necessary!");

	// sets transfer and ptr_ for calculation
	(*Fun_Inner_).transfer_ = &transfer_struct_;
	(*Fun_Outer_).transfer_ = &transfer_struct_;
	(*Fun_Inner_).call_class_ptr = this;
	(*Fun_Outer_).call_class_ptr = this;

	check_dimensions("inner");
	check_dimensions("outer");

	std::vector<std::vector<double>> inner_guesses(pow(base_random_guesses, dimensions_inner), std::vector<double>(dimensions_inner, 0.0)); // vector for inner guesses
	std::vector<std::vector<double>> outer_guesses(pow(base_random_guesses, dimensions_outer), std::vector<double>(dimensions_outer, 0.0)); // vector for outer guesses


	std::vector<double> inner_result_fvals(dimensions_inner, 0.0);				// results of inner zero_function
	std::vector<double> outer_result_fvals(dimensions_outer, 0.0);				// results of outer zero_function

	double inner_result_fval = 0.0;												// holds squared sum of inner_result_fvals
	double outer_result_fval = 0.0;												// holds squared sum of outer_result_fvals

	double inner_best_fval = DBL_MAX;											// best resut of inner problem summed up (for best case)
	double outer_best_fval = DBL_MAX;											// best resut of outer problem summed up (for best case)
	double best_fval = DBL_MAX;													// best resut of inner and outer problem summed up
	double current_fval = 0.0;													// result of inner_result_fval + outer_result_fval (NOT sqared!)

	std::vector<double> best_inner_vec;											// vector which produced best inner result
	std::vector<double> best_outer_vec;											// vector which produced best outer result

	std::default_random_engine engine{ SEED };									// random engine std::mt19937
	std::vector<std::uniform_real_distribution<double>> random_engine_inner;	// vector with engines for inner problem (because boundaries must be set for each dimensions)
	std::vector<std::uniform_real_distribution<double>> random_engine_outer;	// vector with engines for inner problem (because boundaries must be set for each dimensions)

	//int number_of_results = 10;				// number of results stored in the best_few_guesses_inner vector
	//std::vector<std::vector<double>> best_few_guesses_inner(number_of_results, std::vector<double>(3,0.0));
	//std::vector<std::vector<double>> best_few_guesses_outer(number_of_results, std::vector<double>(3, 0.0));
	//std::vector<double> best_few_results(number_of_results, 0.0); // NOT ZERO!!!

	// fills random_engine_inner
	for (int d_i = 0; d_i < dimensions_inner; d_i++) {
		random_engine_inner.push_back(std::uniform_real_distribution<double>{ (*boundaries_inner)[d_i][0], (*boundaries_inner)[d_i][1] });
	}

	// fills random_engine_outer
	for (int d_o = 0; d_o < dimensions_outer; d_o++) {
		random_engine_outer.push_back(std::uniform_real_distribution<double>{ (*boundaries_outer)[d_o][0], (*boundaries_outer)[d_o][1] });
	}

	int iter = 0; // just couts iterations for writing out some results

	// do for all outer guesses
	for (int n_o = 0; n_o < outer_guesses.size(); n_o++) {

		for (int d_o = 0; d_o < dimensions_outer; d_o++) {
			outer_guesses[n_o][d_o] = random_engine_inner[d_o](engine);
		} //finished constructing one outer guess

		transfer_struct_.outer = outer_guesses[n_o]; // sets transfer_struct.outer so it is set for inner and outer zero_function

		// do for all inner guesses
		for (int n_i = 0; n_i < outer_guesses.size(); n_i++) {
			for (int d_i = 0; d_i < dimensions_inner; d_i++) {
				inner_guesses[n_i][d_i] = random_engine_inner[d_i](engine);
			}//finished constructing one inner guess
			transfer_struct_.inner = inner_guesses[n_i]; // sets transfer_struct.inner so it is set for inner and outer zero_function

			// result of inner zero_function
			(*Fun_Inner_).zero_function(&inner_guesses[n_i], &inner_result_fvals);
			inner_result_fval = solution_conversion_n_to_1(&inner_result_fvals); // squares and summs up results

			// result of outer zero_function
			(*Fun_Outer_).zero_function(&outer_guesses[n_o], &outer_result_fvals);
			outer_result_fval = solution_conversion_n_to_1(&outer_result_fvals); // squares and summs up results

			current_fval = inner_result_fval + outer_result_fval; // also squared would be possible but both vals are already squared

			// sets best function values and vectors
			if (current_fval < best_fval) {
				best_inner_vec = inner_guesses[n_i];
				best_outer_vec = outer_guesses[n_o];
				inner_best_fval = inner_result_fval;
				outer_best_fval = outer_result_fval;
				best_fval = current_fval;

				// debugging best values
				/* std::cout << "New best values: ";
				for (int d_i = 0; d_i < dimensions_inner; d_i++) {
					std::cout << best_inner_vec[d_i] << " ";
				}
				std::cout << " | ";
				for (int d_o = 0; d_o < dimensions_outer; d_o++) {
					std::cout << best_outer_vec[d_o] << " ";
				}
				std::cout << " || " << best_fval << std::endl;
				*/
			}

			// debugging guesses
			/*if (iter % 5 == 0) {
				for (int d_i = 0; d_i < dimensions_inner; d_i++) {
					std::cout << inner_guesses[n_i][d_i] << "   \t";
				}
				std::cout << " | ";
				for (int d_o = 0; d_o < dimensions_outer; d_o++) {
					std::cout << outer_guesses[n_o][d_o] << "   \t";
				}
				std::cout << " || " << current_fval << std::endl;
			}*/

			iter++;
		}
	}

	// old version just tries random guesses not every combination
	/*for (int swarm_size = 0; swarm_size < base_random_guesses; swarm_size++) {
		std::cout << iter << ": ";
		for (int d_i = 0; d_i < dimensions_inner; d_i++) {
			inner_guess[d_i] = random_engine_inner[d_i](engine);
			std::cout << inner_guess[d_i] << "  \t";
		}
		transfer_struct_random_search.inner = inner_guess;
		std::cout << " | ";
		for (int d_o = 0; d_o < dimensions_outer; d_o++) {
			outer_guess[d_o] = random_engine_outer[d_o](engine);
			std::cout << outer_guess[d_o] << "  \t";
		}
		transfer_struct_random_search.outer = outer_guess;
		(*Fun_Inner_).zero_function(&inner_guess, &inner_result_fvals);
		inner_result_fval = solution_conversion_n_to_1(&inner_result_fvals);
		(*Fun_Outer_).zero_function(&outer_guess, &outer_result_fvals);
		outer_result_fval = solution_conversion_n_to_1(&outer_result_fvals);
		std::cout << " || " << inner_result_fval << " | " << outer_result_fval;
		std::cout << std::endl;
		current_fval = inner_result_fval + outer_result_fval;
		if (swarm_size == 0) {
			best_inner_vec = inner_guess;
			best_outer_vec = outer_guess;
			best_fval = current_fval;
		}
		else if (current_fval < best_fval) {
			best_inner_vec = inner_guess;
			best_outer_vec = outer_guess;
			best_fval = current_fval;
		}
		iter++;
	}*/

	// sets function values for output
	fval_Inner = inner_best_fval;
	fval_Outer = outer_best_fval;

	// sets initGuess for further calculation
	*initGuess_Inner_ = best_inner_vec;
	*initGuess_Outer_ = best_outer_vec;

	// sets_best_vectors_2 in Class so it can be used further in class
	best_guess_vec_2.clear();
	best_guess_vec_2.push_back(best_inner_vec);
	best_guess_vec_2.push_back(best_outer_vec);

	// sets back inital gues options_powell to random so Inner solution is activated again
	method_Init_Guess_ = "random";
}



void CallClass::set_Boundary_conditions(std::vector<std::vector<double>>* boundaries, std::string problem) {

	// sets standart boundaries if boundaries are empty
	if ((boundary_ == nullptr && boundaries == nullptr) || (boundaries != nullptr && (*boundaries).empty())) {
		std::cout << "Standard boundaries between [0,1] except shift or scaling is applied!" << std::endl << std::endl;

		if (problem == "inner") {
			boundary_Inner.clear();
			for (int d = 0; d < dimensions_Inner; d++) {
				boundary_Inner.push_back({ (0.0 + shift_for_std_boundary) * multiplier_for_std_boundary, (1.0 + shift_for_std_boundary) * multiplier_for_std_boundary });
			}
		}
		else if (problem == "outer") {
			boundary_Outer.clear();
			for (int d = 0; d < dimensions_Outer; d++) {
				boundary_Outer.push_back({ (0.0 + shift_for_std_boundary) * multiplier_for_std_boundary, (1.0 + shift_for_std_boundary) * multiplier_for_std_boundary });
			}
		}
		else {
			std::cout << "Please specify which boundary condition you want to set! the options are: \"inner\" or \"outer\"." << std::endl << std::endl;
			std::cin.get();
			std::cin.get();
		}
	}
	// sets passed boundaries
	else if (boundaries != nullptr) {
		if (problem == "inner") {

			// set boundaries to passed value in function call if size matches exactly
			if ((*boundaries).size() == dimensions_Inner) {
				boundary_Inner = *boundaries;
			}
			// uses one dimentional boundaries for all
			else if ((*boundaries).size() == 1) {
				std::cout << "Warning! Passed boundaries have only one dimesion this boundary will set to all Inner dimensions!" << std::endl;
				boundary_Inner.clear();
				for (int d = 0; d < dimensions_Inner; d++) {
					boundary_Inner.push_back((*boundaries)[0]);
				}
			}
			// set boundaries to passed value in function call if it is large enough (but only the First few are used!)
			else if ((*boundaries).size() >= dimensions_Inner) {
				boundary_Inner.clear();
				for (int i = 0; i < dimensions_Inner; i++) {
					boundary_Inner.push_back((*boundaries)[i]);
				}
				std::cout << "Warning! Passed boundaries for Inner problem are to long! The FIRST " << dimensions_Inner << " elements will be used!" << std::endl;
			}
			else {
				std::cerr << "Error! It was not able to set Inner boundaries to boundaries given in the Function Call!" << std::endl;
				std::cin.get();
				std::cin.get();
			}
		}
		else if (problem == "outer") {
			// set boundaries to passed value in function call if size matches exactly
			if ((*boundaries).size() == dimensions_Outer) {
				boundary_Outer = *boundaries;
			}
			else if ((*boundaries).size() == 1) {
				std::cout << "Warning! Passed boundaries have only one dimesion this boundary will set to all Outer dimensions!" << std::endl;
				boundary_Outer.clear();
				for (int d = 0; d < dimensions_Outer; d++) {
					boundary_Outer.push_back((*boundaries)[0]);
				}
			}
			// set boundaries to passed value in function call if it is large enough (but only the Last few are used!)
			else if ((*boundaries).size() >= dimensions_Outer) {
				for (int d = 0; d < dimensions_Outer; d++) {
					boundary_Outer.push_back((*boundaries)[(*boundaries).size() - dimensions_Outer + d]);
				}
				std::cout << "Warning! Passed boundaries for Outer problem are to long! The LAST " << dimensions_Outer << " elements will be used!" << std::endl;
			}
			else {
				std::cerr << "Error! It was not able to set Outer boundaries to bounaries given in the Function Call!" << std::endl;
				std::cin.get();
				std::cin.get();
			}
		}
		else {
			std::cerr << "Error! It was not able to set Outer boundaries to bounaries given in the Function Call!\n"
				"Please specify which boundary condition you want to set! The options are : \"inner\" or \"outer\"." << std::endl;
			std::cin.get();
			std::cin.get();
		}
	}
	// sets boundary to the one in the Class Call!
	else if (boundary_ != nullptr) {
		if (problem == "inner") {
			// set boundaries to passed value in class call if size matches exactly
			if ((*boundary_).size() == dimensions_Inner) {
				boundary_Inner = *boundary_;
			}
			// uses one dimentional boundaries for all
			else if ((*boundary_).size() == 1) {
				std::cout << "Warning! Class Call boundaries will be set in Inner Boundaries have only one dimesion this boundary will be set to all dimensions!" << std::endl;
				boundary_Inner.clear();
				for (int d = 0; d < dimensions_Inner; d++) {
					boundary_Inner.push_back((*boundary_)[0]);
				}
			}
			// set boundaries to passed value in function call if it is large enough (but only the First few are used!)
			else if ((*boundary_).size() >= dimensions_Inner) {
				boundary_Inner.clear();
				for (int d = 0; d < dimensions_Inner; d++) {
					boundary_Inner.push_back((*boundary_)[d]);
				}
				std::cout << "Boundaries from class call are used! The FIRST " << dimensions_Inner << " elements will be used for the Inner problem!" << std::endl;
			}
			// sets standart boundaries if boundaries are empty
			else if ((*boundary_).empty()) {
				std::cout << "Standard boundaries between [0,1] except shift or scaling is applied!" << std::endl << std::endl;
				if (problem == "inner") {
					boundary_Inner.clear();
					for (int d = 0; d < dimensions_Inner; d++) {
						boundary_Inner.push_back({ (0.0 + shift_for_std_boundary) * multiplier_for_std_boundary, (1.0 + shift_for_std_boundary) * multiplier_for_std_boundary });
					}
				}
			}
			else {
				std::cerr << "Error! It was not able to set Inner boundaries to bounaries given in the Class Call!" << std::endl;
				std::cin.get();
				std::cin.get();
			}
		}
		else if (problem == "outer") {
			//NOT boundaries to passed value in class call if size matches exactly because mix-up is likely to happen!
			//if ((*boundary_).size() == dimensions_Outer) {
			//	boundary_Outer = *boundary_;
			//}
			// uses one dimentional boundaries for all

			// set boundaries to passed values in class call if it is large enough for Inner and Outer (but only the Last few are used!)
			if ((*boundary_).size() == (dimensions_Inner + dimensions_Outer)) {
				for (int d = 0; d < dimensions_Outer; d++) {
					boundary_Outer.push_back((*boundary_)[(*boundary_).size() - dimensions_Outer + d]);
				}
				std::cout << "Boundaries from class call are used! The LAST " << dimensions_Outer << " elements will be used for the Outer problem!" << std::endl;
			}
			// uses one dimentional boundaries from class call for all not sure if this is good!
			else if ((*boundary_).size() == 1) {
				std::cout << "Warning! Class Call boundaries will be set in Outer Boundaries have only one dimesion this boundary will be set to all dimensions!" << std::endl;
				boundary_Outer.clear();
				for (int d = 0; d < dimensions_Outer; d++) {
					boundary_Outer.push_back((*boundary_)[0]);
				}
			}
			// sets standart boundaries if boundaries are empty
			else if ((*boundary_).empty()) {
				std::cout << "Standard boundaries between [0,1] except shift or scaling is applied!" << std::endl << std::endl;
				boundary_Outer.clear();
				for (int d = 0; d < dimensions_Outer; d++) {
					boundary_Outer.push_back({ (0.0 + shift_for_std_boundary) * multiplier_for_std_boundary, (1.0 + shift_for_std_boundary) * multiplier_for_std_boundary });
				}
			}
			else {
				std::cerr << "Error! It was not able to set Outer boundaries to bounaries given in the Class Call!" << std::endl;
				std::cin.get();
				std::cin.get();
			}
		}
		else {
			std::cerr << "Error! It was not able to set Outer boundaries to bounaries given in the Class Call!\n"
				"Please specify which boundary condition you want to set!the options are : \"inner\" or \"outer\"." << std::endl;
			std::cin.get();
			std::cin.get();
		}
	}
	else {

	}

	if (problem == "outer") {
		transfer_struct_.boundary_outer = boundary_Outer;
	}
	else {
		transfer_struct_.boundary_inner = boundary_Inner;
	}
}

void CallClass::set_transfer_struct(std::string problem) {

	if (problem == "outer") {
		// set everything to transfer_struct and class
		set_Boundary_conditions(&boundary_Outer, "outer");
		initGuess_Outer_ = &initGuess_Outer_calc;
		(*initGuess_Outer_).resize(dimensions_Outer);
		transfer_struct_.outer = (*initGuess_Outer_);
		check_dimensions("outer");

		set_Boundary_conditions(&boundary_Inner, "inner");
		transfer_struct_.inner = (*initGuess_Inner_);
	}
	else {
		// set everything to transfer_struct and class
		set_Boundary_conditions(&boundary_Inner, "inner");
		initGuess_Inner_ = &initGuess_Inner_calc;
		(*initGuess_Inner_).resize(dimensions_Inner);
		transfer_struct_.inner.resize(dimensions_Inner);

		if (multiple_problems == 1) {
			// set everything to transfer_struct and class and checks it
			transfer_struct_.outer = (*initGuess_Outer_);
			set_Boundary_conditions(&boundary_Outer, "outer");
		}
	}
}

void CallClass::set_transfer_struct_special_variables(std::vector<double> special_variables) {
	// copies special_variables to transfer struct pointer would also be possible
	transfer_struct_.special_variables = special_variables;
}

void CallClass::set_Inner_initial_guess(std::vector<double>* initGuess_Inner) {
	// sets initial guess for Outer problem when solveing only Inner problem

	initGuess_Inner_ = initGuess_Inner;				// sets inital guess
	transfer_struct_.inner = *initGuess_Inner_;		// sets it also to transfer_struct_
	dimensions_Inner = (*initGuess_Inner).size();	// also sets dimensions
}

void CallClass::set_Outer_initial_guess(std::vector<double>* initGuess_Outer) {
	// sets initial guess for Outer problem when solveing only Inner problem

	initGuess_Outer_ = initGuess_Outer;				// sets inital guess
	transfer_struct_.outer = *initGuess_Outer_;		// sets it also to transfer_struct_
	dimensions_Outer = (*initGuess_Outer).size();	// also sets dimensions
}

void CallClass::set_inital_guess_options(std::string method, std::string options_Init_Guess, void* options_initGuess) {


	// setter for inital guess options_powell
	if (method == "pso") {	// options_inialGuess must be of type options_pso !!!
		method_Init_Guess_ = method;
		// setter for PSO inital guess options_powell
		if (options_initGuess == nullptr) {
			options_initalGuess_ptr_ = &options_pso_std_Inner; // just always set to Inner std options_powell
		}
		else {
			options_initalGuess_ptr_ = options_initGuess;
		}

		check_options(method, ((options_pso*)options_initalGuess_ptr_)->options_typ);
	}

	// setter for inital guess options_powell for Inner problem
	if (method == "inner") {
		
		if (options_Init_Guess == "pso") {	// options_inialGuess must be of type options_pso !!!
			// setter for PSO inital guess options_powell
			method_Inner_ = 0;
			if (options_initGuess == nullptr) {
				options_pso_Inner_ptr = &options_pso_std_Inner; // just always set to Inner std options_powell
			}
			else {
				options_pso_Inner_ptr = (options_pso*)options_initGuess;
			}
			check_options(method, options_pso_Inner_ptr->options_typ);
		}
		else if (options_Init_Guess == "powell") {
			method_Inner_ = 1;
			if (options_initGuess == nullptr) {
				options_powell_Inner_ptr = &options_powell_std; // just always set to Inner std options_powell
			}
			else {
				options_powell_Inner_ptr = (options_powell*)options_initGuess;
			}
			//check_options(method, options_powell_Inner_ptr->options_typ);
		}
		else if (options_Init_Guess == "simplex") {
			method_Inner_ = 2;
			if (options_initGuess == nullptr) {
				options_dhs_Inner_ptr = &options_dhs_std; // just always set to Inner std options_powell
			}
			else {
				options_dhs_Inner_ptr = (options_dhs*)options_initGuess;
			}
			check_options(method, options_dhs_Inner_ptr->options_typ);
		}
		else if (options_Init_Guess == "ga") {	// options_inialGuess must be of type options_pa !!!
			// setter for PSO inital guess options_powell
			method_Inner_ = 3;
			if (options_initGuess == nullptr) {
				options_ga_Inner_ptr = &options_ga_std_Inner; // just always set to Inner std options_powell
			}
			else {
				options_ga_Inner_ptr = (options_ga*)options_initGuess;
			}
			check_options(method, options_ga_Inner_ptr->options_typ);
		}
	}

	// sets method to random
	if (method == "random") {
		method_Init_Guess_ = method;
	}

	// sets method so Inner problem won't be solved
	if (method == "random_method_no_inner_solution_necessary!") { // important for inner and outer random search
		method_Init_Guess_ = method;
	}

	return;
}

void CallClass::set_base_for_random_guess_search(int base) {
	// changes base so more or less guesses are executed
	base_random_guesses = base;
}



std::vector<double> CallClass::sigmoid_forward_Inner(std::vector<double> vec) {

	std::vector<double> vec_for(vec.size(), 0.0);

	for (int i = 0; i < vec.size(); i++) {
		// if both boundaries are "0.0" no transformation is performed
		if ((boundary_Inner[i][1] == 0.0) && (boundary_Inner[i][0] == 0.0)) {		// == dangerous! double is not able to store certain values exactly "-" even more dangerous
			vec_for[i] = vec[i];
		}
		// transformation is performed (0...lower bound, 1...higher bound)
		else {
			vec_for[i] = -std::log(((boundary_Inner[i][1] - boundary_Inner[i][0]) / (vec[i] - boundary_Inner[i][0])) - 1.0);
		}
	}
	return vec_for;
}

std::vector<double> CallClass::sigmoid_forward_Outer(std::vector<double> vec) {

	std::vector<double> vec_for(vec.size(), 0.0);

	for (int i = 0; i < vec.size(); i++) {
		// if both boundaries are "0.0" no transformation is performed!
		if ((boundary_Outer[i][1] == 0.0) && (boundary_Outer[i][0] == 0.0)) {		// == dangerous! double is not able to store certain values exactly "-" even more dangerous
			vec_for[i] = vec[i];
		}
		// transformation is performed (0...lower bound, 1...higher bound)
		else {
			vec_for[i] = -std::log(((boundary_Outer[i][1] - boundary_Outer[i][0]) / (vec[i] - boundary_Outer[i][0])) - 1.0);
		}
	}
	return vec_for;
}

std::vector<double> CallClass::sigmoid_backward_Inner(std::vector<double> vec) {

	std::vector<double> vec_back(vec.size(), 0.0);

	for (int i = 0; i < vec.size(); i++) {
		// if both boundaries are "0.0" no transformation is performed
		if ((boundary_Inner[i][1] == 0.0) && (boundary_Inner[i][0] == 0.0)) {		// == dangerous! double is not able to store certain values exactly "-" even more dangerous
			vec_back[i] = vec[i];
		}
		// transformation is performed (0...lower bound, 1...higher bound)
		else {
			vec_back[i] = boundary_Inner[i][0] + (boundary_Inner[i][1] - boundary_Inner[i][0]) / (1.0 + std::exp(-vec[i]));
		}
	}
	return vec_back;
}

std::vector<double> CallClass::sigmoid_backward_Outer(std::vector<double> vec) {

	std::vector<double> vec_back(vec.size(), 0.0);

	for (int i = 0; i < vec.size(); i++) {
		// if both boundaries are "0.0" no transformation is performed
		if ((boundary_Outer[i][1] == 0.0) && (boundary_Outer[i][0] == 0.0)) {		// == dangerous! double is not able to store certain values exactly "-" even more dangerous
			vec_back[i] = vec[i];
		}
		// transformation is performed (0...lower bound, 1...higher bound)
		else {
			vec_back[i] = boundary_Outer[i][0] + (boundary_Outer[i][1] - boundary_Outer[i][0]) / (1.0 + std::exp(-vec[i]));
		}
	}
	return vec_back;
}

double CallClass::solution_conversion_n_to_1(std::vector<double>* vec) {
	// squares and sums up a vector

	double num = 0.0;

	for (int i = 0; i < (*vec).size(); i++) {
		num += (*vec)[i] * (*vec)[i];
	}
	return num;
}



void CallClass::display_results() {
	// displays results for calculation

	std::string unit = " [ms] ";
	double time = solveing_time;
	if (solveing_time > 1e3) {
		time = time / 1e3;
		unit = " [s] ";
	}
	else if (solveing_time > 3600 * 1e3) {
		time = time / (3600 * 1e3);
		unit = " [h] ";
	}

	// Outputs result for single problem
	if (multiple_problems == 0) {
		std::cout << std::endl << "===================== Results =====================" << std::endl;
		std::cout << "Results Inner problem (" << method_to_string(method_Inner_) << "): " << std::endl;
		for (int i = 0; i < vec_Inner.size() - 1; i++) { std::cout << vec_Inner[i] << " | "; }
		std::cout << vec_Inner.back() << " || " << fval_Inner << std::endl;
		std::cout << "FunEval: " << fun_eval_Inner << " | " << time << unit << std::endl;
		std::cout << "===================================================" << std::endl;
	}
	// Outputs result for multiple problems
	else if (multiple_problems == 1) {
		std::cout << std::endl << "===================== Results =====================" << std::endl;
		std::cout << "Results Inner problem (" << method_to_string(method_Inner_) << "): " << std::endl;
		for (int i = 0; i < vec_Inner.size() - 1; i++) { std::cout << vec_Inner[i] << " | "; }
		std::cout << vec_Inner.back() << " || " << fval_Inner << std::endl;

		std::cout << "Results Outer problem (" << method_to_string(method_Outer_) << "): " << std::endl;
		for (int i = 0; i < vec_Outer.size() - 1; i++) { std::cout << vec_Outer[i] << " | "; }
		std::cout << vec_Outer.back() << " || " << fval_Outer << std::endl;
		std::cout << "FunEval: " << fun_eval_Outer << " | " << time << unit << std::endl;
		std::cout << "===================================================" << std::endl;
	}

	std::cout << std::endl << std::endl;
	return;
}

std::vector<double> CallClass::result() {
	if (multiple_problems == 0) {
		return vec_Inner;
	}
	else if (multiple_problems == 1) {
		return vec_Outer;
	}
}

double CallClass::function_value() {
	if (multiple_problems == 0) {
		return fval_Inner;
	}
	else if (multiple_problems == 1) {
		return fval_Outer;
	}
}

int CallClass::function_evaluations() {
	if (multiple_problems == 0) {
		return fun_eval_Inner;
	}
	else if (multiple_problems == 1) {
		return fun_eval_Outer;
	}
}

void CallClass::display_inital_guess_results() {
	// displays results for Inital guess
	
	std::string method_for_inner_guess;					// holds string for solveing Inner method when outer is searched
	std::vector<double> inner_guess;					// holds solution of Inner problem if it was solved and Outer guess was unknown

	// decides if Inner problem was solved with a solver (nested probelem when Outer guess is searched)
	if (initGuess_case == 2) {
		method_for_inner_guess = method_to_string(method_Inner_);
		inner_guess = vec_Inner;
	}
	else {
		method_for_inner_guess = method_Init_Guess_;
		inner_guess = (*initGuess_Inner_);
	}

	// outputs single inital guess
	if (multiple_problems == 0) {
		std::cout << std::endl << "========= Results Init Guess =========" << std::endl;
		std::cout << "Inner inital guess (" << method_Init_Guess_ << "): " << std::endl;
		for (int i = 0; i < (*initGuess_Inner_).size() - 1; i++) { std::cout << (*initGuess_Inner_)[i] << " | "; }
		std::cout << (*initGuess_Inner_).back() << " || " << fval_Inner << std::endl;
	}
	// outputs result for inital guess
	else if (multiple_problems == 1) {
		std::cout << std::endl << "========= Results Init Guess =========" << std::endl;
		std::cout << "Inner inital guess (" << method_for_inner_guess << "): " << std::endl;
		for (int i = 0; i < inner_guess.size() - 1; i++) { std::cout << inner_guess[i] << " | "; }
		std::cout << inner_guess.back() << " || " << fval_Inner << std::endl;

		std::cout << "Outer inital guess (" << method_Init_Guess_ << "): " << std::endl;
		for (int i = 0; i < (*initGuess_Outer_).size() - 1; i++) { std::cout << (*initGuess_Outer_)[i] << " | "; }
		std::cout << (*initGuess_Outer_).back() << " || " << fval_Outer << std::endl;
	}

	std::cout << std::endl << "Press [enter] to continue calculation..." << std::endl;

	std::cin.get();
	std::cin.get();
	
	return;
}

int CallClass::outflag(std::string problem) {
	
	if (problem == "inner") {
		return outflag_Inner;
	}
	else if (problem == "outer") {
		return outflag_Outer;
	}
	else {
		std::cerr << "No outflag problem selected!" << std::endl;
	}
}

void CallClass::display_outflag(std::string problem) {

	if (problem == "inner") {
		if (outflag_Inner == 1)
		std::cout << "Inner Problem terminated: maximum number of Iterations reached." << std::endl;
		if (outflag_Inner == 2)
			std::cout << "Inner Problem terminated: maximum number of function evaluations reached." << std::endl;
		if (outflag_Inner == 3)
			std::cout << "Inner Problem: best function call was nan." << std::endl;
		if (outflag_Inner == 4)
			std::cout << "Inner Problem: best function call was inf." << std::endl;
		if (outflag_Inner == 10)
			std::cout << "Inner Problem finished: requested stability of the function value reached." << std::endl;
		if (outflag_Inner == 99)
			std::cout << "Inner outflag was not overidden!" << std::endl;
	}
	else if (problem == "outer") {
		if (outflag_Outer == 1)
			std::cout << "Outer Problem terminated: maximum number of Iterations reached." << std::endl;
		if (outflag_Outer == 2)
			std::cout << "Outer Problem terminated: maximum number of function evaluations reached." << std::endl;
		if (outflag_Outer == 3)
			std::cout << "Outer Problem: best function call was nan." << std::endl;
		if (outflag_Outer == 4)
			std::cout << "Outer Problem: best function call was inf." << std::endl;
		if (outflag_Outer == 10)
			std::cout << "Outer Problem finished: requested stability of the function value reached." << std::endl;
		if (outflag_Outer == 99)
			std::cout << "Outer outflag was not overidden!" << std::endl;
	}
	else {
		std::cerr << "No outflag problem selected!" << std::endl;
	}
}

void CallClass::check_dimensions(std::string problem) {
	// checks if everything is set correctly

	//checks Outer dimensions for all necessary variables
	if (problem == "outer") {
		if (dimensions_Outer != (*initGuess_Outer_).size()) {
			std::cerr << std::endl << "Error! Dimensions of Outer inital guess and dimensions_Outer and  does not match!" << std::endl;
		}
		if (dimensions_Outer != boundary_Outer.size()) {
			std::cerr << std::endl << "Error! Dimensions of Outer boundaris and dimensions_Outer and  does not match!" << std::endl;
		}
		if (dimensions_Outer != Fun_Outer_->transfer_->boundary_outer.size()) {
			std::cerr << std::endl << "Error! Dimensions of Outer zero_function (transfer_struct.boundary_outer) and dimensions_Outer and  does not match!" << std::endl;
		}
		if (dimensions_Outer != Fun_Outer_->transfer_->outer.size()) {
			std::cerr << std::endl << "Error! Dimensions of Outer zero_function (transfer_struct.outer) and dimensions_Outer and  does not match!" << std::endl;
		}
	}
	//checks Inner dimensions for all necessary variables
	else {
		if (dimensions_Inner != (*initGuess_Inner_).size()) {
			std::cerr << std::endl << "Error! Dimensions of Inner inital guess and dimensions_Inner and  does not match!" << std::endl;
		}
		if (dimensions_Inner != boundary_Inner.size()) {
			std::cerr << std::endl << "Error! Dimensions of Inner boundaris and dimensions_Inner and  does not match!" << std::endl;
		}
		if (dimensions_Inner != Fun_Inner_->transfer_->boundary_inner.size()) {
			std::cerr << std::endl << "Error! Dimensions of Inner zero_function (transfer_struct.boundary_inner) and dimensions_Inner and  does not match!" << std::endl;
		}
		if (dimensions_Inner != Fun_Inner_->transfer_->inner.size()) {
			std::cerr << std::endl << "Error! Dimensions of Inner zero_function (transfer_struct.inner) and dimensions_Inner and  does not match!" << std::endl;
		}
	}
}

void CallClass::check_options(std::string method, std::string options_) {

	// checks if option type is correct

	if (method != options_) {
		if (method != options_ + "_rl") {
			std::cerr << std::endl << "Error! Selected solver (" << method << ") and passed solver options (" << options_ + "_rl"
				<< ") does not match!" << std::endl << std::endl;
			std::cin.get();
			std::cin.get();
		}
	}
}

std::string CallClass::method_to_string(int method) {

	if (method == 0) {
		return "pso";
	}
	if (method == 10) {
		return "pso_rl";
	}
	else if (method == 1) {
		return "powell";
	}
	else if (method == 2) {
		return "simplex";
	}
	else if (method == 3) {
		return "ga";
	}
}

void CallClass::displayBanner(std::string banner) {

	if (banner == "ThermoFit") {
		std::cout <<
			"  __________________________________________________________________\n" <<
			" |         ________                              __________         |\n" <<
			" |        /_  __/ /_  ___  _________ ___  ____  / ___ /_/ /_        |\n" <<
			" |         / / / __ \\/ _ \\/ ___/ __ `__ \\/ __ \\/ /_  / / __/        |\n" <<
			" |        / / / / / /  __/ /  / / / / / / /_/ / __/ / / /_          |\n" <<
			" |       /_/ /_/ /_/\\___/_/  /_/ /_/ /_/\\____/_/   /_/\\__/          |\n" <<
			" |          The Tool for Optimization in Thermodynamics             |\n" <<
			" |__________________________________________________________________|\n\n" << std::endl;
	}

	if (banner == "OptiForge") {
		// old banner for OptiForge
		std::cout <<
			"  __________________________________________________________________\n" <<
			" |           ____        __  _______                                |\n" <<
			" |          / __ \\____  / /_/_/ ___/___  _________  __              |\n" <<
			" |         / / / / __ \\/ __/ / /_ / __ \\/ ___/ __ `/ _ \\            |\n" <<
			" |        / /_/ / /_/ / /_/ / __// /_/ / /  / /_/ /  __/            |\n" <<
			" |        \\____/ ____/\\__/_/_/   \\____/_/   \\__  /\\___/             |\n" <<
			" |            /_/                          /____/                   |\n" <<
			" |                The Tool for forging new Optimizers               |\n" <<
			" |__________________________________________________________________|\n\n" <<
			std::endl;
	}

	return;
}

void CallClass::printVectorToFile(std::vector<std::vector<double>>* matrix, std::string file_name, std::vector<std::vector<std::string>>* header, std::string extension, char delimiter) {

	std::fstream fout;									// Contains file stream (fileout)
	std::string name = file_name + extension;

	try {
		if (std::filesystem::remove(name)) {
			fout.open(name, std::ios::out | std::ios::app);
		}
		else {
			fout.open(name, std::ios::out | std::ios::app);
		}
	}
	catch (const std::filesystem::filesystem_error& err) {
		std::cout << "filesystem error: " << err.what() << '\n';
	}

	for (int i = 0; i < (*header).size(); i++) {
		for (int j = 0; j < (*header)[0].size(); j++) {
			fout << (*header)[i][j] << delimiter;
		}
		fout << "\n";
	}

	for (int i = 0; i < (*matrix).size(); i++) {
		for (int j = 0; j < (*matrix)[0].size(); j++) {
			fout << (*matrix)[i][j] << delimiter;
		}
		fout << "\n";
	}

	fout.close();
}

void CallClass::printVectorToFile(std::vector<double>* matrix, std::string file_name, std::vector<std::vector<std::string>>* header, std::string extension, char delimiter) {

	std::fstream fout;									// Contains file stream (fileout)
	std::string name = file_name + extension;

	try {
		if (std::filesystem::remove(name)) {
			fout.open(name, std::ios::out | std::ios::app);
		}
		else {
			fout.open(name, std::ios::out | std::ios::app);
		}
	}
	catch (const std::filesystem::filesystem_error& err) {
		std::cout << "filesystem error: " << err.what() << '\n';
	}

	for (int i = 0; i < (*header).size(); i++) {
		for (int j = 0; j < (*header)[0].size(); j++) {
			fout << (*header)[i][j] << delimiter;
		}
		fout << "\n";
	}

	for (int i = 0; i < (*matrix).size(); i++) {
		fout << (*matrix)[i] << delimiter;

		fout << "\n";
	}

	fout.close();
}
