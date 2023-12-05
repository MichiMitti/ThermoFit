#include "init_call.hpp"


int  init_call_function(int case_selection) {

	parameters params;	// Parameters can be set for problem but not necessarily

//************************* Information ***************************\\ 

#pragma region Information
	if (case_selection == 0) {

 

		std::cout <<
			"  __________________________________________________________________\n" <<
			" |         ________                              __________         |\n" <<
			" |        /_  __/ /_  ___  _________ ___  ____  / ___ /_/ /_        |\n" <<
			" |         / / / __ \\/ _ \\/ ___/ __ `__ \\/ __ \\/ /_  / / __/        |\n" <<
			" |        / / / / / /  __/ /  / / / / / / /_/ / __/ / / /_          |\n" <<
			" |       /_/ /_/ /_/\\___/_/  /_/ /_/ /_/\\____/_/   /_/\\__/          |\n" <<
			" |          The Tool for Optimization in Thermodynamics             |\n" <<
			" |__________________________________________________________________|\n\n" << std::endl;

		std::cout <<
			"Welcome to ThermoFit!\n"
			"This is a example function shows how use ThermoFit.\n"
			"The most important class is the Call_Class wich calls and manages all problems and optimizers. \n"
			"The following examples show how this is done. Currently these testcases are implemented: \n\n"
			"1....Solution of a single problem (only Inner)\n"
			"2....Solution of multiple problems (two nested problems e.g. Inner and Outer)\n"
			"11...Benchmark functions\n"
			"21...Only a inital guess for Inner problem searched randomly\n"
			"22...Only a inital guess for Outer problem searched randomly\n"
			"23...A inital guess searched randomly for Inner and Outer probem\n"
			"101..Testing area\n";		
	}

#pragma endregion

//**************** Test case with only one problem ****************\\ 

#pragma region Single Problem

	else if (case_selection == 1) {
		std::vector<double> initGuess_Inner = {  }; // don't know why but {0.0,0.0} is not working for simplex
		std::vector<double> initGuess_Outer = { 1.0, 1.0 };

		// boundaries will be split automatically for Inner and Outer Problem
		std::vector<std::vector<double>> boundaries = {
			{0,10},		// Inner Problem x
			{0,11},		// Inner Problem y
		};

		// PSO Test
		options_pso opts_pso;							// pso options
		opts_pso.swarm_size = 10;
		opts_pso.max_Iter = 100;
		opts_pso.omega_end = opts_pso.omega_start = 0.7;
		opts_pso.tolerance = 1e-08;
		opts_pso.debug = 0;

		options_dhs opts_simplex;				// simplex options

		Call_Test_Inner Single_Problem;					// problem class (derived from Zero_Function)
		Single_Problem.params_ = &params;				// set parameters (also possible when problem needs *params passed in ctor)

		CallClass Call_test_Single(&Single_Problem, &boundaries);						// initialize Probelm, when only one boundary vector exits it can be passed in the ctor other can be set to nullptr or can be left empty

		// Are the solver options_powell the same as specified in the string?
		Call_test_Single.solve_Inner_problem_with("pso", {}, &initGuess_Inner);		// selected method and options_powell for this method must be compatible! It is not possible to catch this error (because a void pointer is passed)
		//Call_test_Single.set_inital_guess_options("pso", "pso", &opts_pso);			// pso can also be used for searching inital Guess with passing a NULL pointer the selected method can be overwritten or std options_powell are set 
		//Call_test_Single.set_inital_guess_options("pso");

		Call_test_Single.set_Outer_initial_guess(&initGuess_Outer);						// sets transfer_struct to right values! (used when Inner problem is dependent on Outer problems starting values)

		Call_test_Single.solve();														// solves problem

		Call_test_Single.display_results();												// displays result

	}

#pragma endregion

//************* Test case with Inner and Outer problem ************\\ 

#pragma region Multiple Problems
	else if (case_selection == 2) {

		std::vector<double> initGuess_Inner = { 1.0, 1.0 };				// (0.0,0.0) not working for Simplex!!!
		std::vector<double> initGuess_Outer = { 1.0, 1.0 };						// (0.0,0.0) not working for Simplex!!!

		// boundaries will be split automatically for Inner and Outer Problem
		std::vector<std::vector<double>> boundaries = {					// combined boundaries but can also be splitted into Inner and Outer boundaries
			{0,10},		// Inner Problem x
			{0,11},		// Inner Problem y
			{0,12},		// Outer Problem x
			{0,13},		// Outer Problem y
		}; 

		options_pso opts_pso;											// pso options
		opts_pso.swarm_size = 10;
		opts_pso.max_Iter = 100;
		opts_pso.omega_end = opts_pso.omega_start = 0.7;
		opts_pso.tolerance = 1e-08;
		opts_pso.debug = 0;
		opts_pso.rl_on = 3;

		options_dhs opts_simplex;								// simplex options
		opts_simplex.debug = 0;
		opts_simplex.log_iter = 0;
		opts_simplex.disp_iter = 100;

		Call_Test_Inner Inner_problem;									// Inner problem class (derived from Zero_Function)
		Call_Test_Outer Outer_problem;									// Outer problem class (derived from Zero_Function)

		Inner_problem.params_, Outer_problem.params_ = &params;			// set parameters (also possible when problem needs *params passed in ctor)

		CallClass Call_test_Multiple(&Inner_problem, &Outer_problem);

		// Are the solver options_powell the same as specified in the passed string?
		Call_test_Multiple.solve_Inner_problem_with("powell", &boundaries, &initGuess_Inner);							// selected method and options for this method must be compatible!
		Call_test_Multiple.solve_Outer_problem_with("pso_rl", &boundaries, &initGuess_Outer,&opts_pso);	// selected method and options_powell for this method must be compatible!
		//Call_test_Multiple.set_inital_guess_options("inner", "simplex", &opts_simplex);						// inner method when calc the outer guess can be set this way
		//Call_test_Multiple.set_inital_guess_options("pso");														// solves problem already completely ;)
		

		Call_test_Multiple.solve();										// solves problem
		Call_test_Multiple.display_results();							// displays result

	}
#pragma endregion

//********** Testing optimizer with Benchmark Functions ***********\\ 

#pragma  region Benchmark Funcions 
	else if (10 < case_selection && case_selection <= 19) {

	// boundaries will set for all dimensions when called with set_Boundary_condition(...)
	//std::vector<std::vector<double>> boundaries = { {-10.0, 10.0} };	// sphere
	std::vector<std::vector<double>> boundaries = { {-5.12, 5.12} };	// rastrigin
	//std::vector<std::vector<double>> boundaries = { {-10.0, 10.0} };	// rosenbrock

	std::vector<double> initGuess;
	//std::vector<double> initGuess = { 0.1, 0.1 ,0.1, 0.1}; // don't know why but {0.0,0.0} is not working for simplex

	std::default_random_engine engine{ SEED };
	std::uniform_real_distribution<double> random_dist{ boundaries[0][0], boundaries[0][1] }; // random inital guess value

	int dim = 2;									// problem dimensions
	double value = 0.2;								// inital guess value
	for (int i = 0; i < dim; i++) {
		//initGuess.push_back(value);				// fixed value
		initGuess.push_back(random_dist(engine));	//random value
	}

	Benchmark_Functions Benchmark;												// problem class (derived from Zero_Function)
	Benchmark.params_ = &params;												// set parameters (also possible when problem needs *params passed in ctor)

	// available functions: "sphere", "rastrigin", "rosenbrock"
	std::string used_function = "sphere";
	Benchmark.set_function(used_function);										// sets function (also possible with int for iteration if multipe should be tested)

	CallClass Call_Test_Optim(&Benchmark);										// initialize Probelm, when only one boundary vector exits it can be passed in the ctor other can be set to nullptr or can be left empty

	options_pso opts_pso;							// pso options
	opts_pso.swarm_size = 10;
	opts_pso.omega_start = opts_pso.omega_end = 0.7;
	opts_pso.tolerance = 1e-8;
	opts_pso.debug = 0;
	opts_pso.write_CSV = 0;
	opts_pso.problem_description = Benchmark.problem_description;
	opts_pso.filename = "log_" + used_function + ".csv";
	opts_pso.rl_on = 1;
	opts_pso.break_tol = 100;

	options_dhs opts_simplex;						// simplex options
	opts_simplex.tol_fun_ = 1e-8;
	opts_simplex.debug = 0;

	options_ga opts_ga;								// ga options
	opts_ga.population_size = 30;
	opts_ga.tolerance = 1e-02;
	opts_ga.debug = 0;
	opts_ga.write_CSV = 1;
	opts_ga.problem_description = Benchmark.problem_description;
	opts_ga.filename = "log_" + used_function + ".csv";

	// Are the solver options the same as specified in the string?
	Call_Test_Optim.solve_Inner_problem_with("pso_rl", &boundaries, &initGuess, &opts_pso);	// selected method and options_powell for this method must be compatible! It is not possible to catch this error (because a void pointer is passed)
	//Call_Test_Optim.solve_Inner_problem_with("ga", &boundaries, &initGuess, &opts_ga);
	//Call_Test_Optim.solve_Inner_problem_with("simplex", &boundaries, &initGuess, &opts_simplex);

	//Call_Test_Optim.set_Boundary_conditions(&boundaries);						// sets boundary conditions to all dimensions possible but not necessary it the boundaries are applied automatic to all dimension

	std::cout << Benchmark.problem_description << std::endl << std::endl;

	

	Call_Test_Optim.solve();													// solves problem

	Call_Test_Optim.display_results();											// displays result

	std::cout << "\nOutflag: " << Call_Test_Optim.outflag() << std::endl;		// displays outflag
	Call_Test_Optim.display_outflag();
	}
#pragma endregion

//*************** Only searches for an inital guess ***************\\ 

#pragma region Inital Guess Caluculation
	else if (20 < case_selection && case_selection <= 29) {

		std::cout << "Only inital guess will be caluculated" << std::endl;

		std::vector<std::vector<double>> boundaries_inner = {		// boundaries for Inner problem
			{0,10},		// Inner Problem x
			{0,11},	// Inner Problem y
		};
		std::vector<std::vector<double>> boundaries_outer = {		// boundaries for Outer problem

			{0,12},		// Outer Problem x
			{0,13}		// Outer Problem y
		};
		std::vector<double> initGuess_Inner = { 1.0, 1.0 };			// Inner inital guess
		std::vector<double> initGuess_Outer = { 2.0, 2.0 };			// Outer inital guess

		Call_Test_Inner Inner_problem;								// Inner problem class (derived from Zero_Function)
		Call_Test_Outer Outer_problem;								// Outer problem class (derived from Zero_Function)

		Inner_problem.params_, Outer_problem.params_ = &params;		// set parameters (also possible when problem needs *params passed in ctor)

	// Only a inital guess for Inner problem searched randomly
		if (case_selection == 21) {
			// Options which must be set when finding inital guess for Inner probelm
			CallClass Call_test_Initial_Guess(&Inner_problem);										// one problem must be specified
			Call_test_Initial_Guess.set_Outer_initial_guess(&initGuess_Outer);						// must be set if transfer_struct_->outer is used in Inner problem (for single problem)
			
			Call_test_Initial_Guess.random_guess_search(&boundaries_inner);
			Call_test_Initial_Guess.display_inital_guess_results();
		}
	// Only a inital guess for Outer problem searched randomly
		else if (case_selection == 22) {
			// Options which must be set when finding only inital guess for Outer probelm
			CallClass Call_test_Initial_Guess(&Inner_problem, &Outer_problem);						// two problems must be specified!
			Call_test_Initial_Guess.set_Inner_initial_guess(&initGuess_Inner);						// just for knowing number of dimensions!
			Call_test_Initial_Guess.set_Boundary_conditions(&boundaries_inner, "inner");			// boundaries for PSO (for inner problem)
			//Call_test_Initial_Guess.set_inital_guess_options("inner","simplex");					// 

			Call_test_Initial_Guess.random_guess_search(&boundaries_outer, "outer");
			Call_test_Initial_Guess.display_inital_guess_results();
		}
	// A inital guess searched randomly for Innerand Outer probem within boundaries
		else if (case_selection == 23) {
			// Options which must be set when finding inital guess for both probelms
			CallClass Call_test_Initial_Guess(&Inner_problem, &Outer_problem);						// two problems must be specified!

			Call_test_Initial_Guess.random_guess_search(&boundaries_inner, &boundaries_outer);
			Call_test_Initial_Guess.display_inital_guess_results();
		}
	}
#pragma  endregion

//************************* Testing Area **************************\\ 

#pragma  region Test Area 
	else if ((100 < case_selection) && (case_selection <= 109)) {

	// sigmoid transform function test
	if (case_selection == 101) {										// test for transformation

		std::cout << "This is a test for variable transformation!\n";

		std::vector<double> initGuess_Inner = { 1.0, 1.0 };				// (0.0,0.0) not working for Simplex!!!
		std::vector<double> initGuess_Outer = { 1.0, 1.0 };				// (0.0,0.0) not working for Simplex!!!

		std::vector<std::vector<double>> boundaries = {					// combined boundaries but can also be splitted into Inner and Outer boundaries
			{-10,10},		// Inner Problem x
			{-11,11},		// Inner Problem y
			{-12,12},		// Outer Problem x
			{-13,13},		// Outer Problem y
			//{0,14}			// Test if boundary error works but only if they are passed in the class call
		};

		Call_Test_Inner Inner_problem;									// Inner problem class (derived from Zero_Function)
		Call_Test_Outer Outer_problem;									// Outer problem class (derived from Zero_Function)

		Inner_problem.params_, Outer_problem.params_ = &params;			// set parameters (also possible when problem needs *params passed in ctor)

		CallClass Test_sigmoid(&Inner_problem, &Outer_problem, &boundaries);

		Test_sigmoid.solve_Inner_problem_with("pso", {}, & initGuess_Inner);
		Test_sigmoid.solve_Outer_problem_with("pso", {}, &initGuess_Outer);

		std::vector<double> vec_test = { -7.2, 1.2 };

		std::vector<double> vec_for;
		std::vector<double> vec_back;

		std::cout << "Input Inner: \n";
		for (int i = 0; i < vec_test.size(); i++) {
			std::cout << vec_test[i] << " | ";
		}
		std::cout << std::endl;

		vec_for = Test_sigmoid.sigmoid_forward_Inner(vec_test);

		std::cout << "Forward Inner: \n";
		for (int i = 0; i < vec_test.size(); i++) {
			std::cout << vec_for[i] << " | ";
		}
		std::cout << std::endl;

		vec_back = Test_sigmoid.sigmoid_backward_Inner(vec_for);

		std::cout << "Backward Inner: \n";
		for (int i = 0; i < vec_test.size(); i++) {
			std::cout << vec_back[i] << " | ";
		}
		std::cout << std::endl;
	}
}
#pragma endregion

	return 0;
}
