// Particle Swarm Optimizer \\ 
//     ____  _____ ____ 
//    / __ \/ ___// __ \
//   / /_/ /\__ \/ / / /
//  / ____/___/ / /_/ / 
// /_/    /____/\____/    
// 
// @github: MichiMitti

#pragma once
#include "pso.hpp"

PSO::PSO() {
	getVal_funcCtr_ = 0;
}

PSO::PSO(Zero_Function* fun, std::vector<double>* initGuess, options_pso* opts) :
	fun_{ fun },
	initGuess_{ *initGuess },
	opts_{ *opts },	
	N_{ (int)initGuess->size() },
	constraints_{ opts->constraints },
	swarm_size_{ opts->swarm_size },
	omega_{ opts->omega_start },
	c1_{ opts->c1 },
	c2_{ opts->c2 },
	velocity_max_{ opts->velocity_max },
	velocity_max_vec{ },	
	Swarm{ },
	best_swarm_value{ DBL_MAX },
	best_swarm_position{ },
	old_best_swarm_value{ DBL_MAX },
	fvals{ },
	fval{ DBL_MAX },
	engine{ opts->seed },
	rand_dbl{ 0.0, 1.0 },
	iter{ 0 },
	i_break{ 0 },
	getVal_funcCtr_{ 1 },
	agent_file{ "" },
	delimiter{ ';' },
	state{ },
	groups{ 1 },
	diversity{ 0.0 },
	swarm_mean{ 0.0 },
	distance_to_mean{ 0.0 },
	iteration{ 0 },
	iter_last_improvement{ 0 },
	no_improvement{ 0.0 },
	Net{ nullptr },
	Net_1{ },
	Net_2{ },
	Net_3{ }
	{
	state.resize(N_, 0.0);
	fvals.resize(N_,0.0);	// Not ideal should be the size of outputdimensions of zero_function!!!

	opts_ptr = &opts_;

	setConstraints();

	agent_file = "lib/pso/" + opts_.filename;	// sets filepath to csv file
}

void PSO::getVals(Zero_Function* fun, std::vector<double>* initGuess, options_pso* opts) {
	//fun_ = fun;
	//initGuess_ = *initGuess;
	//opts_ = *opts;							// just opts_pointer 
	//N_ = (int)(*initGuess).size();
	//constraints_ = (*opts).constraints;
	//swarm_size_ = (*opts).swarm_size;
	//omega_ = (*opts).omega_start;
	//c1_ = (*opts).c1;
	//c2_ = (*opts).c2;
	//velocity_max_ = (*opts).velocity_max;
	//velocity_max_vec = { };
	//Swarm = { };
	//best_swarm_value = { DBL_MAX };
	//best_swarm_position = { };
	//fvals = { };
	//fval = { DBL_MAX };

	//engine.seed( opts->seed);
	//decltype(rand_dbl.param()) new_range(0, 1.0);	// I think it's better to set engine and distribution in default ctor
	//rand_dbl.param(new_range);

	//iter = 0;
	//getVal_funcCtr_ = 1 ;

	//fvals.resize(N_);

	//if (N_ != constraints_.size()) {
	//	// it sets all constraints to the first entry if not enought boundaries are given!!!
	//	std::cout << "Constraints are overriden in the PSO call number of dimensions in the constraints does not match the inital guesses!" << std::endl;
	//	for (int d = 1; d < N_; d++) {
	//		constraints_.push_back(constraints_[0]);
	//	}
	//}

	//opts_.constraints = constraints_; // sets constraints also to options_powell
}

void PSO::calfun(std::vector<double>* x, std::vector<double>* fval, options_pso* opts) { // returns the residual sum of squares of the function
	(*fun_).zero_function(x, fval);
	(*opts).cnt_calfun++;
	return;
}

void PSO::selectNet() {

	if (opts_ptr->rl_on == 1) {
		Net = &Net_1;
	}
	else if (opts_ptr->rl_on == 2) {
		Net = &Net_2;
	}
	else if (opts_ptr->rl_on == 3) {
		Net = &Net_3;
	}
	else {
		std::cout << "\nThe selected Network is not trained!!!\n" << std::endl;
	}
}

void PSO::solve() {

	if (opts_ptr->rl_on == 0) {
		solve_std();	// solves with the standart pso 
	}
	else {
		selectNet();
		solve_rl();		// solves with rl trained networks
	}
}

void PSO::solve_std() {

	// displays header
	if (opts_.debug == 1) {
		std::cout << "\n===========================================\n"
			<< "Particle Swarm Optimizer" <<
			"\n===========================================\n\n" <<
			"maximum iterations: " << opts_.max_Iter << std::endl << "requested tolerance: " << opts_.tolerance <<
			"\n\nrunning PSO..." << std::endl;
	}

	// creates max velocity vector
	for (int d = 0; d < N_; d++) {
		velocity_max_vec.push_back((constraints_[d][1] - constraints_[d][0]) * opts_.velocity_max);
	}

	// creates file
	if (opts_.write_CSV) {
		createCSV();
	}

	best_swarm_position.resize(N_);

	initializeSwarm();

	while (iter < opts_.max_Iter) {

		// prints to file
		if (opts_.write_CSV) {
			writeParticleToCSV();
		}

		// calculate weight of inertia
		omega_ = opts_.omega_start - (opts_.omega_start - opts_.omega_end) / (opts_.max_Iter);

		//calculate velocity for all particles
		for (int p = 0; p < swarm_size_; p++) {

			// moves particle
			updateParticle(p);

			// updates best position and value for the *swarm*
			if (Swarm[p].fvalue < best_swarm_value) {
				best_swarm_value = Swarm[p].fvalue;

				if (std::isnan(best_swarm_value)) {
					outflag_ = 3;
				}
				if (std::isinf(best_swarm_value)) {
					outflag_ = 4;
				}

				best_swarm_position = Swarm[p].position;
			}

		}

		// prints for debugging
		if (opts_.debug == 1) {
			if (iter % 10 == 0) {
				std::cout << "Iter(" << iter << "): ";
				for (int i = 0; i < best_swarm_position.size() - 1; i++) { std::cout << best_swarm_position[i] << " | "; }
				std::cout << best_swarm_position.back() << " || " << best_swarm_value << std::endl;
			}
		}

		// ends PSO if accurate enough
		if (best_swarm_value < opts_.tolerance) {
			iter++;

			outflag_ = 0;

			// prints found solution
			if (opts_.debug == 1) {
				std::cout << "Solution found within " << iter << " iterations!" << std::endl;
			}

			// close file
			if (opts_.write_CSV) {
				writeParticleToCSV();
				closeCSV();
			}

			break;
		}

		// Checks if the best swarm value changed the last iters
		if (std::fabs(old_best_swarm_value - best_swarm_value) < opts_.tolerance) {
			i_break++;
		}
		else {
			i_break = 0;
		}

		if (i_break == opts_.break_tol) {

			outflag_ = 10;
			done = 1;
			break;
		}


		// ends PSO if max function calls are reached
		if (opts_.maxFunEvals <= opts_.cnt_calfun) {
			outflag_ = 2;
			break;
		}

		old_best_swarm_value = best_swarm_value;
		iter++;
	}
	if (iter == opts_.max_Iter) {
		outflag_ = 1;
		std::cout << "No solution found within tolerance!" << std::endl;

		// close file
		if (opts_.write_CSV) {
			writeParticleToCSV();
			closeCSV();
		}
	}

	if (opts_.maxFunEvals <= opts_.cnt_calfun) {
		outflag_ = 2;
		std::cout << "No solution found within max function calls!" << std::endl;

		// close file
		if (opts_.write_CSV) {
			writeParticleToCSV();
			closeCSV();
		}
	}

	return;
}

void PSO::initializeSwarm() {

	i_break = 0;

	opts_.cnt_calfun = 0;

	best_swarm_value = DBL_MAX;

	Swarm.resize(swarm_size_);	// sets right size

	for (int p = 0; p < swarm_size_; p++) {

		// sets positions and velocities
		Swarm[p].position.resize(N_);
		Swarm[p].velocity.resize(N_);

		// initialize position randomly
		for (int d = 0; d < N_; d++) {
			if (opts_.starting_points[0].size() == 0) {
				Swarm[p].position[d] = ((constraints_[d][1] - constraints_[d][0]) * rand_dbl(engine) + constraints_[d][0]);
			}
			else {
				Swarm[p].position[d] = opts_.starting_points[p][d];
			}

			//random velocity for first iteration (Could decrease perfomance)
			Swarm[p].velocity[d] = ((constraints_[d][1] - constraints_[d][0]) * rand_dbl(engine) + constraints_[d][0]);
		}
		
		// calculates the zero_fuction falue and sets it for particle
		calfun(&Swarm[p].position, &fvals, &opts_);

		fval = std::accumulate(fvals.begin(), fvals.end(), 0.0);

		Swarm[p].fvalue = fval;

		Swarm[p].best_position = Swarm[p].position;
		Swarm[p].best_fvalue = Swarm[p].fvalue;

		// also sets bests swarm value
		if (best_swarm_value > Swarm[p].best_fvalue) {

			best_swarm_value = Swarm[p].best_fvalue;

			if (std::isnan(best_swarm_value)) {
				outflag_ = 3;
			}
			if (std::isinf(best_swarm_value)) {
				outflag_ = 4;
			}

			best_swarm_position = Swarm[p].position;
		}
	}
	// initialize velocity for each particle as random possible in upper for loop
	for (int p = 0; p < swarm_size_; p++) {
		for (int d = 0; d < N_; d++) {
			// calculates new *velocity* for current paritcle
			Swarm[p].velocity[d] = (
				(omega_ * Swarm[p].velocity[d]) +
				(rand_dbl(engine) * c1_ * (Swarm[p].best_position[d] - Swarm[p].position[d])) +
				(rand_dbl(engine) * c2_ * (best_swarm_position[d] - Swarm[p].position[d]))
				);

			// clips velocity if to big in positive direction
			if (Swarm[p].velocity[d] > velocity_max_vec[d]) {
				Swarm[p].velocity[d] = velocity_max_vec[d];
				if (opts_.debug == 1) { std::cout << "positive velocity cliped! of particle [" << p << "]" << std::endl; }
			}
			// clips velocity if to big in negative direction
			if (Swarm[p].velocity[d] < -velocity_max_vec[d]) {
				Swarm[p].velocity[d] = -velocity_max_vec[d];
				if (opts_.debug == 1) { std::cout << "negativ velocity cliped! of particle [" << p << "]" << std::endl; }
			}
		}
	}
}

void PSO::setDimensions(std::vector<std::vector<double>> constr, int dim) {
	constraints_ = constr;
	N_ = dim;
	fvals.resize(N_, 0.0);
	initGuess_.resize(N_, 0.0);
	setConstraints();
}

void PSO::setConstraints() {
	if (N_ != constraints_.size()) {
		// it sets all constraints to the first entry if not enought boundaries are given!!!
		//std::cout << "Constraints are overriden in the PSO call number of dimensions in the constraints does not match the inital guesses!" << std::endl;
		for (int d = 1; d < N_; d++) {
			constraints_.push_back(constraints_[0]);
		}
	}
	opts_.constraints = constraints_;			// sets constraints also to options_powell
}

void PSO::updateParticle(int p) {

	// updates Particle position
	calcPosition(p);

	// caluclates function value and sets it to particle
	calfun(&Swarm[p].position, &fvals, &opts_);

	fval = std::accumulate(fvals.begin(), fvals.end(), 0.0);
	Swarm[p].fvalue = fval;

	// updates best position and value for individual *particle*
	if (Swarm[p].fvalue < Swarm[p].best_fvalue) {
		Swarm[p].best_fvalue = Swarm[p].fvalue;

		Swarm[p].best_position = Swarm[p].position;
	}
}

void PSO::calcPosition(int p) {

	for (int d = 0; d < N_; d++) {

		// calculates new *velocity* for current paritcle
		Swarm[p].velocity[d] = (
			(omega_ * Swarm[p].velocity[d]) +
			(rand_dbl(engine) * c1_ * (Swarm[p].best_position[d] - Swarm[p].position[d])) +
			(rand_dbl(engine) * c2_ * (best_swarm_position[d] - Swarm[p].position[d]))
			);

		// clips velocity if to big in positive direction
		if (Swarm[p].velocity[d] > velocity_max_vec[d]) {
			Swarm[p].velocity[d] = velocity_max_vec[d];
			if (opts_.debug == 1) { std::cout << "positive velocity cliped! of particle [" << p << "]" << std::endl; }
		}
		// clips velocity if to big in negative direction
		if (Swarm[p].velocity[d] < -velocity_max_vec[d]) {
			Swarm[p].velocity[d] = -velocity_max_vec[d];
			if (opts_.debug == 1) { std::cout << "negativ velocity cliped! of particle [" << p << "]" << std::endl; }
		}

		// calculates new *position*
		Swarm[p].position[d] += Swarm[p].velocity[d];

		// set constratints to postition?
		// Dampening boundary condition
		if (Swarm[p].position[d] < opts_.constraints[d][0]) {	// first therm to go back to boundary second therm to dampen velocity
			Swarm[p].position[d] = Swarm[p].position[d] - ((Swarm[p].position[d] - opts_.constraints[d][0]) + rand_dbl(engine)*(Swarm[p].position[d] - opts_.constraints[d][0]));
		}
		if (Swarm[p].position[d] > opts_.constraints[d][1]) {	// first therm to go back to boundary second therm to dampen velocity
			Swarm[p].position[d] = Swarm[p].position[d] - ((Swarm[p].position[d] - opts_.constraints[d][1]) + rand_dbl(engine) * (Swarm[p].position[d] - opts_.constraints[d][1]));
		}
	}
}



std::vector<double> PSO::solution() {
	return best_swarm_position;
}

double PSO::solution_function() {
	return best_swarm_value;
}

int PSO::function_calls() {
	return (opts_).cnt_calfun;
}



void PSO::createCSV() {

	// stops time with high resolution
	start_time = std::chrono::high_resolution_clock::now(); 

	// to write date and time in csv file
	auto now = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(now);

	// removes old file and create new one
	try {
		if (std::filesystem::remove(agent_file)) {
			std::cout << "file " << opts_.filename << " deleted.\n";
			fout.open(agent_file, std::ios::out | std::ios::app);
			std::cout << "file " << opts_.filename << " created.\n";
		}
		else {
			fout.open(agent_file, std::ios::out | std::ios::app);
			std::cout << "file " << opts_.filename << " created.\n";
		}
	}
	catch (const std::filesystem::filesystem_error& err) {
		std::cout << "filesystem error: " << err.what() << '\n';
	}

	// writes date and time to csv
	fout << std::ctime(&end_time) /*<< delimiter*/;						// throws error if #define _CRT_SECURE_NO_WARNINGS is not defined (dont know why)

	// header for csv
	fout << "\n===========================================" << delimiter << "\n"
		<< "Particle Swarm Optimizer" << delimiter 
		<< "\n===========================================" << delimiter << "\n\n" 
		<< "maximum iterations: " << delimiter << opts_.max_Iter << delimiter << "\n"
		<< "requested tolerance: " << delimiter << opts_.tolerance << delimiter << "\n";
	
	fout << "Number of dimensions: " << delimiter << N_ << delimiter << "\n";
	fout << "Boundaries: " << delimiter << opts_.constraints[0][0] << delimiter << opts_.constraints[0][1] << delimiter << "\n";
	fout << "Swarm size: " << delimiter << swarm_size_ << delimiter << "\n";
		
	fout << "\n" << opts_.problem_description << delimiter << "\n" << std::endl;
	
	// variable names 
	fout << "Particle " << delimiter;
	for (int i = 0; i < N_; i++) {
		fout << "Position " << i + 1 << delimiter;
	}
	for (int i = 0; i < N_; i++) {
		fout << "Velocity " << i + 1 << delimiter;
	}
	fout << "Function value" << delimiter;

	fout << std::endl;
}

void PSO::writeParticleToCSV() {

	// writes iteration and best current swarm position
	fout << "Iteration: "<< delimiter << iter << delimiter << "Best: " << delimiter;
	for (int i = 0; i < N_; i++) {
		fout << best_swarm_position[i] << delimiter;
	}
	fout << best_swarm_value << delimiter;

	fout << std::endl;

	// writes position and velocity
	for (int p = 0; p < swarm_size_; p++) {

		fout << p << delimiter;

		for (int i = 0; i < N_; i++) {
			fout << Swarm[p].position[i] << delimiter;
		}
		for (int i = 0; i < N_; i++) {
			fout << Swarm[p].velocity[i] << delimiter;
		}
		fout << Swarm[p].fvalue << delimiter;

		fout << "\n";
	}

	fout << std::endl;
}

void PSO::closeCSV() {

	// end time for stopwatch
	end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds> (end_time - start_time);

	// prints if solution was found or not
	if (iter == opts_.max_Iter) {
		fout << "\n1" << delimiter << "!!!NO SOLUTION FOUND !!!" << delimiter << "\n";
	}
	else {
		fout << "\n0" << delimiter << "Problem solved!" << delimiter << "\n";
	}

	// writes performance
	fout << "\nSolved within: " << delimiter << iter << delimiter << " iterations" << delimiter << "\n";
	fout << "Solved with: " << delimiter << opts_.cnt_calfun << delimiter << " function calls" << delimiter << "\n";
	fout << "Solved after: " << delimiter << duration.count() << delimiter << " microseconds" << delimiter << std::endl;

	fout.close();
}

void PSO::solve_rl() {


	solve_rl_step(0);

	iter = 0;

	action = {				// set action matrix to right size
		{0.7, 1.5, 1.5 },
		{0.7, 1.5, 1.5 },
		{0.7, 1.5, 1.5 },
		{0.7, 1.5, 1.5 },
		{0.7, 1.5, 1.5 }
	};

	while (iter < opts_ptr->max_Iter) {

		getAction();

		solve_rl_step(iter);

		iter++;
		// terminates is optimizer is done
		if (done == 1) {
			break;
		}
	}
};

void PSO::solve_rl_step(int step) {
	iter = step;

	// initalisation step could also be done in a function which is called from RL
	if (iter == 0) {

		// first velocities are calculated with std hyperparameters
		opts_ptr->action = { { omega_, c1_, c2_ } };

		done = 0;

		// displays header
		if (opts_.debug == 1) {
			std::cout << "\n===========================================\n"
				<< "Particle Swarm Optimizer" <<
				"\n===========================================\n\n" <<
				"maximum iterations: " << opts_.max_Iter << std::endl << "requested tolerance: " << opts_.tolerance <<
				"\n\nrunning PSO..." << std::endl;
		}

		// creates max velocity vector
		for (int d = 0; d < N_; d++) {
			velocity_max_vec.push_back((constraints_[d][1] - constraints_[d][0]) * opts_.velocity_max);
		}

		// creates file
		if (opts_.write_CSV) {
			createCSV();
		}

		best_swarm_position.resize(N_, 0.0);
		swarm_mean.resize(N_, 0.0);
		iter_last_improvement = 0;

		initializeSwarm();
	}
	else {



		// prints to file
		if (opts_.write_CSV) {
			writeParticleToCSV();
		}

		// calculate weight of inertia
		//omega_ = opts_.omega_start - (opts_.omega_start - opts_.omega_end) / (opts_.max_Iter);

		//calculate velocity for all particles
		for (int p = 0; p < swarm_size_; p++) {

			// Sets hyperparams from options
			omega_ = opts_ptr->action[p % groups][0];
			c1_ = opts_ptr->action[p % groups][1];
			c2_ = opts_ptr->action[p % groups][2];

			// moves particle
			updateParticle(p);

			// updates best position and value for the *swarm*
			if (Swarm[p].fvalue < best_swarm_value) {
				best_swarm_value = Swarm[p].fvalue;

				if (std::isnan(best_swarm_value)) {
					outflag_ = 3;
				}
				if (std::isinf(best_swarm_value)) {
					outflag_ = 4;
				}

				best_swarm_position = Swarm[p].position;

				// updates this iteration for improvement
				iter_last_improvement = iter;
			}
		}

		// prints for debugging
		if (opts_.debug == 1) {
			if (iter % 10 == 0) {
				std::cout << "Iter(" << iter << "): ";
				for (int i = 0; i < best_swarm_position.size() - 1; i++) { std::cout << best_swarm_position[i] << " | "; }
				std::cout << best_swarm_position.back() << " || " << best_swarm_value << std::endl;
			}
		}

		if (best_swarm_value < opts_.tolerance) {

			outflag_ = 0;

			iter++;
			// prints found solution
			if (opts_.debug == 1) {
				std::cout << "Solution found within " << iter << " iterations!" << std::endl;
			}

			// close file
			if (opts_.write_CSV) {
				writeParticleToCSV();
				closeCSV();
			}

			done = 1;	// sets solver to done
		}

		// Checks if the best swarm value changed the last iters
		if (std::fabs(old_best_swarm_value - best_swarm_value) < opts_.tolerance ) {
			i_break++;
		}
		else {
			i_break = 0;
		}

		if (i_break == opts_.break_tol) {

			outflag_ = 10;
			done = 1;
		}
	}

	if (iter == opts_.max_Iter) {

		outflag_ = 1;

		std::cout << "No solution found within tolerance!" << std::endl;

		// close file
		if (opts_.write_CSV) {
			writeParticleToCSV();
			closeCSV();
		}
	}

	old_best_swarm_value = best_swarm_value;

	calcState();

	transformState();

	return;
}

void PSO::calcState() {

	//====== calcs iteration state ======\\ 
	iteration = iter / (double)opts_.max_Iter;
	

	//====== calcs diversity state ======\\ 
	// average position vector
	std::fill(swarm_mean.begin(),swarm_mean.end(),0.0);

	// calculates average position of all particles
	for (int d = 0; d < N_; d++) {
		for (int p = 0; p < swarm_size_; p++) {
			swarm_mean[d] += Swarm[p].position[d]/swarm_size_;
		}
	}

	// summs up distances between particle and mean
	diversity = 0.0;
	for (int p = 0; p < swarm_size_; p++) {
		distance_to_mean = 0.0;
		for (int d = 0; d < N_; d++) {
			distance_to_mean += std::pow(						// square not in paper but I think its necessary
				(Swarm[p].position[d] - swarm_mean[d]) / 	
				(constraints_[d][1] - constraints_[d][0]),		// scale to boundaries
				2);  
		}
		diversity += std::sqrt(distance_to_mean);
	}
	diversity = diversity / swarm_size_;
	
	//====== calcs no_improvement state ======\\ 
	no_improvement = (iter - iter_last_improvement) / (double)opts_.max_Iter;


	if (no_improvement < 0.0) {
		std::cout << "ERROR STATE 2" << std::endl;
	}

	state = { iteration, diversity, no_improvement };
	return;
}

void PSO::transformState() {

	transformed_state.resize(15, 0.0);

	for (int i = 0; i < state.size(); i++) { // 0-2
		for (int j = 0; j < 4; j++) { // 0-4
			transformed_state[i * 5 + j] = std::sin(state[i] * std::pow(2, j));
		}
	}
	return;
}

void PSO::getAction() {

	action_vector = Net->forward(transformed_state);

	convertAction();

	opts_ptr->action = action;
}

void PSO::convertAction() {
	for (int g = 0; g < 5; g++) {
		action[g][0] = action_vector[g * 3 + 0] * 0.5 + 0.5;			// omega
		action[g][1] = action_vector[g * 3 + 1] * 1.0 + 1.5;			// c1
		action[g][2] = action_vector[g * 3 + 2] * 1.0 + 1.5;			// c2
	}
}
