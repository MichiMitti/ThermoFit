// Genetic Algorithm \\
//      ________ 
//    / ____/   |
//   / / __/ /| |
//  / /_/ / ___ |
//  \____/_/  |_|         
//
// @github: MichiMitti

#include "ga.hpp"

bool operator<(const Individual& Ind1, const Individual& Ind2) { return Ind1.fitness < Ind2.fitness; };

GA::GA(Zero_Function* fun, std::vector<double>* initGuess, options_ga* opts) :
	fun_{ fun },
	initGuess_{ *initGuess },
	opts_{ *opts },
	N_{ (int)initGuess->size() },
	constraints_{ opts->constraints },
	tolerance_{ opts->tolerance },
	population_size_{ opts->population_size },
	chromosome_lenght_{ N_ },
	max_generations_{ opts->max_Iter },
	mutation_rate_{ opts->mutation_rate },
	Population{ static_cast<const unsigned __int64>(population_size_) },
	Next_Gen{ static_cast<const unsigned __int64>(population_size_) },
	Combined_Pop{ static_cast<const unsigned __int64>(2*population_size_) },
	best_pop_fitness{ DBL_MAX },
	best_pop_chromosome{ },
	old_best_pop_fitness{ DBL_MAX },
	Parent1{ nullptr },
	Parent2{ nullptr },
	Offspring1{ },
	Offspring2{ },
	engine{ opts->seed },
	rand_dbl{ 0.0, 1.0 },
	rand_int_dim{0, chromosome_lenght_ - 1 },
	rand_int_pop{ 0, population_size_ - 1 },
	fvals{ },
	fval{ DBL_MAX },
	cum_prop_of_pop{ },
	iter{ 0 },
	i_break{ 0 }
	{ 

	opts_ptr = &opts_;

	Offspring1.chromosome.resize(chromosome_lenght_);
	Offspring1.chromosome.resize(chromosome_lenght_);

	fvals.resize(N_);

	cum_prop_of_pop.resize(population_size_);

	if (N_ != constraints_.size()) {
		// it sets all constraints to the first entry if not enought boundaries are given!!!
		std::cout << "Constraints are overriden in the PSO call number of dimensions in the constraints does not match the inital guesses!" << std::endl;
		for (int d = 1; d < N_; d++) {
			constraints_.push_back(constraints_[0]);
		}
	}

	opts_.constraints = constraints_;

	// fills random_engine vector
	for (int d = 0; d < N_; d++) {
		random_engine.push_back(std::uniform_real_distribution<double>{ constraints_[d][0], constraints_[d][1] });
	}	

	agent_file = "lib/ga/" + opts_.filename;	// sets filepath to csv file
}

void GA::calfun(std::vector<double>* x, std::vector<double>* fvals, options_ga* opts) { // returns the residual sum of squares of the function
	(*fun_).zero_function(x, fvals);
	(*opts).cnt_calfun++;
	return;
}

void GA::solve() {

	// displays header
	if (opts_.debug == 1) {
		std::cout << "\n===========================================\n"
			<< "Genetic Algorithm" <<
			"\n===========================================\n\n" <<
			"maximum iterations: " << opts_.max_Iter << std::endl << "requested tolerance: " << opts_.tolerance <<
			"\n\nrunning GA..." << std::endl;
	}

	// creates file
	if (opts_.write_CSV) {
		createCSV();
	}

	initializePopulation();
	evaluatePopulation();
	sortPopulation(&Population);

	best_pop_chromosome = Population[0].chromosome;
	best_pop_fitness = Population[0].fitness;

	while (iter < opts_.max_Iter) {

		// prints to file
		if (opts_.write_CSV) {
			writeIndividualToCSV();
		}

		if (opts_.debug == 1) {
			if (iter % 10 == 0) {
				std::cout << "Gen(" << iter << "): ";
				for (int i = 0; i < Population[0].chromosome.size() - 1; i++) { std::cout << Population[0].chromosome[i] << " | "; }
				std::cout << Population[0].chromosome.back() << " || " << Population[0].fitness << std::endl;
			}
		}

		if (Population[0].fitness < best_pop_fitness) {

			best_pop_fitness = Population[0].fitness;

			if (std::isnan(best_pop_fitness)) {
				outflag_ = 3;
			}
			if (std::isinf(best_pop_fitness)) {
				outflag_ = 4;
			}
			best_pop_chromosome = Population[0].chromosome;
		}


		if (best_pop_fitness < opts_.tolerance) {
			iter++;
			outflag_ = 0;
			// prints found solution
			if (opts_.debug == 1) {
				std::cout << "Solution found within " << iter << " iterations!" << std::endl;
			}

			// close file
			if (opts_.write_CSV) {
				writeIndividualToCSV();
				closeCSV();
			}

			if (opts_.rl_on) {
				done = 1;	// sets solver to done
			}
			break;
		}

		// Checks if the best swarm value changed the last iters
		if (std::fabs(old_best_pop_fitness - best_pop_fitness) < opts_.tolerance) {
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

		old_best_pop_fitness = best_pop_fitness;

		createNextGeneration();
		Population = Next_Gen;

		iter++;
	}

	if (iter == opts_.max_Iter) {
		outflag_ = 1;
		std::cout << "No solution found within tolerance!" << std::endl;

		// close file
		if (opts_.write_CSV) {
			writeIndividualToCSV();
			closeCSV();
		}
	}

	if (opts_.maxFunEvals <= opts_.cnt_calfun) {
		outflag_ = 2;
		std::cout << "No solution found within max function calls!" << std::endl;

		// close file
		if (opts_.write_CSV) {
			writeIndividualToCSV();
			closeCSV();
		}
	}

	return;
}


void GA::initializePopulation() {

	i_break = 0;
	best_pop_fitness = DBL_MAX;

	for (int p = 0; p < population_size_; p++) {						// loops over the population
		Population[p].chromosome.resize(chromosome_lenght_);
		for (int c = 0; c < chromosome_lenght_; c++) {					// loops over all dimensions (genes)

			Population[p].chromosome[c] = random_engine[c](engine);		// initalizes to random value
		}
	}
	return;
}

void GA::evaluatePopulation() {
	for (int p = 0; p < population_size_; p++) {

		calfun(&Population[p].chromosome, &fvals, &opts_);
		fval = std::accumulate(fvals.begin(), fvals.end(), 0.0);
		Population[p].fitness = fval;
	}
	return;
}

void GA::evaluateIndividual(Individual* Ind) {
	calfun(&(Ind->chromosome), &fvals, &opts_);
	fval = std::accumulate(fvals.begin(), fvals.end(), 0.0);
	Ind->fitness = fval;
}

void GA::sortPopulation(std::vector<Individual>* population) {
	std::sort(population->begin(), population->end());		// is certainly faster then own implementation
	return;
}


void GA::createNextGeneration() {

	// creats Probability vector once per generation
	createCumPropVec();	
	
	for (int i = 0; i < population_size_; i += 2) {

		// selects parents
		Parent1 = getParent();
		Parent2 = getParent();

		// do not allow two equal parents
		while (Parent1 == Parent2) {
			Parent2 = getParent();
		}

		// generats offsprings
		getOffspring(Parent1, Parent2);

		// addds offsprings to the Combined Population
		Combined_Pop[i] = Offspring1;
		Combined_Pop[i+1] = Offspring2;
		//Combined_Pop[i] = single_crossover(Parent1, Parent2);
		// Combined_Pop[i] = uniform_crossover(Parent1, Parent2);
	}

	// copies individuals
	std::copy(Population.begin(), Population.end(), Combined_Pop.begin() + population_size_);

	// sorts combined population
	sortPopulation(&Combined_Pop);

	// the best Individuals will survive
	for (int i = 0; i < population_size_; i++) {
		Next_Gen[i] = Combined_Pop[i];
	}

	return;
}


Individual* GA::getParent() {

	if (rand_dbl(engine) > 0.5) {
		// Tournament Selction
		return tournamentSelection();
	}
	else {
		// Biased random
		return biasedRandomSelection();
	}
}

Individual* GA::tournamentSelection() {
	Individual* Ind_1 = &(Population[rand_int_pop(engine)]);
	Individual* Ind_2 = &(Population[rand_int_pop(engine)]);

	while (Ind_1 == Ind_2) {
		Ind_2 = &(Population[rand_int_pop(engine)]);
	}
	if (Ind_1->fitness > Ind_2->fitness) {
		return Ind_1;
	}
	else {
		return Ind_2;
	}
}

void GA::createCumPropVec() {
	double sum = 0.0;

	for (int p = 0; p < population_size_; p++) {
		sum += Population[p].fitness;
	}

	std::vector<double> prop(population_size_);			// probabilities
	double prop_sum = 0.0;								// sum of propbabilies
	std::vector<double> norm_prop(population_size_);	// normalized probabilities

	for (int p = 0; p < population_size_; p++) {
		prop[p] = sum / Population[p].fitness;
		prop_sum += prop[p];
	}
	for (int p = 0; p < population_size_; p++) {
		norm_prop[p] = prop[p] / prop_sum;
	}

	double cum_total = 0.0;

	for (int p = 0; p < population_size_; p++) {
		cum_total += norm_prop[p];
		cum_prop_of_pop[p] = cum_total;
	}
}

Individual* GA::biasedRandomSelection() {

	double selected_value = rand_dbl(engine);

	for (int i = 0; i < population_size_; i++) {

		if (cum_prop_of_pop[i] >= selected_value) {
			return &Population[i];
		}
	}
	
	std::cerr << "Error in biasedRandomSelection()!" << std::endl;
}


void GA::getOffspring(Individual* A, Individual* B) {

	//Offspring1 = performCrossover_splitVec(A, B);
	//Offspring2 = performCrossover_splitVec(B, A);

	Offspring1 = doCrossover(A, B); // also other crossovers would be possible
	Offspring2 = doCrossover(B, A);

	mutate(&Offspring1);
	mutate(&Offspring2);

	evaluateIndividual(&Offspring1);
	evaluateIndividual(&Offspring2);

	return ;
}

Individual GA::doCrossover(Individual* Ind_1, Individual* Ind_2) {

	if (rand_dbl(engine) < 0.4) {
		// Single Point Crossover
		return singleCrossover(Ind_1, Ind_2);
	}
	else if(rand_dbl(engine) < 0.8) {
		// Uniform Crossover
		return uniformCrossover(Ind_1, Ind_2);
	}
	else {
		// Random Crossover (average with random numbers)
		return randomCrossover(Ind_1, Ind_2);
	}
}

Individual GA::randomCrossover(Individual* Ind_1, Individual* Ind_2) {
	// probably not so good because it averages over Parents therefore the solution gets in the middle of the boundary more likely I guess
	// but it will introduce new numbers into Generation

	Individual Offspring;
	Offspring.chromosome.resize(chromosome_lenght_);

	double random_weight;

	for (int i = 0; i < chromosome_lenght_; i++) {
		random_weight = rand_dbl(engine);
		Offspring.chromosome[i] = random_weight * Ind_1->chromosome[i] + (1 - random_weight) * Ind_2->chromosome[i];
	}

	return Offspring;
}

Individual GA::singleCrossover(Individual* Ind_1, Individual* Ind_2) {

	Individual Offspring;
	Offspring.chromosome.resize(chromosome_lenght_);

	int crossover_point = rand_int_dim(engine);

	for (int i = 0; i < crossover_point; i++) {
		Offspring.chromosome[i] = Ind_1->chromosome[i];
	}

	for (int i = crossover_point; i < chromosome_lenght_; i++) {

		Offspring.chromosome[i] = Ind_2->chromosome[i];				// look if any is already there
	}

	return Offspring;
}

Individual GA::uniformCrossover(Individual* Ind_1, Individual* Ind_2) {

	Individual Offspring;
	Offspring.chromosome.resize(chromosome_lenght_);

	double r;

	for (int c = 0; c < chromosome_lenght_; c++) {

		r = rand_dbl(engine);

		if (r < 0.5) {
			Offspring.chromosome[c] = Ind_1->chromosome[c];
		}
		else {
			Offspring.chromosome[c] = Ind_2->chromosome[c];
		}
	}
	return Offspring;
}

void GA::mutate(Individual* Ind) {
	for (int c = 0; c < chromosome_lenght_; c++) {
		if (rand_dbl(engine) < mutation_rate_) {
			Ind->chromosome[c] = constraints_[c][0] + rand_dbl(engine) * (constraints_[c][1] - constraints_[c][0]); // mutates randomly in boundaries
		}
	}
	
	return;
}


std::vector<double> GA::solution() {
	return Population[0].chromosome;
}

double GA::solution_function() {
	return Population[0].fitness;
}

int GA::function_calls() {
	return (opts_).cnt_calfun;
}



void GA::createCSV() {

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
		<< "Genetic Algorithm" << delimiter
		<< "\n===========================================" << delimiter << "\n\n"
		<< "maximum iterations: " << delimiter << opts_.max_Iter << delimiter << "\n"
		<< "requested tolerance: " << delimiter << opts_.tolerance << delimiter << "\n";

	fout << "Number of dimensions: " << delimiter << N_ << delimiter << "\n";
	fout << "Boundaries: " << delimiter << opts_.constraints[0][0] << delimiter << opts_.constraints[0][1] << delimiter << "\n";
	fout << "Population size: " << delimiter << population_size_ << delimiter << "\n";

	fout << "\n" << opts_.problem_description << delimiter << "\n" << std::endl;

	// variable names 
	fout << "Individual " << delimiter;
	for (int i = 0; i < N_; i++) {
		fout << "Chromosome " << i + 1 << delimiter;
	}

	fout << "Function value" << delimiter;

	fout << std::endl;
}

void GA::writeIndividualToCSV() {

	// writes iteration and best current swarm position
	fout << "Iteration: " << delimiter << iter << delimiter << "Best: " << delimiter;
	for (int i = 0; i < N_; i++) {
		fout << Population[0].chromosome[i] << delimiter; // best chromosome is at position 0
	}
	fout << Population[0].fitness << delimiter;

	fout << std::endl;

	// writes position and velocity
	for (int p = 0; p < population_size_; p++) {

		fout << p << delimiter;

		for (int i = 0; i < N_; i++) {
			fout << Population[p].chromosome[i] << delimiter;
		}

		fout << Population[p].fitness << delimiter;

		fout << "\n";
	}

	fout << std::endl;
}

void GA::closeCSV() {

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