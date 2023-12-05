#pragma once
#include <chrono>
#include <iostream>

#include "init_call.hpp"
#include "../lib/pso/networks_psorl.hpp"


int main() {

	//**************** Test Call Function ****************\\ 

	auto time_start = std::chrono::high_resolution_clock::now();

	// Show cases for Call Function
	init_call_function(0);

	auto time_end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds> (time_end - time_start);
	std::cout << std::endl << "Finished";
	std::cout << std::endl << "The run was finished in " << duration.count()/1000.0 << " milliseconds" << std::endl << std::endl;

	system("pause");

	return 0;
}