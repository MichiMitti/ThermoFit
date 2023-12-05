#pragma once
#include "../../lib/structs.hpp"

struct parameters {

	double a_soll = 4.0;
	double b_soll = 4.0;


	//variables for old test cases
	double a_powell;
	double b_powell;

	double x_0;
	double x_1;


};

template <typename T> double sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

struct parameters_old {
	// Includes all parameters settings and variables that are necessary in more than one function 
		//bool test = true;
};

void init(parameters(*params));

//template <typename T> double sgn(T val) {
//	return (T(0) < val) - (val < T(0));
//}
