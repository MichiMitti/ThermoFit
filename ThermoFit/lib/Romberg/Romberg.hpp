/*
	Numerical integration with Romberg's method
*/

#pragma once
#include "Romberg.hpp"
#include "../matrix/matrix.hpp"
#include <vector>
#include <iostream>
#include <ppl.h>
#include <cmath>
#include<functional>
#include "../core/global.hpp"

#ifndef Romberg_h
#define Romberg_h

using namespace concurrency;
using namespace std;

struct para_romberg {
	int max_steps = 20;
	int max_steps_ini = max_steps;
	int min_steps = 10;
	double acc = 1E-6;
	int debug = 0;
	int n_steps = 0;	// number of iterations in last run
};

// romberg function for comparing with class
double romberg_old(double(*func)(double), double *a, double *b);

// working function for ROM, similar to Powell
class func_ROM {
public:
	virtual void func_rom(double* x1, double* fval) { return; };
};

// Romberg class
class ROM {
public:
	ROM(para_romberg* paraROM, func_ROM *fun_ROM);	//ctor
	ROM(func_ROM *fun_ROM);							//ctor
	~ROM();											//dtor
	void set_paraROM(para_romberg* paraROM);
	double romberg(double *a, double *b);
	double romberg2(double *a, double *b);
	double romberg_old(double *a, double *b);
	void calc_fun(double* x1, double *fval);
	double calc_fun_double(double* x1);

private:
	int max_steps_;
	int max_steps_ini_;		// number of max_steps at creation of Romberg-object => array size is fixed to initial size, if max_steps decrease, other values should be set to 0
	int min_steps_;
	double acc_;			// accuracy
	double *h_;				// step size
	double a1;
	double fval_a, fval_b, fval_j;	// function values
	int debug_;
	para_romberg *paraROM_;
	func_ROM* fun_;
	matrix *r_;
};

#endif