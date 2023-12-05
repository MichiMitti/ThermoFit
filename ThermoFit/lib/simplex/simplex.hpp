/* Downhill Simplex Optimizer
    _____ _                 __         
   / ___/(_)___ ___  ____  / /__  _  __
   \__ \/ / __ `__ \/ __ \/ / _ \| |/_/
  ___/ / / / / / / / /_/ / /  __/>  <  
 /____/_/_/ /_/ /_/ ____/_/\___/_/|_|  
                 /_/        

@author: Roland Nagl, TU Graz
*/

#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP

#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <limits>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ppl.h>
#include <vector>
#include <numeric>
#include "../core/global.hpp"

#include "../core/zero_function.hpp"
#include "../core/options.hpp"
#include "../core/optimizer.hpp"

//===== To Do =====\\ 
//* inherit from optimizer ( change Solve_DHS(...) to solve() )
//* 

struct options_dhs : options {
	std::string options_typ = "simplex";								// stores options_powell name for checking if right options_powell are used in CallClass
	int cnt_calfun = 0;		// counter for function calls
	int debug = 0;			// debug mode
	int logmode_ = 0;		// control verbosity of logging; 0 = no logging, 1 = log current best parameter set after each iteration
	int log_iter = 1;		// log parameters and function value every x iterations
	int disp_iter = 1;		// print out function value every x iterations
	double tol_ = 1E8 * std::numeric_limits<double>::epsilon();					//termination criteria
	std::vector<std::vector<double> > x_ = std::vector<std::vector<double> >();	//x: The Simplex
	int iterations_ = 10000;													// max iterations
	double div_ini_x_ = 20;			// divisor for calculation of initial simplex (init_para / divisor)
	int break_tol = 100;			//number of allowed successive violations of the termination criteria fun_val_old = fun_val_new
	double tol_fun_ = 1E-11;
	double a_ = 1.0, b_ = 1.0, g_ = 0.5, h_ = 0.5;	//coefficients
	//a: reflection  -> xr  
	//b: expansion   -> xe 
	//g: contraction -> xc
	//h: full contraction to x1
};

// Zero-function for Simplex
class func_DHS {

public:
	virtual void func_dhs(std::vector<double>(*vec), double* fval, options_dhs* para_sim) {};
};

class DHS {

	/*
	umzusetzende Änderungen:
	*erledigt* Implementierung als Klasse statt als namespace (wie bei Powell)
	*erledigt* statt Übergabe der Abbruchbedingungen (max. iterations und tol) setzen der Optionen über struct options_powell wie bei Powell
	*erledigt* setzen der Koeffizienten für Reflexion etc. über struct options_powell (Umbenennung auf alpha, beta, gamma, sigma)
	-- Optimierung mit Pointern
	-- Dynamische Speicherallozierung
	*/

public:
	DHS(Zero_Function* fun, options_dhs* para_sim);	//constructor
	DHS(Zero_Function* fun);					//constructor
	//~DHS() {};						//destructor

	// Solve method
	//std::vector<double> DHS::Simplex_Solve(std::vector<double> *init);
	std::vector<double> Simplex_Solve(std::vector<double>* init);

	// old versions
	//std::vector<double> Simplex(double f(std::vector<double>*), std::vector<double> *init);
	//std::vector<double> Simplex2(double f(std::vector<double>*), std::vector<double> *init);
	//std::vector<double> Simplex3(double (*f)(std::vector<double>*), std::vector<double> *init);

	// Call Zero-function
	void calfun(std::vector<double>* x, std::vector<double>* fvals, options_dhs* para_sim);

	int function_calls();

	double fval;						// sum of error function values

private:
	Zero_Function* _fun;
	options_dhs* para_sim_;
	int mark_branch_;					// Marker for simplex unit operation that leads to a updated parameter set, which need to be evaluated in the next iteration
	int logmode, debug_;
	int disp_iter_, log_iter_, cnt_;
	double tol;
	std::vector<double> fvals;			// vector of squared function values

	std::vector<std::vector<double>> x;
	int iterations;
	double div_ini_x;
	double tol_fun;		// tolerable difference of the function_value between 2 steps
	int break_tol;		//number of allowed successive violations of the termination criteria tol_fun ( fun_val_old = fun_val_new )
	double a, b, g, h;  //coefficients		   //a: reflection  -> xr  
	//b: expansion   -> xe 
	//g: contraction -> xc
	//h: full contraction to x1
};

#endif // !SIMPLEX_HPP
