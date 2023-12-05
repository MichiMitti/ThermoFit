// Downhill Simplex Optimizer
//    _____ _                 __         
//   / ___/(_)___ ___  ____  / /__  _  __
//   \__ \/ / __ `__ \/ __ \/ / _ \| |/_/
//  ___/ / / / / / / / /_/ / /  __/>  <  
// /____/_/_/ /_/ /_/ ____/_/\___/_/|_|  
//                 /_/                 

#include "simplex.hpp"

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
using namespace concurrency;

//constructor
DHS::DHS(Zero_Function* fun, options_dhs* para_sim) :
	fval(0.0),
	debug_((*para_sim).debug),
	_fun(fun),
	para_sim_(para_sim),
	logmode((*para_sim).logmode_),
	log_iter_((*para_sim).log_iter),
	disp_iter_((*para_sim).disp_iter),
	tol((*para_sim).tol_),
	x((*para_sim).x_),
	iterations((*para_sim).iterations_),
	div_ini_x((*para_sim).div_ini_x_),
	break_tol((*para_sim).break_tol),
	tol_fun((*para_sim).tol_fun_),
	a((*para_sim).a_),
	b((*para_sim).b_),
	g((*para_sim).g_),
	h((*para_sim).h_) {
	cnt_ = 0;
	(*para_sim).cnt_calfun = 0;
}

//constructor overload
DHS::DHS(Zero_Function *fun) {
	options_dhs para_sim;
	debug_ = para_sim.debug;
	_fun = fun;
	para_sim_ = &para_sim;
	logmode = para_sim.logmode_;
	log_iter_ = para_sim.log_iter;
	disp_iter_ = para_sim.disp_iter;
	tol = para_sim.tol_;
	x = para_sim.x_;
	iterations = para_sim.iterations_;
	div_ini_x = para_sim.div_ini_x_;
	break_tol = para_sim.break_tol;
	tol_fun = para_sim.tol_fun_;
	a = para_sim.a_;
	b = para_sim.b_;
	g = para_sim.g_;
	h = para_sim.h_;
	cnt_ = 0;
	para_sim.cnt_calfun = 0;
}

// Call Zero-function
void DHS::calfun(std::vector<double>* x, std::vector<double>* fvals, options_dhs *para_sim) {

	(*_fun).zero_function(x, fvals);
	(*para_sim).cnt_calfun++;
	return;
}

/*
	Simplex Algorithms
*/ 

std::vector<double> DHS::Simplex_Solve(std::vector<double> *init) {
	/*
		Simplex Algorithm
	*/

	// Marker for simplex unit operation that leads to a updated parameter set, which need to be evaluated in the next iteration
	// For reflection, reflection+expansion and contraction only the parameter set at x[xnp1] is updated and needs to be reevaluated in the next iteration
	// For full contraction, all parameter sets except x[x1] are updated and need to be reevaluated in the next iteration
	// 0 = initial value before first iteration cycle
	// 1 = reflection only
	// 2 = contraction only
	// 3 = full contraction
	// 4 = reflection + expansion
	mark_branch_ = 0;

	// Init logfile
	std::ofstream logfile;
	if (logmode) logfile.open("data_output/Simplex_log.txt");

	if (debug_ != 0) {
		std::cout << "\n===========================================\n"
			<< "Downhill Simplex" <<
			"\n===========================================\n\n" <<
			"maximum iterations: " << iterations << std::endl << "requested tolerance: " << tol_fun <<
			"\n\nrunning Simplex..." << std::endl;
	}


	int N = (*init).size();							//space dimension

	// initialize vectors and doubles that hold their respective error function values
	std::vector<double> xcentroid_old(N, 0);		// simplex center * (N+1)
	std::vector<double> xcentroid_new(N, 0);		// simplex center * (N+1)
	std::vector<double> fvals(N, 0.0);				// holds vector with current function values sumed up it is vf
	std::vector<double> vf(N + 1, 0);				// holds values of error function evaluated at all simplex verteces       
	std::vector<double> xg(N, 0);					// central vertex
	std::vector<double> xr(N, 0);					// reflection
	std::vector<double> xc(N, 0);					// contraction
	std::vector<double> xe(N, 0);					// expansion
	double diff = 0;
	double fxr = 0;
	double fxe = 0;
	double fxc = 0;
	double fun_val_new = 0;
	double delta_fun = 0;

	// Init vertex indices for current min. and max. values
	int x1 = 0, xn = 0, xnp1 = 0;         //x1:   f(x1) = min { f(x1), f(x2)...f(x_{n+1} }
										  //xnp1: f(xnp1) = max { f(x1), f(x2)...f(x_{n+1} }
										  //xn:   f(xn)<f(xnp1) && f(xn)> all other f(x_i)
	//iteration step number								  
	cnt_ = 0;
	(*para_sim_).cnt_calfun = 0;

	//initialize function value for termination criteria
	double fun_val_old = 0;

	//initialize counter for termination criteria
	int i_break = 0;

	//if no initial simplex is specified: construct the trial simplex based upon the initial guess parameters
	if (x.size() == 0) {
		if (debug_ != 0) { std::cout << "constructing initial simplex..."; }

		std::vector<double> del(*init);
		std::transform(del.begin(), del.end(), del.begin(),
			std::bind(std::divides<double>(), std::placeholders::_1, div_ini_x));	// divide all values in del by div_ini_x, e.g. 20 is picked if initial trial is assumed to be close to the result

		// create N verteces, each offsets another parameter by the corresponding value of the del-vector
		for (int i = 0; i < N; ++i) {
			std::vector<double> tmp(*init);
			tmp[i] += del[i];					// note: the (N+1)th vertex is the init-vector
			x.push_back(tmp);
		}
		x.push_back(*init);						//x.size()=N+1, x[i].size()=N

		// finished constructing the initial simplex
		if (debug_ != 0) { std::cout << "done" << std::endl; }
	}

	/*
		optimization cycle begins
	*/
	if (debug_ != 0 ) { std::cout << "starting Simplex optimization" << std::endl << std::endl; }
	for (cnt_ = 0; cnt_ < iterations; cnt_++) {
	
		// 0 = initial value before first iteration cycle
		// 1 = reflection only, only f(xnp1) has to be evaluated
		// 2 = contraction only
		// 3 = full contraction
		// 4 = reflection + expansion
		// 5 = reflection + expansion, but fxe > fxr, therefore set xnp1 = xe
		
		if (mark_branch_ == 0) {
		// 0 = evaluate funcction at every vertex at first iteration
			for (int i = 0; i < N + 1; ++i) {
				calfun(&x[i], &fvals, para_sim_);
				vf[i] = std::accumulate(fvals.begin(), fvals.end(), 0.0);
			}
		}

		else if (mark_branch_ == 3) {
		// 3 = full contraction, every vertex except x1 has changed and has to be evaluated
			for (int i = 0; i < N + 1; ++i) {
				if (i != x1)
				calfun(&x[i], &fvals, para_sim_);
				vf[i] = std::accumulate(fvals.begin(), fvals.end(), 0.0);
			}
		}
		
		//// Parallel version
		//parallel_for(0, N + 1, [&](int i) {
		//	vf[i] = f(&x[i]);
		//});

		//find index of max, second max, min of vf (xnp1, xn, x1)
		x1 = 0; xn = 0; xnp1 = 0;
		for (int i = 0; i < vf.size(); ++i) {
			if (vf[i] < vf[x1]) {
				x1 = i;
			}
			if (vf[i] > vf[xnp1]) {
				xnp1 = i;
			}
		}
		xn = x1;
		for (int i = 0; i < vf.size(); ++i) {
			if (vf[i]<vf[xnp1] && vf[i]>vf[xn])
				xn = i;
		}

		// Print iteration, current best function value, branch marker and total No. of function calls 
		if (debug_ != 0 && cnt_ % disp_iter_ == 0) std::cout << "iteration: " << cnt_ << ",\tfunction value: " << vf[x1] << "\tBranch marker: " << mark_branch_ << std::endl << "No. of function calls:\t" << (*para_sim_).cnt_calfun << std::endl << std::endl;

		// Write current best parameters to logfile
		if (logmode && (cnt_ % log_iter_ == 0)) {
			logfile << "------------------------------\niteration: " << cnt_ << ",\tfunction value: " << vf[x1] << "\tBranch marker: " << mark_branch_ << std::endl << "Function calls:\t" << (*para_sim_).cnt_calfun << std::endl << std::endl;
			for (std::vector<double>::iterator i = x[x1].begin(); i < x[x1].end(); i++) {
				logfile << *i << std::endl;
			}
			logfile << std::endl;
		}

		// Write all parameters to logfile
		if (logmode && (cnt_ % log_iter_ == 0)) {
			logfile << std::endl << "iteration: " << cnt_ << std::endl;
			for (int k = 0; k < x.size(); k++) {
				logfile << "function value: " << vf[k] << std::endl;
				for (std::vector<double>::iterator i = x[k].begin(); i < x[k].end(); i++) {
					logfile << *i << std::endl;
				}
				logfile << std::endl;
			}
		}

		// calculate sum of individual parameters at all vertices except the worst one for computing the mean value
		// i: cycles through verteces, j: cycles through parameters at the verteces i
		// xg: centroid of the N best verteces
		for (std::vector<double>::iterator i = xg.begin(); i < xg.end(); i++) {
			*i = 0;
		}
		for (int i = 0; i < x.size(); i++) {
			if (i != xnp1) {
				for (int j = 0; j < xg.size(); j++) {
					xg[j] += x[i][j];
				}
			}
		}

		// (xg + x[xnp1]) / N+1 = xcentroid_new
		std::transform(xg.begin(), xg.end(),
			x[xnp1].begin(), xcentroid_new.begin(), std::plus<double>());
		for (int i = 0; i < xcentroid_new.size(); i++) {
			xcentroid_new[i] /= (N + 1);
		}

		// compute mean by dividing the sum by the number of parameters
		for (int i = 0; i < xg.size(); i++) {
			xg[i] /= N;
		}

		//xg found, xcentroid_new updated
		diff = 0;

		/*
			termination criteria   
		*/   
		//check if the difference is less than the termination criteria
		for (int i = 0; i < N; ++i)
			diff += fabs(xcentroid_old[i] - xcentroid_new[i]);		//calculate the difference of the simplex centers

		// update best function value
		fun_val_new = vf[x1];

		//update simplex center
		xcentroid_old.swap(xcentroid_new);

		// compute difference of old and new best function value
		delta_fun = fabs(fun_val_new - fun_val_old);
		if (delta_fun < tol_fun || fun_val_new < tol) {
			if (fun_val_new < tol) {
				if (debug_ != 0) {
					std::cout << "\nrequested tolerance of the function value (" << tol << ") achieved!" << std::endl <<
						"fun_val_new = " << fun_val_new << std::endl;
					//std::cout << "fun_val_new = " << fun_val_new << std::endl << "fun_val_old = " << fun_val_old << std::endl;
					std::cout << "\noptimized function value: " << vf[x1] << "\niterations: " << cnt_ <<
						"\n\n===========================================\n" << std::endl;
				}
				fval = vf[x1];
				return x[x1];
			}
			i_break++;					// counts number of iteration cycles that didn't change the function value by at least tol_fun
			if (i_break > break_tol) {
				if (debug_ != 0) {
					std::cout << "\nrequested stability of the function value (" << tol_fun << ") achieved!" << std::endl <<
						"delta_fun_val = " << delta_fun << std::endl;
					//std::cout << "fun_val_new = " << fun_val_new << std::endl << "fun_val_old = " << fun_val_old << std::endl;
					std::cout << "\noptimized function value: " << vf[x1] << "\niterations: " << cnt_ <<
						"\n\n===========================================\n" << std::endl;
				}
				fval = vf[x1];
				return x[x1];
			}
		}
		else {
			// set counter to 0 if fun_val_new != fun_val_old
			i_break = 0;
			// update function value
			fun_val_old = fun_val_new;
		}

		// reflection xr, update xnp1 if xr is better
		for (std::vector<double>::iterator i = xr.begin(); i < xr.end(); i++) {
			*i = 0;
		}
		for (int i = 0; i < N; ++i)
			xr[i] = xg[i] + a*(xg[i] - x[xnp1][i]);

		//fxr = f(&xr);
		calfun(&xr, &fvals, para_sim_);
		fxr = std::accumulate(fvals.begin(), fvals.end(), 0.0);

		if (vf[x1] <= fxr && fxr <= vf[xn]) {
			std::copy(xr.begin(), xr.end(), x[xnp1].begin());
			vf[xnp1] = fxr;
			mark_branch_ = 1;
			if (debug_ != 0) std::cout << "reflection, xnp1 updated" << std::endl << std::endl;
		}
		//expansion xe, update xnp1 if xe is better:
		else if (fxr < vf[x1]) {
			for (std::vector<double>::iterator i = xe.begin(); i < xe.end(); i++) {
				*i = 0;
			}
			for (int i = 0; i < N; ++i)
				xe[i] = xr[i] + b*(xr[i] - xg[i]);
			//fxe = f(&xe);
			calfun(&xe, &fvals, para_sim_);
			fxe = std::accumulate(fvals.begin(), fvals.end(), 0.0);

			if (fxe < fxr) {
				std::copy(xe.begin(), xe.end(), x[xnp1].begin());
				vf[xnp1] = fxe;
				mark_branch_ = 4;
				if (debug_ != 0) std::cout << "reflection + expansion, fxr < vf[x1], xnp1 updated" << std::endl << std::endl;
			}
			else {
				std::copy(xr.begin(), xr.end(), x[xnp1].begin());
				vf[xnp1] = fxr;
				mark_branch_ = 5;
				if (debug_ != 0) std::cout << "reflection, fxr < vf[x1], xnp1 updated" << std::endl << std::endl;
			}
		}

		//contraction xc, update xnp1 if xc is better:
		else if (fxr > vf[xn]) {
			for (std::vector<double>::iterator i = xc.begin(); i < xc.end(); i++) {
				*i = 0;
			}
			for (int i = 0; i < N; ++i)
				xc[i] = xg[i] + g*(x[xnp1][i] - xg[i]);
			//fxc = f(&xc);
			calfun(&xc, &fvals, para_sim_);
			fxc = std::accumulate(fvals.begin(), fvals.end(), 0.0);
			if (fxc < vf[xnp1]) {
				std::copy(xc.begin(), xc.end(), x[xnp1].begin());
				vf[xnp1] = fxc;
				mark_branch_ = 2;
				if (debug_ != 0) {
					std::cout << "contraction, xnp1 updated" << std::endl << std::endl;
				}
			}
			else {	// full contraction in the direction of x1
				for (int i = 0; i < x.size(); i++) {
					if (i != x1) {
						for (int j = 0; j < N; j++)
							x[i][j] = x[x1][j] + h * (x[i][j] - x[x1][j]);
					}
				}
				mark_branch_ = 3;
				if (debug_ != 0) { std::cout << "full contraction, start new iteration cycle" << std::endl << std::endl; }
			}
		}	//contraction finished, xc is not used outside the scope
	}	//optimization is finished

	fval = vf[x1]; // sets error function value to be available outside of simplex

	if (cnt_ == iterations) {	//max number of iterations achieved before tol is satisfied
		if (debug_ != 0) { std::cout << "Iteration limit (" << cnt_ << ") achieved, result may not be optimal!" << std::endl; }
		if (logmode) logfile.close();
	}
	if (debug_ != 0) { std::cout << "\noptimized function value = " << vf[x1] << "\n========================================== = \n" << std::endl; }
	if (logmode) logfile.close();

	return x[x1];
}

int DHS::function_calls() {
	return (*para_sim_).cnt_calfun;
}

#pragma region old versions
//
//std::vector<double> DHS::Simplex(double f(std::vector<double>*),	//target function
//	std::vector<double> *init) {									//initial parameters
//
//	const int debug = 0;	//debug mode
//	std::cout << "\n===========================================\n" << "Downhill Simplex" <<
//		"\n===========================================\n\n" <<
//		"maximum iterations: " << iterations << std::endl << "requested tolerance: " << tol_fun <<
//		"\n\nrunning Simplex..." << std::endl;
//	int N = (*init).size();                         //space dimension
//
//	std::vector<double> xcentroid_old(N, 0);   //simplex center * (N+1)
//	std::vector<double> xcentroid_new(N, 0);   //simplex center * (N+1)
//	std::vector<double> vf(N + 1, 0);          //f evaluated at simplex verteces       
//	std::vector<double> xg(N, 0);
//	std::vector<double> xr(N, 0);
//	std::vector<double> xc(N, 0);
//	std::vector<double> xe(N, 0);
//	double diff = 0;
//	double fxr = 0;
//	double fun_val_new = 0;
//	double delta_fun = 0;
//
//	int x1 = 0, xn = 0, xnp1 = 0;         //x1:   f(x1) = min { f(x1), f(x2)...f(x_{n+1} }
//										  //xnp1: f(xnp1) = max { f(x1), f(x2)...f(x_{n+1} }
//										  //xn:   f(xn)<f(xnp1) && f(xn)> all other f(x_i)
//	cnt_ = 0;				//iteration step number
//	double fun_val_old = 0;	//initialize function value for termination criteria
//	int i_break = 0;			//initialize counter for termination criteria
//
//	if (x.size() == 0) {	//if no initial simplex is specified: construct the trial simplex 
//							//based upon the initial guess parameters
//		std::cout << "constructing initial simplex...";
//		std::vector<double> del(*init);
//		std::transform(del.begin(), del.end(), del.begin(),
//			std::bind2nd(std::divides<double>(), div_ini_x));	// divide all values in del by 20, 20 is picked 
//																//assuming initial trial close to true
//
//		for (int i = 0; i<N; ++i) {			// create N verteces, each offsets another parameter
//			std::vector<double> tmp(*init);		// by the corresponding value of the del-vector
//			tmp[i] += del[i];				// note: the (N+1)th vertex is the init-vector
//			x.push_back(tmp);
//		}
//		x.push_back(*init);	//x.size()=N+1, x[i].size()=N
//
//							// no inital xcentroid_old is needed: due to initializiation with {0} the difference
//							// fabs(xcentroid_old-xcentroid_new) is always > 0
//
//							////xcentroid = center of the simplex
//							//for (int i = 0; i < x.size(); i++) {				// i: cycles through verteces
//							//	for (int j = 0; j < xcentroid_old.size(); j++) {	// j: cycles through parameters at the verteces i
//							//		xcentroid_old[j] += (x[i][j])/(N+1);		// calculate sum of individual parameters at all
//							//	}												// vertices except the worst one for computing the mean value							
//							//}
//							//std::cout << "initial xcentroid_old = " << std::endl;
//							//V.display(&xcentroid_old);
//							//for (int i = 0; i < x.size(); i++) {
//							//	for (int j = 0; j < xcentroid_old.size(); j++) {
//							//		xcentroid_old[j] += x[i][j];
//							//		xcentroid_old[j] /= (N+1);
//							//	}
//							//}
//							//std::transform(init.begin(), init.end(), 
//							//	xcentroid_old.begin(), std::bind2nd(std::multiplies<D>(), N+1) );
//							//std::transform(xcentroid_old.begin(), xcentroid_old.end(),
//							//	xcentroid_old.begin(), std::bind2nd(std::divides<D>(), N+1));
//		std::cout << "done" << std::endl;
//	}	//constructing the simplex finished
//
//		//ofstream para_dat;		// purge datfile
//		//para_dat.open("para.dat", ios::trunc);
//		//if (debug == 1) { 
//		//	para_dat << "iteration\t" << "fun_value\t" << "best parameters" << std::endl;
//		//	para_dat.close();
//		//	para_dat.open("para.dat", ios::app);
//		//}
//
//		//optimization begins
//	std::cout << "starting optimization" << std::endl;
//	for (cnt_ = 0; cnt_ < iterations; cnt_++) {
//		if (cnt_ % disp_iter_ == 0) std::cout << "iteration: " << cnt_ << ",\tfunction value: " << f(&x[x1]) << std::endl;
//		//if (debug == 1) {	// write current parameters to logfile => slows program down!
//		//	para_dat << cnt_ << "\t" << fun_val_old << "\t{ ";
//		//	for (std::vector<double>::iterator i = x[x1].begin(); i < x[x1].end(); i++) {
//		//		para_dat << *i << "\t";
//		//	}
//		//	para_dat << "}" << std::endl;
//		//}
//
//		for (int i = 0; i < N + 1; ++i) {
//			vf[i] = f(&x[i]);
//		}
//		//parallel_for (0, N + 1, [&] (int i) {
//		//	vf[i] = f(&x[i]);
//		//});
//		x1 = 0; xn = 0; xnp1 = 0;	//find index of max, second max, min of vf.
//		for (int i = 0; i < vf.size(); ++i) {
//			if (vf[i] < vf[x1]) {
//				x1 = i;
//			}
//			if (vf[i] > vf[xnp1]) {
//				xnp1 = i;
//			}
//		}
//		xn = x1;
//
//		for (int i = 0; i < vf.size(); ++i) {
//			if (vf[i]<vf[xnp1] && vf[i]>vf[xn])
//				xn = i;
//		}
//		//x1, xn, xnp1 are found
//
//		for (std::vector<double>::iterator i = xg.begin(); i < xg.end(); i++) {
//			*i = 0;										//xg: centroid of the N best verteces
//		}
//
//		for (int i = 0; i < x.size(); i++) {			// i: cycles through verteces
//			if (i != xnp1) {
//				for (int j = 0; j < xg.size(); j++) {	// j: cycles through parameters at the verteces i
//					xg[j] += x[i][j];			// calculate sum of individual parameters at all  
//				}								// vertices except the worst one for computing the mean value
//			}
//		}
//		//for (int i = 0; i < x.size(); ++i) {
//		//	if(i!=xnp1)
//		//		std::transform(xg.begin(), xg.end(), x[i].begin(), xg.begin(), std::plus<D>() );
//		//}
//
//		//for (int i = 0; i < xg.size(); i++) {				// i: cycles through parameters
//		//	xcentroid_new[i] = (xg[i] + x[xnp1][i]) / (N + 1);	// store sum of all parameters (including the worst) in xcentroid_new
//		//}
//		std::transform(xg.begin(), xg.end(),
//			x[xnp1].begin(), xcentroid_new.begin(), std::plus<double>());
//
//		for (int i = 0; i < xg.size(); i++) {
//			xg[i] /= N;								// compute mean by dividing the sum by the number of parameters
//		}
//
//		/*			std::transform(xg.begin(), xg.end(), xg.begin(),
//		std::bind2nd(std::divides<D>(), N) );*/			// compute mean value by dividing sum by the number of parameters
//
//		//xg found, xcentroid_new updated
//		diff = 0;
//		//termination condition         
//		//see if the difference is less than the termination criteria
//		for (int i = 0; i < N; ++i)
//			diff += fabs(xcentroid_old[i] - xcentroid_new[i]);		//calculate the difference of the simplex centers
//
//		fun_val_new = f(&x[x1]);	// optimize
//		double fxnp1, fxn, fxg;
//		fxnp1 = f(&x[xnp1]);		// optimize
//		fxn = f(&x[xn]);			// optimize
//		fxg = f(&xg);				// optimize
//									/*	std::cout << "iteration: " << cnt_ << ", diff = " << diff << ", function value: " << fun_val_new
//									<< ", f(xnp1): " << fxnp1 << ", f(xn): " << fxn << ", f(xg): " << fxg << "\nx1:\n";*/
//
//									//if (diff < tol) {
//									//	std::cout << "\nrequested tolerance level (" << tol << ") achieved: diff = " << diff << std::endl
//									//		<< "iterations: " << cnt_ << std::endl;
//									//	return x[x1];;			//terminate the optimizer
//									//}
//									//else xcentroid_old.swap(xcentroid_new); //update simplex center
//		xcentroid_old.swap(xcentroid_new); //update simplex center
//
//		delta_fun = fabs(fun_val_new - fun_val_old);
//		//std::cout << "delta_fun = " << delta_fun << std::endl;
//		if (delta_fun < tol_fun || fun_val_new < tol) {
//			if (fun_val_new < tol) {
//				std::cout << "\nrequested tolerance of the function value (" << tol << ") achieved!" << std::endl <<
//					"fun_val_new = " << fun_val_new << std::endl;
//				//std::cout << "fun_val_new = " << fun_val_new << std::endl << "fun_val_old = " << fun_val_old << std::endl;
//				std::cout << "\noptimized function value: " << f(&x[x1]) << "\niterations: " << cnt_ <<
//					"\n\n===========================================\n" << std::endl;
//				return x[x1];;
//			}
//			i_break++;					// counts number of iteration cycles that didn't change the function value by at least tol_fun
//			if (i_break > break_tol) {
//				std::cout << "\nrequested stability of the function value (" << tol_fun << ") achieved!" << std::endl <<
//					"delta_fun_val = " << delta_fun << std::endl;
//				//std::cout << "fun_val_new = " << fun_val_new << std::endl << "fun_val_old = " << fun_val_old << std::endl;
//				std::cout << "\noptimized function value: " << f(&x[x1]) << "\niterations: " << cnt_ <<
//					"\n\n===========================================\n" << std::endl;
//				return x[x1];;
//			}
//		}
//		else {
//			i_break = 0;			// set counter to 0 if fun_val_new != fun_val_old
//			fun_val_old = fun_val_new; // update function value
//		}
//
//		//reflection:
//		for (std::vector<double>::iterator i = xr.begin(); i < xr.end(); i++) {
//			*i = 0;
//		}
//		for (int i = 0; i < N; ++i)
//			xr[i] = xg[i] + a*(xg[i] - x[xnp1][i]);
//		//reflection, xr found
//
//		fxr = f(&xr);	//record function at xr
//
//		if (vf[x1] <= fxr && fxr <= vf[xn])
//			//x[xnp1] = xr;
//			std::copy(xr.begin(), xr.end(), x[xnp1].begin());
//
//		//expansion:
//		else if (fxr < vf[x1]) {
//			for (std::vector<double>::iterator i = xe.begin(); i < xe.end(); i++) {
//				*i = 0;
//			}
//			for (int i = 0; i < N; ++i)
//				xe[i] = xr[i] + b*(xr[i] - xg[i]);
//			if (f(&xe) < fxr)
//				std::copy(xe.begin(), xe.end(), x[xnp1].begin());
//			else
//				std::copy(xr.begin(), xr.end(), x[xnp1].begin());
//		}	//expansion finished,  xe is not used outside the scope
//
//			//contraction:
//		else if (fxr > vf[xn]) {
//			for (std::vector<double>::iterator i = xc.begin(); i < xc.end(); i++) {
//				*i = 0;
//			}
//			for (int i = 0; i < N; ++i)
//				xc[i] = xg[i] + g*(x[xnp1][i] - xg[i]);
//			if (f(&xc) < vf[xnp1])
//				std::copy(xc.begin(), xc.end(), x[xnp1].begin());
//			else {
//				for (int i = 0; i < x.size(); ++i) {
//					if (i != x1) {
//						for (int j = 0; j < N; ++j)
//							x[i][j] = x[x1][j] + h * (x[i][j] - x[x1][j]);
//					}
//				}
//			}
//		}	//contraction finished, xc is not used outside the scope
//	}	//optimization is finished
//
//	if (cnt_ == iterations) {	//max number of iterations achieved before tol is satisfied
//		std::cout << "Iteration limit (" << cnt_ << ") achieved, result may not be optimal!" << std::endl;
//	}
//	std::cout << "\noptimized function value = " << f(&x[x1]) << "\n========================================== = \n" << std::endl;
//
//	//if (debug == 1) para_dat.close();
//	return x[x1];
//}
//
//std::vector<double> DHS::Simplex2(double f(std::vector<double>*),		//target function
//	std::vector<double> *init) {										//initial parameters
//
//	const int debug = 0;	//debug mode
//
//	std::cout << "\n===========================================\n"
//		<< "Downhill Simplex" <<
//		"\n===========================================\n\n" <<
//		"maximum iterations: " << iterations << std::endl << "requested tolerance: " << tol_fun <<
//		"\n\nrunning Simplex..." << std::endl;
//
//	int N = (*init).size();							//space dimension
//
//													// initialize vectors and doubles
//	std::vector<double> xcentroid_old(N, 0);		//simplex center * (N+1)
//	std::vector<double> xcentroid_new(N, 0);		//simplex center * (N+1)
//	std::vector<double> vf(N + 1, 0);				//f evaluated at simplex verteces       
//	std::vector<double> xg(N, 0);
//	std::vector<double> xr(N, 0);
//	std::vector<double> xc(N, 0);
//	std::vector<double> xe(N, 0);
//	double diff = 0;
//	double fxr = 0;
//	double fxe = 0;
//	double fxc = 0;
//	double fun_val_new = 0;
//	double delta_fun = 0;
//
//	int x1 = 0, xn = 0, xnp1 = 0;         //x1:   f(x1) = min { f(x1), f(x2)...f(x_{n+1} }
//										  //xnp1: f(xnp1) = max { f(x1), f(x2)...f(x_{n+1} }
//										  //xn:   f(xn)<f(xnp1) && f(xn)> all other f(x_i)
//										  //iteration step number								  
//	cnt_ = 0;
//
//	//initialize function value for termination criteria
//	double fun_val_old = 0;
//
//	//initialize counter for termination criteria
//	int i_break = 0;
//
//	//if no initial simplex is specified: construct the trial simplex based upon the initial guess parameters
//	if (x.size() == 0) {
//		std::cout << "constructing initial simplex...";
//		std::vector<double> del(*init);
//		std::transform(del.begin(), del.end(), del.begin(),
//			std::bind2nd(std::divides<double>(), div_ini_x));	// divide all values in del by div_ini_x, e.g. 20 is picked if initial trial is assumed to be close to the result
//
//																// create N verteces, each offsets another parameter by the corresponding value of the del-vector
//		for (int i = 0; i < N; ++i) {
//			std::vector<double> tmp(*init);
//			tmp[i] += del[i];					// note: the (N+1)th vertex is the init-vector
//			x.push_back(tmp);
//		}
//		x.push_back(*init);						//x.size()=N+1, x[i].size()=N
//
//												// no inital xcentroid_old is needed: due to initializiation with {0} the difference fabs(xcentroid_old-xcentroid_new) is always > 0
//
//												//constructing the simplex finished
//		std::cout << "done" << std::endl;
//	}
//
//	//optimization begins
//	std::cout << "starting optimization" << std::endl;
//	for (cnt_ = 0; cnt_ < iterations; cnt_++) {
//		if (cnt_ % disp_iter_ == 0) std::cout << "iteration: " << cnt_ << ",\tfunction value: " << f(&x[x1]) << std::endl;
//
//		//for (int i = 0; i < N + 1; ++i) {
//		//	vf[i] = f(&x[i]);
//		//}
//		parallel_for(0, N + 1, [&](int i) {
//			vf[i] = f(&x[i]);
//		});
//
//		//find index of max, second max, min of vf (xnp1, xn, x1)
//		x1 = 0; xn = 0; xnp1 = 0;
//		for (int i = 0; i < vf.size(); ++i) {
//			if (vf[i] < vf[x1]) {
//				x1 = i;
//			}
//			if (vf[i] > vf[xnp1]) {
//				xnp1 = i;
//			}
//		}
//		xn = x1;
//		for (int i = 0; i < vf.size(); ++i) {
//			if (vf[i]<vf[xnp1] && vf[i]>vf[xn])
//				xn = i;
//		}
//
//		//xg: centroid of the N best verteces
//		for (std::vector<double>::iterator i = xg.begin(); i < xg.end(); i++) {
//			*i = 0;
//		}
//
//		// calculate sum of individual parameters at all vertices except the worst one for computing the mean value
//		// i: cycles through verteces, j: cycles through parameters at the verteces i
//		for (int i = 0; i < x.size(); i++) {
//			if (i != xnp1) {
//				for (int j = 0; j < xg.size(); j++) {
//					xg[j] += x[i][j];
//				}
//			}
//		}
//
//		std::transform(xg.begin(), xg.end(),
//			x[xnp1].begin(), xcentroid_new.begin(), std::plus<double>());
//
//		// compute mean by dividing the sum by the number of parameters
//		for (int i = 0; i < xg.size(); i++) {
//			xg[i] /= N;
//		}
//
//		//xg found, xcentroid_new updated
//		diff = 0;
//
//		//termination condition         
//		//check if the difference is less than the termination criteria
//		for (int i = 0; i < N; ++i)
//			diff += fabs(xcentroid_old[i] - xcentroid_new[i]);		//calculate the difference of the simplex centers
//
//																	// update best function value
//		fun_val_new = vf[x1];
//
//		//update simplex center
//		xcentroid_old.swap(xcentroid_new);
//
//		delta_fun = fabs(fun_val_new - fun_val_old);
//		if (delta_fun < tol_fun || fun_val_new < tol) {
//			if (fun_val_new < tol) {
//				std::cout << "\nrequested tolerance of the function value (" << tol << ") achieved!" << std::endl <<
//					"fun_val_new = " << fun_val_new << std::endl;
//				//std::cout << "fun_val_new = " << fun_val_new << std::endl << "fun_val_old = " << fun_val_old << std::endl;
//				std::cout << "\noptimized function value: " << f(&x[x1]) << "\niterations: " << cnt_ <<
//					"\n\n===========================================\n" << std::endl;
//				return x[x1];;
//			}
//			i_break++;					// counts number of iteration cycles that didn't change the function value by at least tol_fun
//			if (i_break > break_tol) {
//				std::cout << "\nrequested stability of the function value (" << tol_fun << ") achieved!" << std::endl <<
//					"delta_fun_val = " << delta_fun << std::endl;
//				//std::cout << "fun_val_new = " << fun_val_new << std::endl << "fun_val_old = " << fun_val_old << std::endl;
//				std::cout << "\noptimized function value: " << f(&x[x1]) << "\niterations: " << cnt_ <<
//					"\n\n===========================================\n" << std::endl;
//				return x[x1];;
//			}
//		}
//		else {
//			// set counter to 0 if fun_val_new != fun_val_old
//			i_break = 0;
//			// update function value
//			fun_val_old = fun_val_new;
//		}
//
//		// reflection xr, update xnp1 if xr is better
//		for (std::vector<double>::iterator i = xr.begin(); i < xr.end(); i++) {
//			*i = 0;
//		}
//		for (int i = 0; i < N; ++i)
//			xr[i] = xg[i] + a*(xg[i] - x[xnp1][i]);
//
//		fxr = f(&xr);
//		if (vf[x1] <= fxr && fxr <= vf[xn])
//			std::copy(xr.begin(), xr.end(), x[xnp1].begin());
//
//		//expansion xe, update xnp1 if xe is better:
//		else if (fxr < vf[x1]) {
//			for (std::vector<double>::iterator i = xe.begin(); i < xe.end(); i++) {
//				*i = 0;
//			}
//			for (int i = 0; i < N; ++i)
//				xe[i] = xr[i] + b*(xr[i] - xg[i]);
//			fxe = f(&xe);
//			if (fxe < fxr)
//				std::copy(xe.begin(), xe.end(), x[xnp1].begin());
//			else
//				std::copy(xr.begin(), xr.end(), x[xnp1].begin());
//		}
//
//		//contraction xc, update xnp1 if xc is better:
//		else if (fxr > vf[xn]) {
//			for (std::vector<double>::iterator i = xc.begin(); i < xc.end(); i++) {
//				*i = 0;
//			}
//			for (int i = 0; i < N; ++i)
//				xc[i] = xg[i] + g*(x[xnp1][i] - xg[i]);
//			fxc = f(&xc);
//			if (fxc < vf[xnp1])
//				std::copy(xc.begin(), xc.end(), x[xnp1].begin());
//			else {
//				for (int i = 0; i < x.size(); i++) {
//					if (i != x1) {
//						for (int j = 0; j < N; j++)
//							x[i][j] = x[x1][j] + h * (x[i][j] - x[x1][j]);
//					}
//				}
//			}
//		}	//contraction finished, xc is not used outside the scope
//	}	//optimization is finished
//
//	if (cnt_ == iterations) {	//max number of iterations achieved before tol is satisfied
//		std::cout << "Iteration limit (" << cnt_ << ") achieved, result may not be optimal!" << std::endl;
//	}
//	std::cout << "\noptimized function value = " << f(&x[x1]) << "\n========================================== = \n" << std::endl;
//
//	return x[x1];
//}
//
//std::vector<double> DHS::Simplex3(double(*f)(std::vector<double>*),		//target function
//	std::vector<double> *init) {											//initial parameters
//
//	const int debug = 0;	//debug mode
//	std::ofstream logfile;
//	if (logmode) logfile.open("Simplex_log.txt");
//
//	std::cout << "\n===========================================\n"
//		<< "Downhill Simplex" <<
//		"\n===========================================\n\n" <<
//		"maximum iterations: " << iterations << std::endl << "requested tolerance: " << tol_fun <<
//		"\n\nrunning Simplex..." << std::endl;
//
//	int N = (*init).size();							//space dimension
//
//													// initialize vectors and doubles
//	std::vector<double> xcentroid_old(N, 0);		//simplex center * (N+1)
//	std::vector<double> xcentroid_new(N, 0);		//simplex center * (N+1)
//	std::vector<double> vf(N + 1, 0);				//f evaluated at simplex verteces       
//	std::vector<double> xg(N, 0);
//	std::vector<double> xr(N, 0);
//	std::vector<double> xc(N, 0);
//	std::vector<double> xe(N, 0);
//	double diff = 0;
//	double fxr = 0;
//	double fxe = 0;
//	double fxc = 0;
//	double fun_val_new = 0;
//	double delta_fun = 0;
//
//	int x1 = 0, xn = 0, xnp1 = 0;         //x1:   f(x1) = min { f(x1), f(x2)...f(x_{n+1} }
//										  //xnp1: f(xnp1) = max { f(x1), f(x2)...f(x_{n+1} }
//										  //xn:   f(xn)<f(xnp1) && f(xn)> all other f(x_i)
//										  //iteration step number								  
//	cnt_ = 0;
//
//	//initialize function value for termination criteria
//	double fun_val_old = 0;
//
//	//initialize counter for termination criteria
//	int i_break = 0;
//
//	//if no initial simplex is specified: construct the trial simplex based upon the initial guess parameters
//	if (x.size() == 0) {
//		std::cout << "constructing initial simplex...";
//		std::vector<double> del(*init);
//		std::transform(del.begin(), del.end(), del.begin(),
//			std::bind2nd(std::divides<double>(), div_ini_x));	// divide all values in del by div_ini_x, e.g. 20 is picked if initial trial is assumed to be close to the result
//
//																// create N verteces, each offsets another parameter by the corresponding value of the del-vector
//		for (int i = 0; i < N; ++i) {
//			std::vector<double> tmp(*init);
//			tmp[i] += del[i];					// note: the (N+1)th vertex is the init-vector
//			x.push_back(tmp);
//		}
//		x.push_back(*init);						//x.size()=N+1, x[i].size()=N
//
//												// no inital xcentroid_old is needed: due to initializiation with {0} the difference fabs(xcentroid_old-xcentroid_new) is always > 0
//
//												//constructing the simplex finished
//		std::cout << "done" << std::endl;
//	}
//
//	//optimization begins
//	std::cout << "starting optimization" << std::endl;
//	for (cnt_ = 0; cnt_ < iterations; cnt_++) {
//
//		for (int i = 0; i < N + 1; ++i) {
//			vf[i] = f(&x[i]);
//		}
//		//parallel_for(0, N + 1, [&](int i) {
//		//	vf[i] = f(&x[i]);
//		//});
//
//		//find index of max, second max, min of vf (xnp1, xn, x1)
//		x1 = 0; xn = 0; xnp1 = 0;
//		for (int i = 0; i < vf.size(); ++i) {
//			if (vf[i] < vf[x1]) {
//				x1 = i;
//			}
//			if (vf[i] > vf[xnp1]) {
//				xnp1 = i;
//			}
//		}
//		xn = x1;
//		for (int i = 0; i < vf.size(); ++i) {
//			if (vf[i]<vf[xnp1] && vf[i]>vf[xn])
//				xn = i;
//		}
//
//		// Print iteration and current best function value
//		if (cnt_ % disp_iter_ == 0) std::cout << "iteration: " << cnt_ << ",\tfunction value: " << vf[x1] << std::endl;
//
//		// Write current best parameters to logfile
//		if (logmode && (cnt_ % log_iter_ == 0)) {
//			logfile << "------------------------------\niteration: " << cnt_ << ",\tfunction value: " << vf[x1] << std::endl;
//			for (std::vector<double>::iterator i = x[x1].begin(); i < x[x1].end(); i++) {
//				logfile << *i << std::endl;
//			}
//			logfile << std::endl;
//		}
//
//		// Write all parameters to logfile
//		if (logmode && (cnt_ % log_iter_ == 0)) {
//			logfile << std::endl << "iteration: " << cnt_ << std::endl;
//			for (int k = 0; k < x.size(); k++) {
//				logfile << "function value: " << vf[k] << std::endl;
//				for (std::vector<double>::iterator i = x[k].begin(); i < x[k].end(); i++) {
//					logfile << *i << std::endl;
//				}
//				logfile << std::endl;
//			}
//		}
//
//		// calculate sum of individual parameters at all vertices except the worst one for computing the mean value
//		// i: cycles through verteces, j: cycles through parameters at the verteces i
//
//		//xg: centroid of the N best verteces
//		for (std::vector<double>::iterator i = xg.begin(); i < xg.end(); i++) {
//			*i = 0;
//		}
//		for (int i = 0; i < x.size(); i++) {
//			if (i != xnp1) {
//				for (int j = 0; j < xg.size(); j++) {
//					xg[j] += x[i][j];
//				}
//			}
//		}
//
//		// (xg + x[xnp1]) / N+1 = xcentroid_new
//		std::transform(xg.begin(), xg.end(),
//			x[xnp1].begin(), xcentroid_new.begin(), std::plus<double>());
//		for (int i = 0; i < xcentroid_new.size(); i++) {
//			xcentroid_new[i] /= (N + 1);
//		}
//
//		// compute mean by dividing the sum by the number of parameters
//		for (int i = 0; i < xg.size(); i++) {
//			xg[i] /= N;
//		}
//
//		//xg found, xcentroid_new updated
//		diff = 0;
//
//		//termination criterion         
//		//check if the difference is less than the termination criteria
//		for (int i = 0; i < N; ++i)
//			diff += fabs(xcentroid_old[i] - xcentroid_new[i]);		//calculate the difference of the simplex centers
//
//																	// update best function value
//		fun_val_new = vf[x1];
//
//		//update simplex center
//		xcentroid_old.swap(xcentroid_new);
//
//		delta_fun = fabs(fun_val_new - fun_val_old);
//		if (delta_fun < tol_fun || fun_val_new < tol) {
//			if (fun_val_new < tol) {
//				std::cout << "\nrequested tolerance of the function value (" << tol << ") achieved!" << std::endl <<
//					"fun_val_new = " << fun_val_new << std::endl;
//				//std::cout << "fun_val_new = " << fun_val_new << std::endl << "fun_val_old = " << fun_val_old << std::endl;
//				std::cout << "\noptimized function value: " << f(&x[x1]) << "\niterations: " << cnt_ <<
//					"\n\n===========================================\n" << std::endl;
//				return x[x1];;
//			}
//			i_break++;					// counts number of iteration cycles that didn't change the function value by at least tol_fun
//			if (i_break > break_tol) {
//				std::cout << "\nrequested stability of the function value (" << tol_fun << ") achieved!" << std::endl <<
//					"delta_fun_val = " << delta_fun << std::endl;
//				//std::cout << "fun_val_new = " << fun_val_new << std::endl << "fun_val_old = " << fun_val_old << std::endl;
//				std::cout << "\noptimized function value: " << f(&x[x1]) << "\niterations: " << cnt_ <<
//					"\n\n===========================================\n" << std::endl;
//				return x[x1];;
//			}
//		}
//		else {
//			// set counter to 0 if fun_val_new != fun_val_old
//			i_break = 0;
//			// update function value
//			fun_val_old = fun_val_new;
//		}
//
//		// reflection xr, update xnp1 if xr is better
//		for (std::vector<double>::iterator i = xr.begin(); i < xr.end(); i++) {
//			*i = 0;
//		}
//		for (int i = 0; i < N; ++i)
//			xr[i] = xg[i] + a*(xg[i] - x[xnp1][i]);
//
//		fxr = f(&xr);
//		if (vf[x1] <= fxr && fxr <= vf[xn]) {
//			std::copy(xr.begin(), xr.end(), x[xnp1].begin());
//			std::cout << "reflection successful, xnp1 updated" << std::endl;
//		}
//		//expansion xe, update xnp1 if xe is better:
//		else if (fxr < vf[x1]) {
//			for (std::vector<double>::iterator i = xe.begin(); i < xe.end(); i++) {
//				*i = 0;
//			}
//			for (int i = 0; i < N; ++i)
//				xe[i] = xr[i] + b*(xr[i] - xg[i]);
//			fxe = f(&xe);
//			if (fxe < fxr) {
//				std::copy(xe.begin(), xe.end(), x[xnp1].begin());
//				std::cout << "expansion successful, xnp1 updated" << std::endl;
//			}
//			else
//				std::copy(xr.begin(), xr.end(), x[xnp1].begin());
//		}
//
//		//contraction xc, update xnp1 if xc is better:
//		else if (fxr > vf[xn]) {
//			for (std::vector<double>::iterator i = xc.begin(); i < xc.end(); i++) {
//				*i = 0;
//			}
//			for (int i = 0; i < N; ++i)
//				xc[i] = xg[i] + g*(x[xnp1][i] - xg[i]);
//			fxc = f(&xc);
//			if (fxc < vf[xnp1]) {
//				std::copy(xc.begin(), xc.end(), x[xnp1].begin());
//				std::cout << "contraction successful, xnp1 updated" << std::endl;
//			}
//			else {	// full contraction in the direction of x1
//				for (int i = 0; i < x.size(); i++) {
//					if (i != x1) {
//						for (int j = 0; j < N; j++)
//							x[i][j] = x[x1][j] + h * (x[i][j] - x[x1][j]);
//					}
//				}
//				std::cout << "full contraction successful, start new iteration cycle" << std::endl;
//			}
//		}	//contraction finished, xc is not used outside the scope
//	}	//optimization is finished
//
//	if (cnt_ == iterations) {	//max number of iterations achieved before tol is satisfied
//		std::cout << "Iteration limit (" << cnt_ << ") achieved, result may not be optimal!" << std::endl;
//		if (logmode) logfile.close();
//	}
//	std::cout << "\noptimized function value = " << f(&x[x1]) << "\n========================================== = \n" << std::endl;
//	if (logmode) logfile.close();
//
//	return x[x1];
//}
#pragma endregion