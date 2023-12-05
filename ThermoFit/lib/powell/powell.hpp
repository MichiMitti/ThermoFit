/**
CS-11, powell.hpp
shared code for various calling functions
Purpose:
Solver class containing a Powell's Method based on the method described in
Powell, Michael J.D. „A Fortran Subroutine for Solving Systems of Nonlinear Algebraic Equations.“ In Numerical Methods for Nonlinear Algebraic Equations, Herausgeber: Philip Rabinowitz, 115-161. London: Gordon and Breach, Science Publishers Ltd., 1970.
and
Powell, Michael J.D. „A Hybrid Method for Nonlinear Equations.“ In Numerical Methods for Nonlinear Algebraic Equations, Herausgeber: Philip Rabinowitz, 87-114. London: Gordon and Breach, Science Publishers Ltd., 1970.

@authors:
Patrick Zimmermann, TU Graz
Patrick Krenn, TU Graz
@year 2018

revised 2023:
Gottfried Segner, TU Graz
*/

#ifndef POWELL_HPP

#define POWELL_HPP

#pragma once
//#include "matrix.hpp"
//#ifndef matrix_hpp // in case that its in an other folder
//#include "./matrix.hpp"
//#endif
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include<fstream>
#include<sstream>
#include<functional>
#include<vector>
#include "../core/global.hpp"

#include "../core/zero_function.hpp"
#include "../core/options.hpp"
#include "../core/optimizer.hpp"

using namespace std;

struct options_powell : options {
	std::string options_typ = "powell";	// stores options_powell name for checking if right options_powell are used in CallClass
	int dispOpt = 0;
	int maxIter = 1000;
	int maxFunEvals = 2000;
	double acc = 1e-6;
	double dstep = 1e-6;
	double dmax = 10;
};

//multi-dimensional

class powell {
public:
	powell();
	//powell(Zero_Function *fun, matrix* initGuess, options *opts);	//ctor		
	powell(Zero_Function* fun, std::vector<double>* initGuess, options_powell* opts);	//ctor	
	powell(const powell& rhs) {};						// cctor			   
	~powell();											// dtor	

	//friend void evalJacobian1(double *vec, double *temp_vec, matrix *Jacobian, int n, int m, Zero_Function *A, double *fval, double *temp_fval, double* dstep);
	friend void evalJacobian1(double* vec, double* temp_vec, double* Jacobian, int n, int m, Zero_Function* A, double* fval, double* temp_fval, double* dstep);
	friend class Zero_Function;

	//void get_Vals(Zero_Function *fun, matrix* initGuess, options *opts);
	void get_Vals(Zero_Function* fun, std::vector<double>* initGuess, options_powell* opts);
	//void calfun(double* x, double* f, double* Fsq);
	void calfun(std::vector<double>* x, std::vector<double>* f, double* Fsq);
	void solve();
	//void solve(matrix* initGuess);
	void solve(std::vector<double>* initGuess);

	//matrix solution();
	std::vector<double> solution();
	double solution_function();
	//void solution(matrix* solution);
	void solution(std::vector<double>* solution);
	int outflg();

private:
	Zero_Function* _fun;

	//matrix *Jac_, *invJac_, *W_, *invfield_;
	//matrix *Jacm_, *invJacm_, *Wm_, *invfieldm_;
	double* Jac_, * invJac_, * W_, * invfield_;
	double* Jacm_, * invJacm_, * Wm_, * invfieldm_;
	double* fval_, * lastfval_, * tempfval_, * x_, * lastx_, * tempx_, * v_, * jac_x_, * jac_fval_, * g_, * step_, * gtemp_, * phi_, * gamma_, * sig_, * initGuess_, * b_;
	double* initGuessm_, * xm_, * lastxm_, * tempxm_, * fvalm_, * lastfvalm_, * tempfvalm_, * vm_, * jac_xm_, * jac_fvalm_, * g1_, * stepm_, * gtempm_, * phim_, * gammam_, * sigm_, * bm_;

	int* w_;
	bool debug_, run_, flag_, Wflag_;
	int dispOpt_, N_, maxc_, maxfun_, maxiter_, maxi_, is_, ntest_, nt_, kk_, outflg_;
	int getVal_funcCtr_;
	double acc_, dmax_, dstep_, Fsq_, Fmin_, dd_, mu_, sp_, ds_, dss_, dn_, dm_, dmm_, dtest_, gm_, anmult_, dmult_, Phi_, tinc_, ss_, pj_, alpha_, sval_;
	double* a_, * solvec_;
	std::vector<double> x_vec, tempx_vec, fval_vec, tempfval_vec;

	void initialize();
	void reviseDelta();
	void reviseJac();
	void resetStep();
	void interchange();
	void setSpecialStep();
	void mainstep();
	void updateOmega();
	void setNewPoint();

	char separator;
	int nameWidth;
	int numWidth;
};

//void evalJacobian1(double *vec, double *temp_vec, matrix *Jacobian, int n, int m, Zero_Function *A, double *fval, double *temp_fval, double* dstep);
void evalJacobian1(double* vec, std::vector<double>* temp_vec, double* Jacobian, int n, int m, Zero_Function* A, double* fval, std::vector<double>* temp_fval, double* dstep);
void inverseMatrixArray(double* mat, double* invMat, double* invField, int dim);
void dispMatrix(double* mat, int col, int row);


//one-dimensional

class powell_1D {
public:
	powell_1D();
	powell_1D(Zero_Function* fun, double* initGuess, options_powell* opts);	//ctor								
	powell_1D(const powell& rhs) {};						// cctor			   

	friend void evalJacobian1(double*, double*, Zero_Function* fun, double*, double*);
	friend class Zero_Function;

	void get_Vals(Zero_Function* fun, double* initGuess, options_powell* opts);
	void calfun(double* x, double* f, double* Fsq);
	void solve();
	void solve(double* initGuess);
	double solution();
	int outflg();

private:
	Zero_Function* _fun;

	double fval_, lastfval_, tempfval_, x_, lastx_, tempx_, v_, g_, step_, Jac_, invJac_, gtemp_, W_, phi_, gamma_, sig_, initGuess_;
	int w_;
	bool debug_, run_, flag_, Wflag_;
	int dispOpt_, N_, maxc_, maxfun_, maxiter_, maxi_, is_, ntest_, nt_, kk_, outflg_;
	int getVal_funcCtr_;
	double acc_, dmax_, dstep_, Fsq_, Fmin_, dd_, mu_, sp_, ds_, dss_, dn_, dm_, dmm_, dtest_, gm_, anmult_, dmult_, Phi_, tinc_, ss_, pj_, alpha_, sval_;
	double a_, solvec_;

	void initialize();
	void reviseDelta();
	void reviseJac();
	void resetStep();
	void interchange();
	void setSpecialStep();
	void mainstep();
	void updateOmega();
	void setNewPoint();

	char separator;
	int nameWidth;
	int numWidth;
};

void evalJacobian1(double* vec, double* Jacobian, Zero_Function* A, double* fval, double* dstep);
#endif // !POWELL_HPP
