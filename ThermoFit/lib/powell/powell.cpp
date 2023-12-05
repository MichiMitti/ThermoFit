#pragma once
#include "powell.hpp"

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()


//multi-dimensional powell (also works for 1D-problems but 1D-powell is faster)

	// ctor
powell::powell() {
	getVal_funcCtr_ = 0;
};

//powell::powell(Zero_Function* fun, matrix* initGuess, options* opts) :
//	N_((*initGuess).rows()),
//	_fun(fun),
//	dispOpt_((*opts).dispOpt),
//	acc_((*opts).acc),
//	maxiter_((*opts).maxIter),
//	dstep_ ((*opts).dstep),
//	dmax_((*opts).dmax),
//	maxfun_((*opts).maxFunEvals) {
//	debug_ = 0;
//	solvec_ = new double[N_];
//	for (int i = 0; i < N_; i++)
//		solvec_[i] = (*initGuess)(i, 0);
//	getVal_funcCtr_ = 1;
//	initGuessm_ = new double [N_];
//	initGuess_ = initGuessm_;
//	xm_ = new double [N_];
//	x_ = xm_;
//	lastxm_ = new double [N_];
//	lastx_ = lastxm_;
//	tempxm_ = new double [N_];
//	tempx_ = tempxm_;
//	fvalm_ = new double [N_];
//	fval_ = fvalm_;
//	lastfvalm_ = new double [N_];
//	lastfval_ = lastfvalm_;
//	tempfvalm_ = new double [N_];
//	tempfval_ = tempfvalm_;
//	vm_ = new double [N_];
//	v_ = vm_;
//	jac_xm_ = new double[N_];
//	jac_x_ = jac_xm_;	
//	jac_fvalm_ = new double[N_];
//	jac_fval_ = jac_fvalm_;
//	g1_ = new double [N_];
//	g_ = g1_;
//	stepm_ = new double [N_];
//	step_ = stepm_;
//	//Jacm_ = new matrix((size_t)N_, (size_t)N_);
//	Jacm_ = new double[N_ * N_];
//	Jac_ = Jacm_;
//	//invJacm_ = new matrix((size_t)N_, (size_t)N_);
//	invJacm_ = new double[N_ * N_];
//	invJac_ = invJacm_;
//	gtempm_ = new double [N_];
//	gtemp_ = gtempm_;
//	//Wm_ = new matrix((size_t)N_, (size_t)N_);
//	Wm_ = new double[N_ * N_];
//	W_ = Wm_;
//	//invfieldm_ = new matrix((size_t)N_, (size_t)2*N_);
//	invfieldm_ = new double[2 * N_ * N_];
//	invfield_ = invfieldm_;
//	phim_ = new double [N_];
//	(phi_) = phim_;
//	gammam_ = new double [N_];
//	gamma_ = gammam_;
//	sigm_ = new double [N_];
//	sig_ = sigm_;
//	bm_ = new double[N_*N_];
//	b_ = bm_;
//	w_ = new int[N_];
//	a_ = new double[N_];
//};
powell::powell(Zero_Function* fun, std::vector<double>* initGuess, options_powell* opts) :
	N_((*initGuess).size()),
	_fun(fun),
	dispOpt_((*opts).dispOpt),
	acc_((*opts).acc),
	maxiter_((*opts).maxIter),
	dstep_((*opts).dstep),
	dmax_((*opts).dmax),
	maxfun_((*opts).maxFunEvals) {
	debug_ = 0;
	solvec_ = new double[N_];
	for (int i = 0; i < N_; i++)
		solvec_[i] = (*initGuess)[i];
	getVal_funcCtr_ = 1;
	initGuessm_ = new double[N_];
	initGuess_ = initGuessm_;
	xm_ = new double[N_];
	x_ = xm_;
	lastxm_ = new double[N_];
	lastx_ = lastxm_;
	tempxm_ = new double[N_];
	tempx_ = tempxm_;
	fvalm_ = new double[N_];
	fval_ = fvalm_;
	lastfvalm_ = new double[N_];
	lastfval_ = lastfvalm_;
	tempfvalm_ = new double[N_];
	tempfval_ = tempfvalm_;
	vm_ = new double[N_];
	v_ = vm_;
	jac_xm_ = new double[N_];
	jac_x_ = jac_xm_;
	jac_fvalm_ = new double[N_];
	jac_fval_ = jac_fvalm_;
	g1_ = new double[N_];
	g_ = g1_;
	stepm_ = new double[N_];
	step_ = stepm_;
	//Jacm_ = new matrix((size_t)N_, (size_t)N_);
	Jacm_ = new double[N_ * N_];
	Jac_ = Jacm_;
	//invJacm_ = new matrix((size_t)N_, (size_t)N_);
	invJacm_ = new double[N_ * N_];
	invJac_ = invJacm_;
	gtempm_ = new double[N_];
	gtemp_ = gtempm_;
	//Wm_ = new matrix((size_t)N_, (size_t)N_);
	Wm_ = new double[N_ * N_];
	W_ = Wm_;
	//invfieldm_ = new matrix((size_t)N_, (size_t)2*N_);
	invfieldm_ = new double[2 * N_ * N_];
	invfield_ = invfieldm_;
	phim_ = new double[N_];
	(phi_) = phim_;
	gammam_ = new double[N_];
	gamma_ = gammam_;
	sigm_ = new double[N_];
	sig_ = sigm_;
	bm_ = new double[N_ * N_];
	b_ = bm_;
	w_ = new int[N_];
	a_ = new double[N_];
	x_vec.resize(N_);
	fval_vec.resize(N_);



};

// dtor
powell::~powell() {
	if (getVal_funcCtr_ != 0) {
		//double arrays
		delete[] a_;
		delete[] solvec_;
		delete[] w_;
		delete[] initGuessm_;
		delete[] xm_;
		delete[] lastxm_;
		delete[] tempxm_;
		delete[] fvalm_;
		delete[] lastfvalm_;
		delete[] tempfvalm_;
		delete[] vm_;
		delete[] jac_xm_;
		delete[] jac_fvalm_;
		delete[] g1_;
		delete[] stepm_;
		delete[] gtempm_;
		delete[] phim_;
		delete[] gammam_;
		delete[] sigm_;
		delete[] bm_;


		//matrices
		delete[] Jacm_;
		delete[] invJacm_;
		delete[] Wm_;
		delete[] invfieldm_;

	}
};

//void powell::get_Vals(Zero_Function *fun, matrix* initGuess, options *opts) {
//	N_ = (*initGuess).rows();
//	_fun = fun;
//	if (getVal_funcCtr_ == 0) {
//		w_ = new int[N_];
//		a_ = new double[N_];
//		solvec_ = new double[N_];
//		for (int i = 0; i < N_; i++)
//			solvec_[i] = (*initGuess)(i, 0);
//		getVal_funcCtr_ = 1;
//		initGuessm_ = new double[N_];
//		initGuess_ = initGuessm_;
//		xm_ = new double [N_];
//		x_ = xm_;
//		lastxm_ = new double[N_];
//		lastx_ = lastxm_;
//		tempxm_ = new double[N_];
//		tempx_ = tempxm_;
//		fvalm_ = new double[N_];
//		fval_ = fvalm_;
//		lastfvalm_ = new double[N_];
//		lastfval_ = lastfvalm_;
//		tempfvalm_ = new double[N_];
//		tempfval_ = tempfvalm_;
//		vm_ = new double[N_];
//		v_ = vm_;
//		jac_xm_ = new double[N_];
//		jac_x_ = jac_xm_;
//		jac_fvalm_ = new double[N_];
//		jac_fval_ = jac_fvalm_;
//		g1_ = new double[N_];
//		g_ = g1_;
//		stepm_ = new double[N_];
//		step_ = stepm_;
//		//Jacm_ = new matrix((size_t)N_, (size_t)N_);
//		Jacm_ = new double[N_ * N_];
//		Jac_ = Jacm_;
//		//invJacm_ = new matrix((size_t)N_, (size_t)N_);
//		invJacm_ = new double[N_ * N_];
//		invJac_ = invJacm_;
//		gtempm_ = new double[N_];
//		gtemp_ = gtempm_;
//		//Wm_ = new matrix((size_t)N_, (size_t)N_);
//		Wm_ = new double[N_ * N_];
//		W_ = Wm_;
//		//invfieldm_ = new matrix((size_t)N_, (size_t)2 * N_);
//		invfieldm_ = new double[2 * N_ * N_];
//		invfield_ = invfieldm_;
//		phim_ = new double[N_];
//		phi_ = phim_;
//		gammam_ = new double[N_];
//		gamma_ = gammam_;
//		sigm_ = new double[N_];
//		sig_ = sigm_;
//		bm_ = new double[N_*N_];
//		b_ = bm_;
//	}
//	for (int i = 0; i < N_; i++)
//		solvec_[i] = (*initGuess)(i, 0);
//	dispOpt_ = (*opts).dispOpt;
//	acc_ = (*opts).acc;
//	dstep_ = (*opts).dstep;
//	dmax_ = (*opts).dmax;
//	maxiter_ = (*opts).maxIter;
//	maxfun_ = (*opts).maxFunEvals;
//	debug_ = 0;
//	getVal_funcCtr_++;
//};
void powell::get_Vals(Zero_Function* fun, std::vector<double>* initGuess, options_powell* opts) {
	N_ = (*initGuess).size();
	_fun = fun;
	if (getVal_funcCtr_ == 0) {
		w_ = new int[N_];
		a_ = new double[N_];
		solvec_ = new double[N_];
		for (int i = 0; i < N_; i++)
			solvec_[i] = (*initGuess)[i];
		getVal_funcCtr_ = 1;
		initGuessm_ = new double[N_];
		initGuess_ = initGuessm_;
		xm_ = new double[N_];
		x_ = xm_;
		lastxm_ = new double[N_];
		lastx_ = lastxm_;
		tempxm_ = new double[N_];
		tempx_ = tempxm_;
		tempx_vec.resize(N_);
		fvalm_ = new double[N_];
		fval_ = fvalm_;
		tempfval_vec.resize(N_);
		lastfvalm_ = new double[N_];
		lastfval_ = lastfvalm_;
		tempfvalm_ = new double[N_];
		tempfval_ = tempfvalm_;
		vm_ = new double[N_];
		v_ = vm_;
		jac_xm_ = new double[N_];
		jac_x_ = jac_xm_;
		jac_fvalm_ = new double[N_];
		jac_fval_ = jac_fvalm_;
		g1_ = new double[N_];
		g_ = g1_;
		stepm_ = new double[N_];
		step_ = stepm_;
		//Jacm_ = new matrix((size_t)N_, (size_t)N_);
		Jacm_ = new double[N_ * N_];
		Jac_ = Jacm_;
		//invJacm_ = new matrix((size_t)N_, (size_t)N_);
		invJacm_ = new double[N_ * N_];
		invJac_ = invJacm_;
		gtempm_ = new double[N_];
		gtemp_ = gtempm_;
		//Wm_ = new matrix((size_t)N_, (size_t)N_);
		Wm_ = new double[N_ * N_];
		W_ = Wm_;
		//invfieldm_ = new matrix((size_t)N_, (size_t)2 * N_);
		invfieldm_ = new double[2 * N_ * N_];
		invfield_ = invfieldm_;
		phim_ = new double[N_];
		phi_ = phim_;
		gammam_ = new double[N_];
		gamma_ = gammam_;
		sigm_ = new double[N_];
		sig_ = sigm_;
		bm_ = new double[N_ * N_];
		b_ = bm_;
		x_vec.resize(N_);
		fval_vec.resize(N_);
	}
	for (int i = 0; i < N_; i++)
		solvec_[i] = (*initGuess)[i];
	dispOpt_ = (*opts).dispOpt;
	acc_ = (*opts).acc;
	dstep_ = (*opts).dstep;
	dmax_ = (*opts).dmax;
	maxiter_ = (*opts).maxIter;
	maxfun_ = (*opts).maxFunEvals;
	debug_ = 0;
	getVal_funcCtr_++;

};

void powell::solve() {
	// initialization of matrices/arrays
	//matrix initGuessm_(N_, 1);
	//initGuess_ = &initGuessm_;
	//matrix xm_(N_, 1);

	for (int i = 0; i < N_; i++)
		initGuessm_[i] = solvec_[i];
	x_ = xm_;
	for (int i = 0; i < N_; i++)
		x_[i] = initGuess_[i];
	//matrix lastxm_(N_, 1);
	//lastx_ = &lastxm_; 
	//matrix tempxm_(N_, 1); 
	//tempx_ = &tempxm_;
	/*matrix fvalm_(N_, 1);
	fval_ = &fvalm_;
	matrix lastfvalm_(N_, 1);
	lastfval_ = &lastfvalm_;
	matrix tempfvalm_(N_, 1);
	tempfval_ = &tempfvalm_;
	matrix vm_(N_, 1);
	v_ = &vm_;
	matrix gm_(N_, 1);
	g_ = &gm_;
	matrix stepm_(N_, 1);
	step_ = &stepm_;
	matrix Jacm_(N_, N_);
	Jac_ = &Jacm_;
	matrix invJacm_(N_, N_);
	invJac_ = &invJacm_;
	matrix gtempm_(N_, 1);
	gtemp_ = &gtempm_;
	matrix Wm_(N_, N_);
	W_ = &Wm_;
	matrix phim_(N_, 1);
	(phi_) = &phim_;
	matrix steptransposedm_(1, N_);
	steptransposed_ = &steptransposedm_;
	matrix gammam_(N_, 1);
	gamma_ = &gammam_;
	matrix sigm_(N_, 1);
	sig_ = &sigm_;*/

	initialize();
	if (dispOpt_ != 0)
		std::cout << "Start solving... Requested accuracy = " << acc_ << std::endl;
	if (dispOpt_ == 1) {
		std::cout << right << setw(numWidth) << setfill(separator) << "Iteration";
		std::cout << right << setw(numWidth) << setfill(separator) << "FunEvals";
		for (int i = 0; i < N_; i++)
			//std::cout << right << setw(numWidth) << setfill(separator) << SSTR("x" << i + 1);
			for (int i = 0; i < N_; i++)
				//std::cout << right << setw(numWidth) << setfill(separator) << SSTR("f" << i + 1);
				std::cout << right << setw(numWidth) << setfill(separator) << "F(x)";
		std::cout << right << setw(numWidth) << setfill(separator) << "Delta" << std::endl;
	}
	while (run_) {

		x_vec.assign(x_, x_ + N_);
		fval_vec.assign(fval_, fval_ + N_);

		calfun(&x_vec, &fval_vec, &Fsq_);

		x_ = x_vec.data();
		fval_ = fval_vec.data();

		maxc_++;
		maxi_++;
		if (Fsq_ > acc_) {
			if (is_ == 1 || is_ == 4) {
				if (Fsq_ >= Fmin_) {
					if (dd_ - dss_ <= 0) {
						ntest_--;
						if (ntest_ < 0) {
							if (dispOpt_ != 0)
								std::cout << "Powell terminated: failed to decrease Fsq_ using a new Jacobian." << std::endl;
							run_ = 0;
							outflg_ = 4;
							resetStep();
						}
						else if (ntest_ == 0) {
							if (dispOpt_ != 0)
								std::cout << "Powell terminated: " << nt_ << " calls of calfun ineffective." << std::endl;
							run_ = 0;
							outflg_ = 3;
							resetStep();
						}
					}
				}
				else {
					ntest_ = nt_;
				}
			}
			if (run_) {
				if (maxc_ <= maxfun_ && maxi_ <= maxiter_) {
					if (dispOpt_ == 1) {
						std::cout << right << setw(numWidth) << setfill(separator) << maxi_;
						std::cout << right << setw(numWidth) << setfill(separator) << maxc_;
						for (int i = 0; i < N_; i++)
							std::cout << right << setw(numWidth) << setfill(separator) << x_[i];
						for (int i = 0; i < N_; i++)
							std::cout << right << setw(numWidth) << setfill(separator) << fval_[i];
						std::cout << right << setw(numWidth) << setfill(separator) << Fsq_;
						std::cout << right << setw(numWidth) << setfill(separator) << sqrt(dd_) << std::endl;
					}
					//std::cout << maxi_ << "\t\t\t" << maxc_ << "\t" << (*x_)(0, 0) << "\t" << (*x_)(1, 0) << "\t" << (*fval_)(0, 0) << "\t" << (*fval_)(1, 0) << "\t" << Fsq_ << "\t" << sqrt(dd_) << std::endl;
					switch (is_) {
					case 1:

						reviseDelta();
						break;
					case 2:

						reviseJac();
						break;
					case 3:

						tempx_vec.assign(jac_x_, jac_x_ + N_);
						tempfval_vec.assign(jac_fval_, jac_fval_ + N_);

						evalJacobian1(x_, &tempx_vec, Jac_, N_, N_, _fun, fval_, &tempfval_vec, &dstep_);

						inverseMatrixArray(Jac_, invJac_, invfield_, N_);
						//(*Jac_).mat_inverse(invJac_, invfield_);
						if (debug_) {
							std::cout << "J = " << std::endl;
							dispMatrix(Jac_, N_, N_);
							std::cout << "J-1 = " << std::endl;
							dispMatrix(invJac_, N_, N_);
						}
						maxc_ += N_;
						mainstep();
						break;
					case 4:
						if (Fsq_ > Fmin_) {
							interchange();
						}
						else {

							is_ = 2;
							setSpecialStep();
						}
						break;
					case 5:

						Fmin_ = Fsq_;
						for (int i = 0; i < N_; i++) {
							lastx_[i] = x_[i];
							lastfval_[i] = fval_[i];
						}
						is_ = 3;
						break;
					}
				}
				else {
					if (maxc_ >= maxfun_) {
						if (dispOpt_ != 0)
							std::cout << "Powell terminated: maximum number of function evaluations reached." << std::endl;
						run_ = 0;
						outflg_ = 2;
					}
					else {
						if (dispOpt_ != 0)
							std::cout << "Powell terminated: maximum number of Iterations reached." << std::endl;
						run_ = 0;
						outflg_ = 1;
					}
					if (Fsq_ > Fmin_) {
						resetStep();
					}
				}
			}
		}
		else {
			if (dispOpt_ == 1) {
				std::cout << right << setw(numWidth) << setfill(separator) << maxi_;
				std::cout << right << setw(numWidth) << setfill(separator) << maxc_;
				for (int i = 0; i < N_; i++)
					std::cout << right << setw(numWidth) << setfill(separator) << x_[i];
				for (int i = 0; i < N_; i++)
					std::cout << right << setw(numWidth) << setfill(separator) << fval_[i];
				std::cout << right << setw(numWidth) << setfill(separator) << Fsq_;
				std::cout << right << setw(numWidth) << setfill(separator) << sqrt(dd_) << std::endl;
			}
			for (int i = 0; i < N_; i++)
				solvec_[i] = x_[i];
			if (isnan(Fsq_) == 1) {
				if (dispOpt_ != 0)
					std::cout << "Powell terminated: Function value is a complex number." << std::endl;
				run_ = 0;
				outflg_ = 6;
			}
			else {
				if (dispOpt_ != 0)
					std::cout << "Powell terminated: Solution found." << std::endl;
				run_ = 0;
				outflg_ = 0;
			}

		}
		if (debug_)
			std::cin.get();
	}
};

//void powell::solve(matrix * initGuess) {
//	// initialization of matrices/arrays
//	//matrix initGuessm_(N_, 1);
//	//initGuess_ = &initGuessm_;
//	//matrix xm_(N_, 1);
//
//	for (int i = 0; i < N_; i++)
//	initGuessm_[i] = solvec_[i] = (*initGuess)(i, 0);
//	
//	x_ = xm_;
//
//	for (int i = 0; i < N_; i++)
//	x_[i] = initGuess_[i];
//
//	/*matrix lastxm_(N_, 1);
//	lastx_ = &lastxm_; 
//	matrix tempxm_(N_, 1); 
//	tempx_ = &tempxm_;
//	matrix fvalm_(N_, 1);
//	fval_ = &fvalm_;
//	matrix lastfvalm_(N_, 1);
//	lastfval_ = &lastfvalm_;
//	matrix tempfvalm_(N_, 1);
//	tempfval_ = &tempfvalm_;
//	matrix vm_(N_, 1);
//	v_ = &vm_;
//	matrix gm_(N_, 1);
//	g_ = &gm_;
//	matrix stepm_(N_, 1);
//	step_ = &stepm_;
//	matrix Jacm_(N_, N_);
//	Jac_ = &Jacm_;
//	matrix invJacm_(N_, N_);
//	invJac_ = &invJacm_;
//	matrix gtempm_(N_, 1);
//	gtemp_ = &gtempm_;
//	matrix Wm_(N_, N_);
//	W_ = &Wm_;
//	matrix phim_(N_, 1);
//	(phi_) = &phim_;
//	matrix steptransposedm_(1, N_);
//	steptransposed_ = &steptransposedm_;
//	matrix gammam_(N_, 1);
//	gamma_ = &gammam_;
//	matrix sigm_(N_, 1);
//	sig_ = &sigm_;*/
//
//	initialize();
//	if (dispOpt_ != 0)
//		std::cout << "Start solving... Requested accuracy = " << acc_ << std::endl;
//	if (dispOpt_ == 1) {
//		std::cout << right << setw(numWidth) << setfill(separator) << "Iteration";
//		std::cout << right << setw(numWidth) << setfill(separator) << "FunEvals";
//		for (int i = 0; i < N_; i++)
//			//std::cout << right << setw(numWidth) << setfill(separator) << SSTR("x" << i + 1);
//		for (int i = 0; i < N_; i++)
//			//std::cout << right << setw(numWidth) << setfill(separator) << SSTR("f" << i + 1);
//		std::cout << right << setw(numWidth) << setfill(separator) << "F(x)";
//		std::cout << right << setw(numWidth) << setfill(separator) << "Delta" << std::endl;
//	}
//	while (run_) {
//		calfun(x_, fval_, &Fsq_);
//		maxc_++;
//		maxi_++;
//		if (Fsq_ > acc_) {
//			if (is_ == 1 || is_ == 4) {
//				if (Fsq_ >= Fmin_) {
//					if (dd_ - dss_ <= 0) {
//						ntest_--;
//						if (ntest_ < 0) {
//							if (dispOpt_ != 0)
//								std::cout << "Powell terminated: failed to decrease Fsq_ using a new Jacobian." << std::endl;
//							run_ = 0;
//							outflg_ = 4;
//							resetStep();
//						}
//						else if (ntest_ == 0) {
//							if (dispOpt_ != 0)
//								std::cout << "Powell terminated: " << nt_ << " calls of calfun ineffective." << std::endl;
//							run_ = 0;
//							outflg_ = 3;
//							resetStep();
//						}
//					}
//				}
//				else {
//					ntest_ = nt_;
//				}
//			}
//			if (run_) {
//				if (maxc_ <= maxfun_ && maxi_ <= maxiter_) {
//					if (dispOpt_ == 1) {
//						std::cout << right << setw(numWidth) << setfill(separator) << maxi_;
//						std::cout << right << setw(numWidth) << setfill(separator) << maxc_;
//						for (int i = 0; i < N_; i++)
//							std::cout << right << setw(numWidth) << setfill(separator) << x_[i];
//						for (int i = 0; i < N_; i++)
//							std::cout << right << setw(numWidth) << setfill(separator) << fval_[i];
//						std::cout << right << setw(numWidth) << setfill(separator) << Fsq_;
//						std::cout << right << setw(numWidth) << setfill(separator) << sqrt(dd_) << std::endl;
//
//						if (N_ == 2)
//							int test=0;
//					}
//					//std::cout << maxi_ << "\t\t\t" << maxc_ << "\t" << (*x_)(0, 0) << "\t" << (*x_)(1, 0) << "\t" << (*fval_)(0, 0) << "\t" << (*fval_)(1, 0) << "\t" << Fsq_ << "\t" << sqrt(dd_) << std::endl;
//					switch (is_) {
//					case 1:
//					
//						reviseDelta();
//						break;
//					case 2:
//					
//						reviseJac();
//						break;
//					case 3:
//						evalJacobian1(x_, jac_x_, Jac_, N_, N_, _fun, fval_, jac_fval_, &dstep_);
//						inverseMatrixArray(Jac_, invJac_,invfield_,N_);
//						if (debug_) {
//							std::cout << "J = " << std::endl;
//							dispMatrix(Jac_, N_, N_);
//							std::cout << "J-1 = " << std::endl;
//							dispMatrix(invJac_, N_, N_);
//						}
//						maxc_ += N_;
//						mainstep();
//						break;
//					case 4:
//						if (Fsq_ > Fmin_) {
//							interchange();
//						}
//						else {
//						
//							is_ = 2;
//							setSpecialStep();
//						}
//						break;
//					case 5:
//					
//						Fmin_ = Fsq_;
//						for (int i = 0; i < N_; i++){
//							lastx_[i] = x_[i];
//							lastfval_[i] = fval_[i];
//						}
//						is_ = 3;
//						break;
//					}
//				}
//				else {
//					if (maxc_ >= maxfun_) {
//						if (dispOpt_ != 0)
//							std::cout << "Powell terminated: maximum number of function evaluations reached." << std::endl;
//						run_ = 0;
//						outflg_ = 2;
//					}
//					else {
//						if (dispOpt_ != 0)
//							std::cout << "Powell terminated: maximum number of Iterations reached." << std::endl;
//						run_ = 0;
//						outflg_ = 1;
//					}
//					if (Fsq_ > Fmin_) {
//						resetStep();
//					}
//				}
//			}
//		}
//		else {
//			if (dispOpt_ == 1) {
//				std::cout << right << setw(numWidth) << setfill(separator) << maxi_;
//				std::cout << right << setw(numWidth) << setfill(separator) << maxc_;
//				for (int i = 0; i < N_; i++)
//					std::cout << right << setw(numWidth) << setfill(separator) << x_[i];
//				for (int i = 0; i < N_; i++)
//					std::cout << right << setw(numWidth) << setfill(separator) << fval_[i];
//				std::cout << right << setw(numWidth) << setfill(separator) << Fsq_;
//				std::cout << right << setw(numWidth) << setfill(separator) << sqrt(dd_) << std::endl;
//			}
//			for (int i = 0; i < N_; i++)
//				solvec_[i] = x_[i];
//			if (isnan(Fsq_) == 1) {
//				if (dispOpt_ != 0)
//					std::cout << "Powell terminated: Function value is a complex number." << std::endl;
//				run_ = 0;
//				outflg_ = 6;
//			}
//			else {
//				if (dispOpt_ != 0)
//					std::cout << "Powell terminated: Solution found." << std::endl;
//				run_ = 0;
//				outflg_ = 0;
//			}
//
//		}
//		if (debug_)
//			std::cin.get();
//	}
//};

void powell::solve(std::vector<double>* initGuess) {
	// initialization of matrices/arrays
	//matrix initGuessm_(N_, 1);
	//initGuess_ = &initGuessm_;
	//matrix xm_(N_, 1);

	for (int i = 0; i < N_; i++)
		initGuessm_[i] = solvec_[i] = (*initGuess)[i];

	x_ = xm_;

	for (int i = 0; i < N_; i++)
		x_[i] = initGuess_[i];

	/*matrix lastxm_(N_, 1);
	lastx_ = &lastxm_;
	matrix tempxm_(N_, 1);
	tempx_ = &tempxm_;
	matrix fvalm_(N_, 1);
	fval_ = &fvalm_;
	matrix lastfvalm_(N_, 1);
	lastfval_ = &lastfvalm_;
	matrix tempfvalm_(N_, 1);
	tempfval_ = &tempfvalm_;
	matrix vm_(N_, 1);
	v_ = &vm_;
	matrix gm_(N_, 1);
	g_ = &gm_;
	matrix stepm_(N_, 1);
	step_ = &stepm_;
	matrix Jacm_(N_, N_);
	Jac_ = &Jacm_;
	matrix invJacm_(N_, N_);
	invJac_ = &invJacm_;
	matrix gtempm_(N_, 1);
	gtemp_ = &gtempm_;
	matrix Wm_(N_, N_);
	W_ = &Wm_;
	matrix phim_(N_, 1);
	(phi_) = &phim_;
	matrix steptransposedm_(1, N_);
	steptransposed_ = &steptransposedm_;
	matrix gammam_(N_, 1);
	gamma_ = &gammam_;
	matrix sigm_(N_, 1);
	sig_ = &sigm_;*/

	initialize();
	if (dispOpt_ != 0)
		std::cout << "Start solving... Requested accuracy = " << acc_ << std::endl;
	if (dispOpt_ == 1) {
		std::cout << right << setw(numWidth) << setfill(separator) << "Iteration";
		std::cout << right << setw(numWidth) << setfill(separator) << "FunEvals";
		for (int i = 0; i < N_; i++)
			//std::cout << right << setw(numWidth) << setfill(separator) << SSTR("x" << i + 1);
			for (int i = 0; i < N_; i++)
				//std::cout << right << setw(numWidth) << setfill(separator) << SSTR("f" << i + 1);
				std::cout << right << setw(numWidth) << setfill(separator) << "F(x)";
		std::cout << right << setw(numWidth) << setfill(separator) << "Delta" << std::endl;
	}
	while (run_) {

		//x_vec.insert(x_vec.begin(), x_[0], x_[N_-1]);
		//fval_vec.insert(fval_vec.begin(), fval_[0], fval_[N_-1]);

		//x_vec.insert(x_vec.begin(), x_, x_ + N_ - 1);
		//fval_vec.insert(fval_vec.begin(), fval_, fval_ + N_ - 1);

		x_vec.assign(x_, x_ + N_);
		fval_vec.assign(fval_, fval_ + N_);


		calfun(&x_vec, &fval_vec, &Fsq_);

		x_ = x_vec.data();
		fval_ = fval_vec.data();

		//x_ = &x_vec[0];
		//fval_ = &fval_vec[0];

		maxc_++;
		maxi_++;
		if (Fsq_ > acc_) {
			if (is_ == 1 || is_ == 4) {
				if (Fsq_ >= Fmin_) {
					if (dd_ - dss_ <= 0) {
						ntest_--;
						if (ntest_ < 0) {
							if (dispOpt_ != 0)
								std::cout << "Powell terminated: failed to decrease Fsq_ using a new Jacobian." << std::endl;
							run_ = 0;
							outflg_ = 4;
							resetStep();
						}
						else if (ntest_ == 0) {
							if (dispOpt_ != 0)
								std::cout << "Powell terminated: " << nt_ << " calls of calfun ineffective." << std::endl;
							run_ = 0;
							outflg_ = 3;
							resetStep();
						}
					}
				}
				else {
					ntest_ = nt_;
				}
			}
			if (run_) {
				if (maxc_ <= maxfun_ && maxi_ <= maxiter_) {
					if (dispOpt_ == 1) {
						std::cout << right << setw(numWidth) << setfill(separator) << maxi_;
						std::cout << right << setw(numWidth) << setfill(separator) << maxc_;
						for (int i = 0; i < N_; i++)
							std::cout << right << setw(numWidth) << setfill(separator) << x_[i];
						for (int i = 0; i < N_; i++)
							std::cout << right << setw(numWidth) << setfill(separator) << fval_[i];
						std::cout << right << setw(numWidth) << setfill(separator) << Fsq_;
						std::cout << right << setw(numWidth) << setfill(separator) << sqrt(dd_) << std::endl;

						if (N_ == 2)
							int test = 0;
					}
					//std::cout << maxi_ << "\t\t\t" << maxc_ << "\t" << (*x_)(0, 0) << "\t" << (*x_)(1, 0) << "\t" << (*fval_)(0, 0) << "\t" << (*fval_)(1, 0) << "\t" << Fsq_ << "\t" << sqrt(dd_) << std::endl;
					switch (is_) {
					case 1:

						reviseDelta();
						break;
					case 2:

						reviseJac();
						break;
					case 3:
						tempx_vec.assign(jac_x_, jac_x_ + N_);
						tempfval_vec.assign(jac_fval_, jac_fval_ + N_);

						evalJacobian1(x_, &tempx_vec, Jac_, N_, N_, _fun, fval_, &tempfval_vec, &dstep_);
						inverseMatrixArray(Jac_, invJac_, invfield_, N_);
						if (debug_) {
							std::cout << "J = " << std::endl;
							dispMatrix(Jac_, N_, N_);
							std::cout << "J-1 = " << std::endl;
							dispMatrix(invJac_, N_, N_);
						}
						maxc_ += N_;
						mainstep();
						break;
					case 4:
						if (Fsq_ > Fmin_) {
							interchange();
						}
						else {

							is_ = 2;
							setSpecialStep();
						}
						break;
					case 5:

						Fmin_ = Fsq_;
						for (int i = 0; i < N_; i++) {
							lastx_[i] = x_[i];
							lastfval_[i] = fval_[i];
						}
						is_ = 3;
						break;
					}
				}
				else {
					if (maxc_ >= maxfun_) {
						if (dispOpt_ != 0)
							std::cout << "Powell terminated: maximum number of function evaluations reached." << std::endl;
						run_ = 0;
						outflg_ = 2;
					}
					else {
						if (dispOpt_ != 0)
							std::cout << "Powell terminated: maximum number of Iterations reached." << std::endl;
						run_ = 0;
						outflg_ = 1;
					}
					if (Fsq_ > Fmin_) {
						resetStep();
					}
				}
			}
		}
		else {
			if (dispOpt_ == 1) {
				std::cout << right << setw(numWidth) << setfill(separator) << maxi_;
				std::cout << right << setw(numWidth) << setfill(separator) << maxc_;
				for (int i = 0; i < N_; i++)
					std::cout << right << setw(numWidth) << setfill(separator) << x_[i];
				for (int i = 0; i < N_; i++)
					std::cout << right << setw(numWidth) << setfill(separator) << fval_[i];
				std::cout << right << setw(numWidth) << setfill(separator) << Fsq_;
				std::cout << right << setw(numWidth) << setfill(separator) << sqrt(dd_) << std::endl;
			}
			for (int i = 0; i < N_; i++)
				solvec_[i] = x_[i];
			if (isnan(Fsq_) == 1) {
				if (dispOpt_ != 0)
					std::cout << "Powell terminated: Function value is a complex number." << std::endl;
				run_ = 0;
				outflg_ = 6;
			}
			else {
				if (dispOpt_ != 0)
					std::cout << "Powell terminated: Solution found." << std::endl;
				run_ = 0;
				outflg_ = 0;
			}

		}
		if (debug_)
			std::cin.get();
	}
};

void powell::initialize() {

	for (int i = 0; i < N_; i++) {


		lastx_[i] = 0;
		tempx_[i] = 0;
		fval_[i] = 0;
		lastfval_[i] = 0;
		tempfval_[i] = 0;
		v_[i] = 0;
		jac_x_[i] = 0;
		jac_fval_[i] = 0;
		g_[i] = 0;
		step_[i] = 0;
		gtemp_[i] = 0;

		w_[i] = N_ - i;
		W_[i * N_ + i] = 1;
		phi_[i] = 0;
		step_[i] = 0;
	}

	// doubles
	tinc_ = 1;
	//acc_ = 1e-6;
	/*dmax_ = 10.0;*/
	dm_ = dmax_ * dmax_;
	dmm_ = 4 * dm_;
	/*dstep_ = (1e-6);*/
	dd_ = 0;
	dss_ = dstep_ * dstep_;
	dtest_ = 2 * double(N_) - 0.5;
	Fsq_ = 0;
	mu_ = 0;
	sp_ = 0;
	ds_ = 0;
	dn_ = 0;
	gm_ = 0;
	anmult_ = 0;
	dmult_ = 0;
	Phi_ = 0;
	tinc_ = 1;
	ss_ = 0;
	pj_ = 0;
	alpha_ = 0;
	sval_ = 0;

	// bools/ints
	nt_ = N_ + 4;
	ntest_ = nt_;
	maxc_ = 0;
	//maxiter_ = 1000;
	maxi_ = 0;
	//maxfun_ = 10000;
	is_ = 5;
	run_ = 1;
	flag_ = 0;
	is_ = 5;
	kk_ = 0;
	outflg_ = -1;
	// outflg_ = 0 -> solution found
	// outflg_ = 1 -> maxiter_ reached
	// outflg_ = 2 -> maxfun_ reached
	// outflg_ = 3 -> nt_ steps ineffective
	// outflg_ = 4 -> failed to decrease Fsq_ using a new Jacobian
	// outflg_ = 5 -> nearby stationary point predicted
	// outflg_ = 6 -> function value is a complex number

	separator = ' ';
	nameWidth = 6;
	numWidth = 14;
};

//void powell::calfun(double* x, double* f, double* Fsq) { // returns the residual sum of squares of the function
//	
//	/*x_vec(&x[0], &x[N_]);
//	fval_vec(&f[0], &f[N_]);*/
//
//
//	(*_fun).Zero_Functionction(x_vec, fval_vec);
//
//	(*Fsq) = 0;
//	for (int i = 0; i < N_; i++) {
//		(*Fsq) += f[i] * f[i];
//	};
//	/*int j = -1;
//	int i = -1;
//	for (int k = 0; k < (*f).sizeA(); k++) {
//	j++;
//	if (!(j % (*f).cols())) {
//	i++;
//	j = 0;
//	}
//	(*Fsq) += (*f)(i,j) * (*f)(i,j);
//	};*/
//	return;
//};	
void powell::calfun(std::vector<double>* x, std::vector<double>* f, double* Fsq) { // returns the residual sum of squares of the function

	/*x_vec(&x[0], &x[N_]);
	fval_vec(&f[0], &f[N_]);*/


	(*_fun).zero_function(x, f);

	(*Fsq) = 0.0;
	for (int i = 0; i < N_; i++) {
		(*Fsq) += (*f)[i] * (*f)[i];
	};
	/*int j = -1;
	int i = -1;
	for (int k = 0; k < (*f).sizeA(); k++) {
	j++;
	if (!(j % (*f).cols())) {
	i++;
	j = 0;
	}
	(*Fsq) += (*f)(i,j) * (*f)(i,j);
	};*/
	return;
};

void powell::reviseDelta() {
	if (debug_)
		std::cout << "reviseDelta(), is_ = " << is_ << std::endl;
	dmult_ = 0.9 * Fmin_ + 0.1 * Phi_;
	if (Fsq_ > dmult_) {
		dd_ = max(0.25 * dd_, dss_);
		tinc_ = 1;
		if (Fsq_ >= Fmin_) {
			reviseJac();
		}
		else {
			interchange();
		}
	}
	else {
		dmult_ = dmult_ - Fsq_;
		sp_ = 0; ds_ = 0;  ss_ = 0;
		for (int i = 0; i < N_; i++) {
			sp_ += fabs(fval_[i] * (fval_[i] - phi_[i]));
			ds_ = (fval_[i] - phi_[i]);
			ss_ += ds_ * ds_;
		}
		pj_ = 1 + dmult_ / (sp_ + sqrt(sp_ * sp_ + dmult_ * ss_));
		pj_ = sqrt(pj_);
		mu_ = min(2.0, min(pj_, tinc_));
		tinc_ = pj_ / mu_;
		dd_ = min(mu_ * sqrt(dd_), dmax_);
		dd_ *= dd_;
		interchange();
	}
	if (debug_)
		std::cout << "new Delta = " << sqrt(dd_) << std::endl;
	if (N_ == 2)
		int test = 0;
};

void powell::reviseJac() {

	/*	(*gamma_) = (*fval_) - (lastfval_);
		(*step_) = (*x_) - (lastx_);*/
		/*(*step_).mat_transpose(steptransposed_);*/ //7

	for (int j = 0; j < N_; j++) {
		gamma_[j] = fval_[j] - lastfval_[j];
		step_[j] = x_[j] - lastx_[j];
	}

	if (debug_) {
		std::cout << "reviseJac(), is_ = " << is_ << std::endl;
		std::cout << "gamma = " << std::endl;
		for (int i = 0; i < N_; i++)
			std::cout << gamma_[i] << std::endl;
	}
	//matrix singtest_mat(1, 1);
	//singtest_mat = (*steptransposed_) * &((*invJac_) * (gamma_));
	//double singtest_;
	//singtest_ = fabs(singtest_mat(0, 0));
	//ds_ = (*step_).euclidMetric(); //8

	ds_ = 0;
	double singtest = 0;
	double singtest_ = 0;
	for (int j = 0; j < N_; j++) {
		for (int i = 0; i < N_; i++)
			singtest += step_[j] * (invJac_[j * N_ + i] * gamma_[i]);
		ds_ += pow(step_[j], 2);
	}
	singtest_ = fabs(singtest);

	if (singtest_ >= 0.1 * ds_) {
		alpha_ = 1;
	}
	else {
		alpha_ = 0.8;
	};

	/*(*Jac_) = (*Jac_) + &(alpha_ * ((*gamma_) - &((*Jac_)*(step_)))*&((*steptransposed_) / ds_)); */
	/*(*invJac_) = (*invJac_) + &(alpha_ * ((((*step_) - &((*invJac_) * (gamma_))) * (steptransposed_)) * (invJac_)) / (alpha_*singtest_ + (1 - alpha_)*ds_));*/ //9

	for (int j = 0; j < N_; j++) {

		double a = 0;
		for (int k = 0; k < N_; k++) {
			a += Jac_[j * N_ + k] * step_[k];
		}

		for (int i = 0; i < N_; i++) {
			Jac_[j * N_ + i] = Jac_[j * N_ + i] + (alpha_ * (gamma_[j] - a) * (step_[i] / ds_));
		}

	}

	for (int j = 0; j < N_; j++) {

		double a = 0;
		for (int k = 0; k < N_; k++) {
			a += (invJac_[j * N_ + k] * gamma_[k]);
		}

		for (int i = 0; i < N_; i++) {
			b_[j * N_ + i] = 0;
			for (int k = 0; k < N_; k++) {
				b_[j * N_ + i] += ((step_[j] - a) * step_[k]) * invJac_[k * N_ + i];
			}
		}
	}

	for (int j = 0; j < N_; j++) {
		for (int i = 0; i < N_; i++) {
			invJac_[j * N_ + i] = invJac_[j * N_ + i] + (alpha_ * b_[j * N_ + i] / (alpha_ * singtest_ + (1 - alpha_) * ds_));
		}
	}

	if (debug_) {
		std::cout << "alpha = " << alpha_ << std::endl;
		std::cout << "J = " << std::endl;
		dispMatrix(Jac_, N_, N_);
		std::cout << "J-1 = " << std::endl;
		dispMatrix(invJac_, N_, N_);
	}

	mainstep();

};

void powell::resetStep() {
	if (debug_)
		std::cout << "resetStept(), is_ = " << is_ << std::endl;
	*x_ = *lastx_;
	*fval_ = *lastfval_;
	Fsq_ = Fmin_;
};

void powell::interchange() {
	if (debug_)
		std::cout << "interchange(), is_ = " << is_ << std::endl;
	Fmin_ = Fsq_;
	for (int i = 0; i < N_; i++) {
		tempx_[i] = x_[i];
		x_[i] = lastx_[i];
		lastx_[i] = tempx_[i];
		tempfval_[i] = fval_[i];
		fval_[i] = lastfval_[i];
		lastfval_[i] = tempfval_[i];
	}
	/*(*phi_) = &((-1)*(*phi_));*/ //10

	for (int j = 0; j < N_; j++)
		phi_[j] = (-1) * phi_[j];

	if (is_ <= 1) {
		reviseJac();
	}
	else {
		is_ = 2;
		setSpecialStep();
	}
};

void powell::setSpecialStep() {
	if (debug_)
		std::cout << "setSpecialStep(), is_ = " << is_ << std::endl;
	for (int i = 0; i < N_; i++) {
		step_[i] = dstep_ * W_[i * N_ + 0];
		if (i < (N_ - 1))
			w_[i] = 1 + w_[i + 1];
		else
			w_[N_ - 1] = 1;
		for (int j = 0; j < N_; j++) {
			if (j < (N_ - 1))
				W_[i * N_ + j] = W_[i * N_ + j + 1];
			else
				W_[i * N_ + j] = step_[i] / 2 / ds_;
		}
	}
};

void powell::mainstep() {
	if (debug_)
		std::cout << "mainstep(), is_ = " << is_ << std::endl;
	for (int i = 0; i < N_; i++) {
		g_[i] = 0;
		for (int j = 0; j < N_; j++)
			g_[i] -= (Jac_[j * N_ + i] * lastfval_[j]);
	}

	/*(*v_) = (-1)*(*invJac_)*(lastfval_);*/ //1

	for (int j = 0; j < N_; j++) {
		v_[j] = 0;
	}

	for (int j = 0; j < N_; j++) {
		for (int i = 0; i < N_; i++)
			v_[j] += (-1) * (invJac_[j * N_ + i] * lastfval_[i]);
	}

	/*	ds_ = (*g_).euclidMetric();
		dn_ = (*v_).euclidMetric();*/ //11
	ds_ = 0;
	dn_ = 0;
	for (int i = 0; i < N_; i++) {
		ds_ += pow(g_[i], 2);
		dn_ += pow(v_[i], 2);
	}
	sp_ = 0;
	for (int i = 0; i < N_; i++) {
		sp_ += g_[i] * v_[i];
	}

	if (Fmin_ * Fmin_ - dmm_ * ds_ <= 0) {
		is_ = 2;
		if (dn_ > dd_) {
			flag_ = 1;
			while (flag_ == 1) {
				//std::cout << "flag_ = " << flag_ << std::endl;

				/*(*gtemp_) = (*Jac_)*(g_);*/ //2
				/*dmult_ = (*gtemp_).euclidMetric();*/ //3

				dmult_ = 0;
				for (int j = 0; j < N_; j++) {
					gtemp_[j] = 0;
				}
				for (int j = 0; j < N_; j++) {
					for (int i = 0; i < N_; i++)
						gtemp_[j] += (Jac_[j * N_ + i] * g_[i]);
					dmult_ += pow(gtemp_[j], 2);
				}

				dmult_ = ds_ / dmult_; //mu
				ds_ = ds_ * dmult_ * dmult_;
				if (ds_ - dd_ < 0) {
					//std::cout << "setstep3, ds_ =" << ds_ << std::endl;
					flag_ = 0;
					sp_ = sp_ * dmult_;
					if ((sp_ - dd_) * (sp_ - dd_) + (dn_ - dd_) * (dd_ - ds_) < 0)
						anmult_ = 0;
					else
						anmult_ = (dd_ - ds_) / ((sp_ - ds_) + sqrt((sp_ - dd_) * (sp_ - dd_) + (dn_ - dd_) * (dd_ - ds_)));

					if (isinf(anmult_) == 1) {
						anmult_ = 0;
					}

					dmult_ = dmult_ * (1 - anmult_);

				}
				else {
					if (dd_ <= 0) {
						dd_ = max(dss_, min(dm_, ds_));
						ds_ = ds_ / (dmult_ * dmult_);
						//std::cout << "initial dd_ = " << sqrt(dd_) << std::endl;
					}
					else {
						flag_ = 0;
						//std::cout << "setstep2" << std::endl;
						anmult_ = 0;
						dmult_ = dmult_ * sqrt(dd_ / ds_);
					}
				}
				if (flag_ == 0) {

					/*(*step_) = dmult_ * (*g_) + &(anmult_ * (*v_));*/ //4

					for (int j = 0; j < N_; j++)
						step_[j] = dmult_ * g_[j] + anmult_ * v_[j];

					if (debug_) {
						double steplength = 0;
						std::cout << "step = " << std::endl;
						for (int j = 0; j < N_; j++) {
							std::cout << step_[j] << std::endl;
							steplength += pow(step_[j], 2);
						}
						steplength = pow(steplength, 0.5);
						std::cout << "steplength = " << steplength << std::endl;
					}

					/*dn_ = (*step_).euclidMetric();*/ //5

					dn_ = 0;
					for (int j = 0; j < N_; j++)
						dn_ += pow(step_[j], 2);

					ds_ = 0.25 * dn_;
					sp_ = 0;
					for (int i = 0; i < N_; i++) {
						sp_ += step_[i] * W_[i * N_ + 0];
					}
					if (w_[0] - dtest_ <= 0) {
						updateOmega();
					}
					else {
						if (sp_ * sp_ - ds_ >= 0) {
							updateOmega();
						}
						else {
							setSpecialStep();
						}
					}
				}
			}
		}
		else {
			//std::cout << "setstep1" << std::endl;
			for (int i = 0; i < N_; i++)
				step_[i] = v_[i];
			if (debug_) {
				double steplength = 0;
				std::cout << "step = " << std::endl;
				for (int j = 0; j < N_; j++) {
					std::cout << step_[j] << std::endl;
					steplength += pow(step_[j], 2);
				}
				steplength = pow(steplength, 0.5);
				std::cout << "steplength = " << steplength << std::endl;
			}
			dd_ = max(dn_, dss_);
			ds_ = 0.25 * dn_;
			tinc_ = 1;
			if (dn_ - dss_ >= 0) {
				updateOmega();
			}
			else {
				is_ = 4;
				setNewPoint();
			}
		}
	}
	else {
		if (is_ == 1 || is_ == 2) {
			ntest_ = 0;
			(x_) = (lastx_);
			is_ = 3;
		}
		else {
			if (dispOpt_ != 0)
				std::cout << "Powell terminated: nearby stationary point predicted." << std::endl;
			run_ = 0;
			outflg_ = 5;
			resetStep();
		}
	}
};

void powell::updateOmega() {
	double tempVal;// temporary value for swapping of the columns
	if (debug_)
		std::cout << "updateOmega(), is_ = " << is_ << std::endl;
	sp_ = 0;
	for (int i = 0; i < N_; i++) {
		a_[i] = 0;
		for (int j = 0; j < N_; j++) {
			a_[i] += step_[j] * W_[j * N_ + i];
		}
		if (is_ == 1) {
			w_[i] = w_[i + 1] + 1;
		}
		else {
			w_[i] = 1 + w_[i];
			sp_ += a_[i] * a_[i];
			if (sp_ - ds_ <= 0) //*ds_
				;//Wflag_ = 0;
			else {
				is_ = 1;
				kk_ = i;
			}
		}
	}

	for (int i = 0; i < N_ - 1; i++) {
		if (i < kk_)
			w_[i] = 1 + w_[i];
		else
			w_[i] = w_[i + 1] + 1;
	}
	w_[N_ - 1] = 1;

	for (int i = kk_; i > 0; i--) {// swapping of the columns i and i-1
		for (int sc = 0; sc < N_; sc++) {
			tempVal = W_[sc * N_ + i];
			W_[sc * N_ + i] = W_[sc * N_ + i - 1];
			W_[sc * N_ + i - 1] = tempVal;
		}
	}
	sval_ = a_[0] * a_[0];
	for (int j = 0; j < N_; j++)
		sig_[j] = 0;
	for (int i = 1; i < N_; i++) {
		for (int j = 0; j < N_; j++) {
			sig_[j] = sig_[j] + a_[i - 1] * W_[j * N_ + i - 1];
			W_[j * N_ + i - 1] = (sval_ * W_[j * N_ + i] - a_[i] * sig_[j]) / (sqrt(sval_ * (sval_ + a_[i] * a_[i])));
		}
		sval_ += a_[i] * a_[i];
	}
	for (int i = 0; i < N_; i++)
		W_[i * N_ + N_ - 1] = step_[i] / 2 / ds_;
	if (debug_) {
		std::cout << "W = " << std::endl;
		dispMatrix(W_, N_, N_);
	}
	setNewPoint();
};

void powell::setNewPoint() {
	if (debug_)
		std::cout << "setNewPoint(), is_ = " << is_ << std::endl;

	/*(*x_) = (*lastx_) + (step_);
	(*phi_) = (*lastfval_) + &((*Jac_) * (step_));*/
	/*Phi_ = (*phi_).euclidMetric();*/ //6

	Phi_ = 0;
	for (int j = 0; j < N_; j++) {
		phi_[j] = 0;
	}
	for (int j = 0; j < N_; j++) {
		x_[j] = lastx_[j] + step_[j];

		for (int i = 0; i < N_; i++)
			phi_[j] += Jac_[j * N_ + i] * step_[i];

		phi_[j] += lastfval_[j];
		Phi_ += pow(phi_[j], 2);
	}

};

//matrix powell::solution() {
//	matrix sol(N_, 1);
//	if (debug_)
//		std::cout << "Solution(), is_ = " << is_ << std::endl;
//	for (int i = 0; i < N_; i++)
//		sol(i, 0) = solvec_[i];
//	return sol;
//};

std::vector<double> powell::solution() {
	std::vector<double> sol(N_, 0);
	if (debug_)
		std::cout << "Solution(), is_ = " << is_ << std::endl;
	for (int i = 0; i < N_; i++)
		sol[i] = solvec_[i];
	return sol;
};

double powell::solution_function() {
	return Fsq_;
}

//void powell::solution(matrix *solution) {
//	if (debug_)
//		std::cout << "Solution(), is_ = " << is_ << std::endl;
//	for (int i = 0; i < N_; i++)
//		(*solution)(i, 0) = solvec_[i];
//	return;
//};

void powell::solution(std::vector<double>* solution) {
	if (debug_)
		std::cout << "Solution(), is_ = " << is_ << std::endl;
	for (int i = 0; i < N_; i++)
		(*solution)[i] = solvec_[i];
	return;
};

int powell::outflg() {

	if (debug_)
		std::cout << "outflg(), is_ = " << is_ << std::endl;
	return outflg_;
};

//void evalJacobian1(double *vec, double *temp_vec, matrix *Jacobian, int n, int m, Zero_Function *A, double *fval, double *temp_fval, double* dstep)
//{
//	double dx = *dstep;
//
//	//double *temp_vec_test;
//	//double *temp_fval_test;
//
//	//temp_vec_test = new double[n];
//	//temp_fval_test = new double[n];
//
//	for (int j = 0; j < m; j++) 
//		temp_vec[j] = vec[j];
//
//	for (int i = 0; i < n; i++) {
//		temp_vec[i] = vec[i] + dx;
//		(*A).Zero_Functionction(temp_vec, temp_fval);
//		for (int j = 0; j < m; j++) {
//			(*Jacobian)(j, i) = (temp_fval[j] - fval[j]) / dx;
//		}
//		temp_vec[i] = vec[i];
//	}
//
//
//	//if (n == 2)
//	//	std::cout << "test" << std::endl;
//
//	return;
//};

void evalJacobian1(double* vec, std::vector<double>* temp_vec, double* Jacobian, int n, int m, Zero_Function* A, double* fval, std::vector<double>* temp_fval, double* dstep)
{
	double dx = *dstep;

	//double *temp_vec_test;
	//double *temp_fval_test;

	//temp_vec_test = new double[n];
	//temp_fval_test = new double[n];


	for (int j = 0; j < m; j++)
		(*temp_vec)[j] = vec[j];

	for (int i = 0; i < n; i++) {
		(*temp_vec)[i] = vec[i] + dx;
		(*A).zero_function(temp_vec, temp_fval);
		for (int j = 0; j < m; j++) {
			Jacobian[j * n + i] = ((*temp_fval)[j] - fval[j]) / dx;
		}
		(*temp_vec)[i] = vec[i];
	}


	//if (n == 2)
	//	std::cout << "test" << std::endl;

	return;
};

void inverseMatrixArray(double* mat, double* invMat, double* invField, int dim) {

	int is;
	double tempval; // temporary value used for row swap and row addition in the Gauß Jordan Algorithm
	double det = 1.0;


	// Creating the field for the inversion
	for (int i = 0; i < dim; ++i) // loop over the rows
	{
		for (int j = 0; j < dim; ++j) // loop over first half of columns
			invField[i * 2 * dim + j] = mat[i * dim + j];
		for (int j = dim; j < 2 * dim; ++j) // loop over second half of columns
			invField[i * 2 * dim + j] = (i == j - dim) ? 1.0 : 0.0;
	}


	for (int i = 0; i < dim - 1; i++) { // loop over the rows
		is = i;
		// row swap 
		bool swap = 0;
		//for (int j = is + 1; j < dim; j++) {
		//	if (fabs(mat[j*dim + i]) > fabs(mat[i*dim+i])) {
		//		is = j;
		//		swap = 1;
		//	}
		//}

		//// if swap condition is meet, rows are swapped
		//if (swap) {
		//	for (int k = 0; k < 2 * dim; k++) {
		//		tempval = invField[i * 2 * dim + k];
		//		invField[i * 2 * dim + k] = invField[is * 2 * dim + k];
		//		invField[is * 2 * dim + k] = tempval;
		//	}
		//	
		//}

		for (int j = i + 1; j < dim; j++) {
			tempval = invField[j * 2 * dim + i] / invField[i * 2 * dim + i];
			for (int k = i; k < 2 * dim; k++) {
				invField[j * 2 * dim + k] -= tempval * invField[i * 2 * dim + k];
			}
		}
	}

	// Calculationg the determinant of the first half

	for (int i = 0; i < dim; i++)
		det += invField[i * 2 * dim + i];

	if (det != 0.0) {
		for (int i = dim - 1; i > 0; i--) {
			for (int j = i - 1; j >= 0; j--) {
				tempval = invField[j * 2 * dim + i] / invField[i * 2 * dim + i];
				for (int k = i; k < 2 * dim; k++) {
					invField[j * 2 * dim + k] -= tempval * invField[i * 2 * dim + k];
				}
			}
		}
		for (int i = 0; i < dim; i++) {
			tempval = invField[i * 2 * dim + i];
			for (int j = 0; j < dim; j++) {
				invMat[i * dim + j] = invField[i * 2 * dim + j + dim] / tempval;
			}
		}
	}
	else {
		std::cout << "A nicht invertierbar";
	};



}

void dispMatrix(double* mat, int col, int row) {
	// displays the matrix in Array form in the console

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if ((j == 0) && (j == col - 1))
				std::cout << "(" << "\t" << mat[i * col + j] << "\t" << ")" << std::endl;
			else if (j == 0)
				std::cout << "(" << "\t" << mat[i * col + j] << "\t";
			else if (j == col - 1)
				std::cout << mat[i * col + j] << "\t" << ")" << std::endl;
			else
				std::cout << mat[i * col + j] << "\t";
		}
	}
}

//one-dimensional powell

	// ctor
powell_1D::powell_1D() {
	getVal_funcCtr_ = 0;
};

powell_1D::powell_1D(Zero_Function* fun, double* initGuess, options_powell* opts) :
	_fun(fun),
	dispOpt_((*opts).dispOpt),
	acc_((*opts).acc),
	maxiter_((*opts).maxIter),
	dstep_((*opts).dstep),
	dmax_((*opts).dmax),
	maxfun_((*opts).maxFunEvals) {
	N_ = 1;
	debug_ = 0;
	solvec_ = (*initGuess);
	getVal_funcCtr_ = 1;
};

void powell_1D::get_Vals(Zero_Function* fun, double* initGuess, options_powell* opts) {

	_fun = fun;
	solvec_ = (*initGuess);
	dispOpt_ = (*opts).dispOpt;
	acc_ = (*opts).acc;
	dstep_ = (*opts).dstep;
	maxiter_ = (*opts).maxIter;
	maxfun_ = (*opts).maxFunEvals;
	debug_ = 0;
	getVal_funcCtr_++;
};

void powell_1D::solve() {
	// initialization of matrices/arrays
	//double initGuessm_(N_, 1);
	//initGuess_ = &initGuessm_;
	//double xm_(N_, 1);

	//x_ = &xm_;
	x_ = initGuess_ = solvec_;
	//double lastxm_(N_, 1);
	//lastx_ = &lastxm_; 
	//double tempxm_(N_, 1); 
	//tempx_ = &tempxm_;
	/*double fvalm_(N_, 1);
	fval_ = &fvalm_;
	double lastfvalm_(N_, 1);
	lastfval_ = &lastfvalm_;
	double tempfvalm_(N_, 1);
	tempfval_ = &tempfvalm_;
	double vm_(N_, 1);
	v_ = &vm_;
	double gm_(N_, 1);
	g_ = &gm_;
	double stepm_(N_, 1);
	step_ = &stepm_;
	double Jacm_(N_, N_);
	Jac_ = &Jacm_;
	double invJacm_(N_, N_);
	invJac_ = &invJacm_;
	double gtempm_(N_, 1);
	gtemp_ = &gtempm_;
	double Wm_(N_, N_);
	W_ = &Wm_;
	double phim_(N_, 1);
	(phi_) = &phim_;
	double steptransposedm_(1, N_);
	steptransposed_ = &steptransposedm_;
	double gammam_(N_, 1);
	gamma_ = &gammam_;
	double sigm_(N_, 1);
	sig_ = &sigm_;*/

	initialize();
	if (dispOpt_ != 0)
		std::cout << "Start solving... Requested accuracy = " << acc_ << std::endl;
	if (dispOpt_ == 1) {
		std::cout << right << setw(numWidth) << setfill(separator) << "Iteration";
		std::cout << right << setw(numWidth) << setfill(separator) << "FunEvals";
		//std::cout << right << setw(numWidth) << setfill(separator) << SSTR("x" << 1);
		//std::cout << right << setw(numWidth) << setfill(separator) << SSTR("f" << 1);
		std::cout << right << setw(numWidth) << setfill(separator) << "F(x)";
		std::cout << right << setw(numWidth) << setfill(separator) << "Delta" << std::endl;
	}
	while (run_) {
		calfun(&x_, &fval_, &Fsq_);
		maxc_++;
		maxi_++;
		if (Fsq_ > acc_) {
			if (is_ == 1 || is_ == 4) {
				if (Fsq_ >= Fmin_) {
					if (dd_ - dss_ <= 0) {
						ntest_--;
						if (ntest_ < 0) {
							if (dispOpt_ != 0)
								std::cout << "Powell terminated: failed to decrease Fsq_ using a new Jacobian." << std::endl;
							run_ = 0;
							outflg_ = 4;
							resetStep();
						}
						else if (ntest_ == 0) {
							if (dispOpt_ != 0)
								std::cout << "Powell terminated: " << nt_ << " calls of calfun ineffective." << std::endl;
							run_ = 0;
							outflg_ = 3;
							resetStep();
						}
					}
				}
				else {
					ntest_ = nt_;
				}
			}
			if (run_) {
				if (maxc_ <= maxfun_ && maxi_ <= maxiter_) {
					if (dispOpt_ == 1) {
						std::cout << right << setw(numWidth) << setfill(separator) << maxi_;
						std::cout << right << setw(numWidth) << setfill(separator) << maxc_;
						std::cout << right << setw(numWidth) << setfill(separator) << x_;
						std::cout << right << setw(numWidth) << setfill(separator) << fval_;
						std::cout << right << setw(numWidth) << setfill(separator) << Fsq_;
						std::cout << right << setw(numWidth) << setfill(separator) << sqrt(dd_) << std::endl;
					}
					//std::cout << maxi_ << "\t\t\t" << maxc_ << "\t" << (*x_)(0, 0) << "\t" << (*x_)(1, 0) << "\t" << (*fval_)(0, 0) << "\t" << (*fval_)(1, 0) << "\t" << Fsq_ << "\t" << sqrt(dd_) << std::endl;
					switch (is_) {
					case 1:
						reviseDelta();
						break;
					case 2:
						reviseJac();
						break;
					case 3:
						evalJacobian1(&x_, &Jac_, _fun, &fval_, &dstep_);
						invJac_ = 1 / Jac_;
						if (debug_) {
							std::cout << "J = " << std::endl;
							std::cout << Jac_ << std::endl;
							std::cout << "J-1 = " << std::endl;
							std::cout << invJac_ << std::endl;
						}
						maxc_ += N_;
						mainstep();
						break;
					case 4:
						if (Fsq_ > Fmin_) {
							interchange();
						}
						else {
							is_ = 2;
							setSpecialStep();
						}
						break;
					case 5:
						Fmin_ = Fsq_;
						lastx_ = x_;
						lastfval_ = fval_;
						is_ = 3;
						break;
					}
				}
				else {
					if (maxc_ >= maxfun_) {
						if (dispOpt_ != 0)
							std::cout << "Powell terminated: maximum number of function evaluations reached." << std::endl;
						run_ = 0;
						outflg_ = 2;
					}
					else {
						if (dispOpt_ != 0)
							std::cout << "Powell terminated: maximum number of Iterations reached." << std::endl;
						run_ = 0;
						outflg_ = 1;
					}
					if (Fsq_ > Fmin_) {
						resetStep();
					}
				}
			}
		}
		else {
			if (dispOpt_ == 1) {
				std::cout << right << setw(numWidth) << setfill(separator) << maxi_;
				std::cout << right << setw(numWidth) << setfill(separator) << maxc_;
				std::cout << right << setw(numWidth) << setfill(separator) << x_;
				std::cout << right << setw(numWidth) << setfill(separator) << fval_;
				std::cout << right << setw(numWidth) << setfill(separator) << Fsq_;
				std::cout << right << setw(numWidth) << setfill(separator) << sqrt(dd_) << std::endl;
			}
			solvec_ = x_;
			if (isnan(Fsq_) == 1) {
				if (dispOpt_ != 0)
					std::cout << "Powell terminated: Function value is a complex number." << std::endl;
				run_ = 0;
				outflg_ = 6;
			}
			else {
				if (dispOpt_ != 0)
					std::cout << "Powell terminated: Solution found." << std::endl;
				run_ = 0;
				outflg_ = 0;
			}

		}
		if (debug_)
			std::cin.get();
	}
};

void powell_1D::solve(double* initGuess) {
	// initialization of matrices/arrays
	//double initGuessm_(N_, 1);
	//initGuess_ = &initGuessm_;
	//double xm_(N_, 1);

	x_ = solvec_ = (*initGuess);

	/*double lastxm_(N_, 1);
	lastx_ = &lastxm_;
	double tempxm_(N_, 1);
	tempx_ = &tempxm_;
	double fvalm_(N_, 1);
	fval_ = &fvalm_;
	double lastfvalm_(N_, 1);
	lastfval_ = &lastfvalm_;
	double tempfvalm_(N_, 1);
	tempfval_ = &tempfvalm_;
	double vm_(N_, 1);
	v_ = &vm_;
	double gm_(N_, 1);
	g_ = &gm_;
	double stepm_(N_, 1);
	step_ = &stepm_;
	double Jacm_(N_, N_);
	Jac_ = &Jacm_;
	double invJacm_(N_, N_);
	invJac_ = &invJacm_;
	double gtempm_(N_, 1);
	gtemp_ = &gtempm_;
	double Wm_(N_, N_);
	W_ = &Wm_;
	double phim_(N_, 1);
	(phi_) = &phim_;
	double steptransposedm_(1, N_);
	steptransposed_ = &steptransposedm_;
	double gammam_(N_, 1);
	gamma_ = &gammam_;
	double sigm_(N_, 1);
	sig_ = &sigm_;*/

	initialize();
	if (dispOpt_ != 0)
		std::cout << "Start solving... Requested accuracy = " << acc_ << std::endl;
	if (dispOpt_ == 1) {
		std::cout << right << setw(numWidth) << setfill(separator) << "Iteration";
		std::cout << right << setw(numWidth) << setfill(separator) << "FunEvals";
		//std::cout << right << setw(numWidth) << setfill(separator) << SSTR("x" << 1);
		//std::cout << right << setw(numWidth) << setfill(separator) << SSTR("f" << 1);
		std::cout << right << setw(numWidth) << setfill(separator) << "F(x)";
		std::cout << right << setw(numWidth) << setfill(separator) << "Delta" << std::endl;
	}
	while (run_) {
		calfun(&x_, &fval_, &Fsq_);
		maxc_++;
		maxi_++;
		if (Fsq_ > acc_) {
			if (is_ == 1 || is_ == 4) {
				if (Fsq_ >= Fmin_) {
					if (dd_ - dss_ <= 0) {
						ntest_--;
						if (ntest_ < 0) {
							if (dispOpt_ != 0)
								std::cout << "Powell terminated: failed to decrease Fsq_ using a new Jacobian." << std::endl;
							run_ = 0;
							outflg_ = 4;
							resetStep();
						}
						else if (ntest_ == 0) {
							if (dispOpt_ != 0)
								std::cout << "Powell terminated: " << nt_ << " calls of calfun ineffective." << std::endl;
							run_ = 0;
							outflg_ = 3;
							resetStep();
						}
					}
				}
				else {
					ntest_ = nt_;
				}
			}
			if (run_) {
				if (maxc_ <= maxfun_ && maxi_ <= maxiter_) {
					if (dispOpt_ == 1) {
						std::cout << right << setw(numWidth) << setfill(separator) << maxi_;
						std::cout << right << setw(numWidth) << setfill(separator) << maxc_;
						std::cout << right << setw(numWidth) << setfill(separator) << x_;
						std::cout << right << setw(numWidth) << setfill(separator) << fval_;
						std::cout << right << setw(numWidth) << setfill(separator) << Fsq_;
						std::cout << right << setw(numWidth) << setfill(separator) << sqrt(dd_) << std::endl;
					}
					//std::cout << maxi_ << "\t\t\t" << maxc_ << "\t" << (*x_)(0, 0) << "\t" << (*x_)(1, 0) << "\t" << (*fval_)(0, 0) << "\t" << (*fval_)(1, 0) << "\t" << Fsq_ << "\t" << sqrt(dd_) << std::endl;
					switch (is_) {
					case 1:

						reviseDelta();
						break;
					case 2:

						reviseJac();
						break;
					case 3:
						evalJacobian1(&x_, &Jac_, _fun, &fval_, &dstep_);
						invJac_ = 1 / Jac_;
						if (debug_) {
							std::cout << "J = " << std::endl;
							std::cout << Jac_ << std::endl;
							std::cout << "J-1 = " << std::endl;
							std::cout << invJac_ << std::endl;
						}
						maxc_ += N_;
						mainstep();
						break;
					case 4:
						if (Fsq_ > Fmin_) {
							interchange();
						}
						else {

							is_ = 2;
							setSpecialStep();
						}
						break;
					case 5:

						Fmin_ = Fsq_;
						lastx_ = x_;
						lastfval_ = fval_;
						is_ = 3;
						break;
					}
				}
				else {
					if (maxc_ >= maxfun_) {
						if (dispOpt_ != 0)
							std::cout << "Powell terminated: maximum number of function evaluations reached." << std::endl;
						run_ = 0;
						outflg_ = 2;
					}
					else {
						if (dispOpt_ != 0)
							std::cout << "Powell terminated: maximum number of Iterations reached." << std::endl;
						run_ = 0;
						outflg_ = 1;
					}
					if (Fsq_ > Fmin_) {
						resetStep();
					}
				}
			}
		}
		else {
			if (dispOpt_ == 1) {
				std::cout << right << setw(numWidth) << setfill(separator) << maxi_;
				std::cout << right << setw(numWidth) << setfill(separator) << maxc_;
				std::cout << right << setw(numWidth) << setfill(separator) << x_;
				std::cout << right << setw(numWidth) << setfill(separator) << fval_;
				std::cout << right << setw(numWidth) << setfill(separator) << Fsq_;
				std::cout << right << setw(numWidth) << setfill(separator) << sqrt(dd_) << std::endl;
			}
			solvec_ = x_;
			if (isnan(Fsq_) == 1) {
				if (dispOpt_ != 0)
					std::cout << "Powell terminated: Function value is a complex number." << std::endl;
				run_ = 0;
				outflg_ = 6;
			}
			else {
				if (dispOpt_ != 0)
					std::cout << "Powell terminated: Solution found." << std::endl;
				run_ = 0;
				outflg_ = 0;
			}

		}
		if (debug_)
			std::cin.get();
	}
};

void powell_1D::initialize() {

	w_ = 1;
	W_ = 1;
	phi_ = 0;
	step_ = 0;

	// doubles
	tinc_ = 1;
	//acc_ = 1e-6;
	/*dmax_ = 10.0;*/
	dm_ = dmax_ * dmax_;
	dmm_ = 4 * dm_;
	/*dstep_ = (1e-6);*/
	dd_ = 0;
	dss_ = dstep_ * dstep_;
	dtest_ = 2 * double(1) - 0.5;
	Fsq_ = 0;
	mu_ = 0;
	sp_ = 0;
	ds_ = 0;
	dn_ = 0;
	gm_ = 0;
	anmult_ = 0;
	dmult_ = 0;
	Phi_ = 0;
	tinc_ = 1;
	ss_ = 0;
	pj_ = 0;
	alpha_ = 0;
	sval_ = 0;

	// bools/ints
	nt_ = 1 + 4;
	ntest_ = nt_;
	maxc_ = 0;
	//maxiter_ = 1000;
	maxi_ = 0;
	//maxfun_ = 10000;
	is_ = 5;
	run_ = 1;
	flag_ = 0;
	is_ = 5;
	kk_ = 0;
	outflg_ = -1;
	// outflg_ = 0 -> solution found
	// outflg_ = 1 -> maxiter_ reached
	// outflg_ = 2 -> maxfun_ reached
	// outflg_ = 3 -> nt_ steps ineffective
	// outflg_ = 4 -> failed to decrease Fsq_ using a new Jacobian
	// outflg_ = 5 -> nearby stationary point predicted
	// outflg_ = 6 -> function value is a complex number

	separator = ' ';
	nameWidth = 6;
	numWidth = 14;
};

void powell_1D::calfun(double* x, double* f, double* Fsq) {
	//(*_fun).Zero_Functionction(x, f);
	(*Fsq) = 0;
	(*Fsq) += (*f) * (*f);

	return;
};

void powell_1D::reviseDelta() {
	if (debug_)
		std::cout << "reviseDelta(), is_ = " << is_ << std::endl;
	dmult_ = 0.9 * Fmin_ + 0.1 * Phi_;
	if (Fsq_ > dmult_) {
		dd_ = max(0.25 * dd_, dss_);
		tinc_ = 1;
		if (Fsq_ >= Fmin_) {
			reviseJac();
		}
		else {
			interchange();
		}
	}
	else {
		dmult_ = dmult_ - Fsq_;
		sp_ = 0; ds_ = 0;  ss_ = 0;
		sp_ += fabs(fval_ * (fval_ - phi_));
		ds_ = (fval_ - phi_);
		ss_ += ds_ * ds_;

		pj_ = 1 + dmult_ / (sp_ + sqrt(sp_ * sp_ + dmult_ * ss_));
		pj_ = sqrt(pj_);
		mu_ = min(2.0, min(pj_, tinc_));
		tinc_ = pj_ / mu_;
		dd_ = min(mu_ * sqrt(dd_), dmax_);
		dd_ *= dd_;
		interchange();
	}
	if (debug_)
		std::cout << "new Delta = " << sqrt(dd_) << std::endl;

};

void powell_1D::reviseJac() {

	gamma_ = fval_ - lastfval_;
	step_ = x_ - lastx_;
	if (debug_) {
		std::cout << "reviseJac(), is_ = " << is_ << std::endl;
		std::cout << "gamma = " << std::endl;
		std::cout << gamma_ << std::endl;
	}
	double singtest;
	double singtest_;
	singtest = step_ * (invJac_ * gamma_);
	singtest_ = fabs(singtest);
	ds_ = pow(step_, 2);
	if (singtest_ >= 0.1 * ds_) {
		alpha_ = 1;
	}
	else {
		alpha_ = 0.8;
	};
	Jac_ = Jac_ + (alpha_ * (gamma_ - (Jac_ * step_)) * (step_ / ds_));

	invJac_ = invJac_ + (alpha_ * (((step_ - (invJac_ * gamma_)) * step_) * invJac_) / (alpha_ * singtest_ + (1 - alpha_) * ds_));
	if (debug_) {
		std::cout << "alpha = " << alpha_ << std::endl;
		std::cout << "J = " << std::endl;
		std::cout << Jac_ << std::endl;
		std::cout << "J-1 = " << std::endl;
		std::cout << invJac_ << std::endl;;
	}
	mainstep();

};

void powell_1D::resetStep() {
	if (debug_)
		std::cout << "resetStept(), is_ = " << is_ << std::endl;
	x_ = lastx_;
	fval_ = lastfval_;
	Fsq_ = Fmin_;
};

void powell_1D::interchange() {
	if (debug_)
		std::cout << "interchange(), is_ = " << is_ << std::endl;
	Fmin_ = Fsq_;
	tempx_ = x_;
	x_ = lastx_;
	lastx_ = tempx_;
	tempfval_ = fval_;
	fval_ = lastfval_;
	lastfval_ = tempfval_;
	phi_ = ((-1) * phi_);
	if (is_ <= 1) {
		reviseJac();
	}
	else {
		is_ = 2;
		setSpecialStep();
	}
};

void powell_1D::setSpecialStep() {
	if (debug_)
		std::cout << "setSpecialStep(), is_ = " << is_ << std::endl;

	step_ = dstep_ * W_;
	w_ = 1;
	W_ = step_ / 2 / ds_;
};

void powell_1D::mainstep() {
	if (debug_)
		std::cout << "mainstep(), is_ = " << is_ << std::endl;
	g_ = 0;
	g_ -= (Jac_ * lastfval_);

	v_ = (-1) * invJac_ * lastfval_;
	ds_ = pow(g_, 2);
	dn_ = pow(v_, 2);
	sp_ = 0;
	sp_ += g_ * v_;

	if (Fmin_ * Fmin_ - dmm_ * ds_ <= 0) {
		is_ = 2;
		if (dn_ > dd_) {
			flag_ = 1;
			while (flag_ == 1) {
				//std::cout << "flag_ = " << flag_ << std::endl;
				gtemp_ = Jac_ * g_;
				dmult_ = pow(gtemp_, 2);
				dmult_ = ds_ / dmult_; //mu
				ds_ = ds_ * dmult_ * dmult_;
				if (ds_ - dd_ < 0) {
					//std::cout << "setstep3, ds_ =" << ds_ << std::endl;
					flag_ = 0;
					sp_ = sp_ * dmult_;
					if ((sp_ - dd_) * (sp_ - dd_) + (dn_ - dd_) * (dd_ - ds_) < 0)
						anmult_ = 0;
					else
						anmult_ = (dd_ - ds_) / ((sp_ - ds_) + sqrt((sp_ - dd_) * (sp_ - dd_) + (dn_ - dd_) * (dd_ - ds_)));

					if (isinf(anmult_) == 1) {
						anmult_ = 0;
					}

					dmult_ = dmult_ * (1 - anmult_);

				}
				else {
					if (dd_ <= 0) {
						dd_ = max(dss_, min(dm_, ds_));
						ds_ = ds_ / (dmult_ * dmult_);
						//std::cout << "initial dd_ = " << sqrt(dd_) << std::endl;
					}
					else {
						flag_ = 0;
						//std::cout << "setstep2" << std::endl;
						anmult_ = 0;
						dmult_ = dmult_ * sqrt(dd_ / ds_);
					}
				}
				if (flag_ == 0) {
					step_ = dmult_ * g_ + (anmult_ * v_);
					if (debug_) {
						std::cout << "step = " << std::endl;
						std::cout << step_ << std::endl;
						std::cout << "steplength = " << step_ << std::endl;
					}
					dn_ = pow(step_, 2);
					ds_ = 0.25 * dn_;
					sp_ = 0;
					sp_ += step_ * W_;

					if (w_ - dtest_ <= 0) {
						updateOmega();
					}
					else {
						if (sp_ * sp_ - ds_ >= 0) {
							updateOmega();
						}
						else {
							setSpecialStep();
						}
					}
				}
			}
		}
		else {
			//std::cout << "setstep1" << std::endl;
			step_ = v_;
			if (debug_) {
				std::cout << "step = " << std::endl;
				std::cout << step_ << std::endl;
				std::cout << "steplength = " << step_ << std::endl;
			}
			dd_ = max(dn_, dss_);
			ds_ = 0.25 * dn_;
			tinc_ = 1;
			if (dn_ - dss_ >= 0) {
				updateOmega();
			}
			else {
				is_ = 4;
				setNewPoint();
			}
		}
	}
	else {
		if (is_ == 1 || is_ == 2) {
			ntest_ = 0;
			(x_) = (lastx_);
			is_ = 3;
		}
		else {
			if (dispOpt_ != 0)
				std::cout << "Powell terminated: nearby stationary point predicted." << std::endl;
			run_ = 0;
			outflg_ = 5;
			resetStep();
		}
	}
};

void powell_1D::updateOmega() {
	if (debug_)
		std::cout << "updateOmega(), is_ = " << is_ << std::endl;
	sp_ = 0;
	a_ = 0;
	a_ += step_ * W_;
	if (is_ == 1) {
		w_ = w_ + 1;
	}
	else {
		w_ = 1 + w_;
		sp_ += a_ * a_;
		if (sp_ - ds_ <= 0) //*ds_
			;//Wflag_ = 0;
		else {
			is_ = 1;
			kk_ = 0;
		}
	}

	w_ = 1;

	sval_ = a_ * a_;

	sig_ = 0;

	sval_ += a_ * a_;

	W_ = step_ / 2 / ds_;

	if (debug_) {
		std::cout << "W = " << std::endl;
		std::cout << W_ << std::endl;
	}
	setNewPoint();
};

void powell_1D::setNewPoint() {
	if (debug_)
		std::cout << "setNewPoint(), is_ = " << is_ << std::endl;
	Phi_ = 0;
	x_ = lastx_ + step_;
	phi_ = lastfval_ + (Jac_ * step_);
	Phi_ = pow(phi_, 2);
};

double powell_1D::solution() {
	double sol;
	if (debug_)
		std::cout << "Solution(), is_ = " << is_ << std::endl;
	sol = solvec_;
	return sol;
};

int powell_1D::outflg() {

	if (debug_)
		std::cout << "outflg(), is_ = " << is_ << std::endl;
	return outflg_;
};

void evalJacobian1(double* vec, double* Jacobian, Zero_Function* A, double* fval, double* dstep)
{
	double dx = *dstep;
	int sizeVec = sizeof(vec) / sizeof(vec[0]);

	std::vector<double> temp_vec;
	std::vector<double> temp_fval;

	temp_vec.assign((*vec), (*vec) + sizeVec);

	for (int i = 0; i < sizeVec; i++) {
		temp_vec[i] += dx;
	}

	(*A).zero_function(&temp_vec, &temp_fval);

	for (int i = 0; i < sizeVec; i++) {
		Jacobian[i] = (temp_fval[i] - fval[i]) / dx;
	}

	temp_vec.assign((*vec), (*vec) + sizeVec);

	return;
};
