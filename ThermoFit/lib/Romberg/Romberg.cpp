/*
Numerical integration with Romberg's method
*/
#include "Romberg.hpp"



double romberg_old(double (*func)(double), double *a, double *b) {
	const int N = 22;
	double h[N + 1], r[N + 1][N + 1];
	double coeff = 0;
	for (int i = 1; i <= N; ++i) {
		h[i] = ((*b) - (*a)) / pow(2, i - 1);
	}

	r[1][1] = h[1] / 2 * (func((*a)) + func((*b)));
	for (int i = 2; i < N + 1; ++i) {
		coeff = 0;
		for (int k = 1; k <= pow(2, i - 2); ++k) {
			coeff += func((*a) + (2 * k - 1) * h[i]);
		}
		r[i][1] = 0.5 * (r[i - 1][1] + h[i - 1] * coeff);
	}

	for (int i = 2; i < N + 1; ++i) {
		for (int j = 2; j <= i; ++j) {
			r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (pow(4, j - 1) - 1);
		}
	}
	return r[N][N];
}

// constructors
ROM::ROM(para_romberg* paraROM, func_ROM *fun_ROM) :
	paraROM_(paraROM),
	fun_(fun_ROM),
	max_steps_((*paraROM).max_steps),
	max_steps_ini_((*paraROM).max_steps_ini),
	min_steps_((*paraROM).min_steps),
	acc_((*paraROM).acc),
	debug_((*paraROM).debug) {
	r_ = new matrix (max_steps_ini_, max_steps_ini_);	// matrix that contains values of each row
	h_ = new double[max_steps_ini_];				// not initialized with = 0
	for (int i = 0; i < max_steps_ini_; i++) {		// initialize with = 0
		h_[i] = 0.0;
	}
};

ROM::ROM(func_ROM *fun_ROM) : 
	fun_(fun_ROM) {
	para_romberg paraROM;
	paraROM_ = &paraROM;
	max_steps_ = paraROM.max_steps;
	max_steps_ini_ = max_steps_;
	min_steps_ = paraROM.min_steps;
	acc_ = paraROM.acc,
	debug_ = paraROM.debug;
	r_ = new matrix(max_steps_ini_, max_steps_ini_);	// matrix that contains values of each row
	h_ = new double[max_steps_ini_];				// not initialized with = 0
	for (int i = 0; i < max_steps_ini_; i++) {		// initialize with = 0
		h_[i] = 0.0;
	}
};
// destructor
ROM::~ROM() {			
	delete[] h_;
	if (debug_)std::cout << "\ndestructor ROM called\n" << std::endl;
	//delete[] r_;
};

// Update parameters from paraROM struct
void ROM::set_paraROM(para_romberg* paraROM) {
	for (int i = 0; i < max_steps_ini_; i++) {		// initialize with = 0
		h_[i] = 0.0;
	}
	paraROM_ = paraROM;
	max_steps_ = (*paraROM).max_steps;
	min_steps_ = (*paraROM).min_steps;
	acc_ = (*paraROM).acc,
	debug_ = (*paraROM).debug;
};

/* 
	eventuell auf Vektoren umschreiben und Größe dynamisch verändern, statt Matrix mit max_step*max_step zu verwenden
*/

void ROM::calc_fun(double* x1, double *fval) {
	//(*fval) = (*fun_).func_rom((*x1));
	(*fun_).func_rom(x1, fval);
	return;
}

double ROM::calc_fun_double(double* x1) {
	double fval;
	(*fun_).func_rom(x1, &fval);
	return fval;
}

double ROM::romberg_old(double *a, double *b) {
	const int N = 16;
	double h[N + 1], r[N + 1][N + 1];
	double coeff = 0;
	for (int i = 1; i <= N; ++i) {
		h[i] = ((*b) - (*a)) / pow(2, i - 1);
	}

	r[1][1] = h[1] / 2 * (calc_fun_double(a) + calc_fun_double(b));
	for (int i = 2; i < N + 1; ++i) {
		coeff = 0;
		for (int k = 1; k <= pow(2, i - 2); ++k) {
			double bnd_ = (*a) + (2 * k - 1) * h[i];
			coeff += calc_fun_double(&bnd_);
		}
		r[i][1] = 0.5 * (r[i - 1][1] + h[i - 1] * coeff);
	}

	for (int i = 2; i < N + 1; ++i) {
		for (int j = 2; j <= i; ++j) {
			r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (pow(4, j - 1) - 1);
		}
	}
	return r[N][N];
}

double ROM::romberg(double *a, double *b) {
	double d_coeff = 0;							// coefficient  for calculation R(n, 0)
	double diff_ = 1E10;
	if (debug_ == 1) std::cout << "\nStarting Romberg....\n" << "requested accuracy: " << acc_ << std::endl;
	if (debug_ == 1) std::cout << "pow(2, (max_steps_ - 1)) = " << pow(2, (max_steps_ - 1)) << std::endl;

	/// Romberg Sequence
	for (int i = 0; i < max_steps_; i++) {		// step size with respect to n, h[0] = b-a
		h_[i] = ((*b) - (*a)) / pow(2, (i + 1));
	}
	// first trapezoid R(0, 0)
	calc_fun(a, &fval_a);
	calc_fun(b, &fval_b);
	(*r_)(0, 0) = h_[0] * (fval_a + fval_b);
	//(*r_)(0, 0) = h_[0] * ((*fun_).func_rom((*a)) + (*fun_).func_rom(*b));
	if (debug_ == 1) std::cout << "fval(" << (*a) << ") = " << fval_a << std::endl;
	if (debug_ == 1) std::cout << "fval_1(" << (*b) << ") = " << fval_b << std::endl;
	if (debug_ == 1) std::cout << "(*r_)(0,0) = " << (*r_)(0, 0) << std::endl;
	//if (debug_ == 1) std::cin.get();

	// compute R(n, 0)
	//parallel_for(1, max_steps_, [&](int n) {
	for (int n = 1; n < max_steps_; n++) {
		d_coeff = 0;							// set coefficient to 0 at each iteration cycle
		for (int k = 1; k <= pow(2, (n - 1)); k++) {
			//int upper = pow(2, (n - 1)) + 1;
			a1 = (*a) + (2 * k - 1) * h_[n - 1];
			calc_fun(&a1, &fval_a);
			//d_coeff += (*fun_).func_rom((*a) + (2 * k - 1) * h_[n - 1]);
			d_coeff += fval_a;
		}
		(*r_)(n, 0) = 0.5*(*r_)(n - 1, 0) + h_[n - 1] * d_coeff;
	}
	//);
	// compute R(n, m)
	for (int n = 1; n < max_steps_; n++) {
		for (int m = 1; m <= n; m++) {
			(*r_)(n, m) = (pow(4, m)*(*r_)(n, m - 1) - (*r_)(n - 1, m - 1)) / (pow(4, m) - 1);
		}
		diff_ = fabs((*r_)(n, n) - ((*r_)(n - 1, n - 1))) / (*r_)(n, n);
		if (fabs(diff_)< acc_ && n > min_steps_) {
			if (debug_ == 1) std::cout << "achieved requested accuracy, diff = " << diff_ << std::endl << "iterations: " << n << std::endl;
			(*paraROM_).n_steps = n;		// store number of necessary iterations
			return (*r_)(n, n);
		}
	}
	if (debug_ == 1) {
		(*r_).disp();
		std::cout << "reached maximum number of iteration steps (" << max_steps_ << ")\n" << std::endl;
	}
	(*paraROM_).n_steps = max_steps_;		// store number of necessary iterations
	return (*r_)(max_steps_ - 1, max_steps_ - 1);
}

double ROM::romberg2(double *a, double *b) {
	double d_coeff = 0;							// coefficient  for calculation R(n, 0)
	double diff_ = 1E10;
	if (debug_ == 1) std::cout << "\nStarting Romberg....\n" << "requested accuracy: " << acc_ << std::endl;
	
	// Bulirsch Sequence
	h_[0] = ((*b) - (*a));
	for (int i = 1; i < max_steps_/2+1; i++) {		// step size with respect to n, h[0] = b-a
		if (2*i-1 < max_steps_)h_[2*i-1] = ((*b) - (*a)) / pow(2, (i));
		if (2*i < max_steps_) h_[2*i] = ((*b) - (*a)) / (3*pow(2, i-1));
	}
	if (debug_ == 1) std::cout << "h_[0] / h_[max_steps_ - 1] - 1 = " << h_[0] / h_[max_steps_ - 1] - 1 << std::endl;

	//// Romberg Sequence
	//for (int i = 0; i < max_steps_; i++) {
	//	h_[i] = ((*b) - (*a)) / pow(2, (i));
	//}
	//if (debug_ == 1) std::cout << "pow(2, (max_steps_ - 1)) = " << pow(2, (max_steps_ - 1)) << std::endl;

	// Print h
	if (debug_ == 2) {
		for (int i = 0; i < max_steps_; i++) {

		}
	}
	// first trapezoid R(0, 0)
	calc_fun(a, &fval_a);
	calc_fun(b, &fval_b);
	(*r_)(0, 0) = 0.5*h_[0] * (fval_a + fval_b);
	//(*r_)(0, 0) = h_[0] * ((*fun_).func_rom((*a)) + (*fun_).func_rom(*b));
	if (debug_ == 1) std::cout << "fval(" << (*a) << ") = " << fval_a << std::endl;
	if (debug_ == 1) std::cout << "fval_1(" << (*b) << ") = " << fval_b << std::endl;
	if (debug_ == 1) std::cout << "(*r_)(0,0) = " << (*r_)(0, 0) << std::endl;
	//if (debug_ == 1) std::cin.get();

	// compute R(n, 0)
	for (int n = 2; n <= max_steps_; n++) {
		// Trapezoidal sum, compute R(n, 0)
		d_coeff = 0;
		for (int j = 1; j < h_[0]/ h_[n-1]-1; j++) {
			a1 = (*a) + j*h_[n - 1];
			calc_fun(&a1, &fval_j);
			d_coeff += fval_j;
		}
		(*r_)(n-1, 0) = 0.5*h_[n-1]*(fval_a+ fval_b+2.0*d_coeff);
	}
	//// compute R(n, m)
	//for (int k = 1; k <= max_steps_; k++) {
	//	// Richardson extrapolation, compute R(n, k)
	//	(*r_)(n - 1, k - 1) = (*r_)(n - 1, k - 2) + ((*r_)(n - 1, k - 2) - (*r_)(n - 2, k - 1)) / (pow(4, k-1)-1);
	//}

	//// compute R(n, m)
	for (int k = 1; k < max_steps_; k++) {
		for (int n = 1; n < max_steps_; n++) {
		//for (int n = 1; n < k; n++) {
			//(*r_)(n, k) = (*r_)(n + 1, k - 1) + ((*r_)(n + 1, k - 1) - (*r_)(n, k - 1)) / (pow(h_[n] / h_[n + 1], 2 * k) - 1);	// Access violation for n = max_steps_, but wrong result if for-loop with n < max_steps_-1, 
			(*r_)(n, k) = (*r_)(n, k - 1) + ((*r_)(n , k - 1) - (*r_)(n - 1, k - 1)) / (pow(h_[n-k+1-1] / h_[n],2) - 1);			// For Bullirsch-sequence; Access violation for k>n => no error, but needs ~ 3x more function evals than with Romberg sequence to achieve stable result (tested with test function 1/(1+25*x^2) and interfacial concentration profile)
			//(*r_)(n, k) = (*r_)(n, k - 1) + ((*r_)(n, k - 1) - (*r_)(n - 1, k - 1)) / (pow(4,k) - 1);								// For Romberg-sequence
			//(*r_)(k, n) = (*r_)(k, n - 1) + ((*r_)(k, n - 1) - (*r_)(k - 1, n - 1)) / ( pow(h_[k] / h_[k-n+1],2) - 1 );
		}
		diff_ = fabs((*r_)(k, k) - ((*r_)(k - 1, k - 1))) / (*r_)(k, k);
		//if (debug_ == 1) std::cout << "diff = " << diff_ << std::endl;
		if (diff_< acc_ && k > min_steps_) {
			if (debug_ == 1) std::cout << "achieved requested accuracy, diff = " << diff_ << std::endl << "iterations: " << k << std::endl;
			(*paraROM_).n_steps = k;		// store number of necessary iterations
			//return (*r_)(k, k);
		}
	}
	if (debug_ == 1) {
		(*r_).disp();
		std::cout << "reached maximum number of iteration steps (" << max_steps_ << ")\n" << std::endl;
	}
	//(*r_).disp();
	(*paraROM_).n_steps = max_steps_;		// store number of necessary iterations
	return (*r_)(max_steps_ - 1, max_steps_ - 1);
}

//	std::vector<std::vector<double>> (*r_);
//	double d_coeff = 0;							// coefficient  for calculation R(n, 0)
//	double diff_ = 1E10;
//	if (debug_ == 1) std::cout << "Starting Romberg....\n" << "requested accuracy: " << acc_ << "\n" << std::endl;
//
//	for (int i = 0; i < max_steps_; i++) {		// step size with respect to n, h[0] = b-a
//		h_[i] = ((*b) - (*a)) / pow(2, (i + 1));
//	}
//	// first trapezoid R(0, 0)
//	(*r_)[0][0] = h_[0] * (func((*a)) + func(*b));
//
//	// compute R(n, 0)
//	for (int n = 1; n < max_steps_; n++) {
//		d_coeff = 0;							// set coefficient to 0 at each iteration cycle
//		for (int k = 1; k <= pow(2, (n - 1)); k++) {
//			d_coeff += func((*a) + (2 * k - 1) * h_[n - 1]);
//		}
//		(*r_)[n][0] = 0.5*(*r_)[n - 1][0] + h_[n - 1] * d_coeff;
//	}
//
//	// compute R(n, m)
//	for (int n = 1; n < max_steps_; n++) {
//		for (int m = 1; m <= n; m++) {
//			(*r_)[n][m] = (pow(4, m)*(*r_)[n][m - 1] - (*r_)[n - 1][m - 1]) / (pow(4, m) - 1);
//		}
//		diff_ = fabs((*r_)[n][n] - ((*r_)[n - 1] [n - 1]));
//		if (diff_< acc_ && n > min_steps_) {
//			if (debug_ == 1) std::cout << "achieved requested accuracy, diff = " << diff_ << std::endl << "iterations: " << n << std::endl;
//			return (*r_)[n][n];
//		}
//	}
//	if (debug_ == 1) {
//		//(*r_).disp();
//		std::cout << "reached maximum number of iteration steps (" << max_steps_ << ")\n" << std::endl;
//	}
//	return (*r_)[max_steps_ - 1][max_steps_ - 1];
//}

//double ROM::romberg2(double *a, double *b) {
//	double d_coeff = 0;							// coefficient  for calculation R(n, 0)
//	double diff_ = 1E10;
//	if (debug_ == 1) std::cout << "\nStarting Romberg....\n" << "requested accuracy: " << acc_ << std::endl;
//
//	for (int i = 0; i < max_steps_; i++) {		// step size with respect to n, h[0] = b-a
//		h_[i] = ((*b) - (*a)) / pow(2, (i + 1));
//	}
//	// first trapezoid R(0, 0)
//
//	(*r_)(0, 0) = h_[0] * ((*fun_).func_rom((*a)) + (*fun_).func_rom(*b));
//	//if (debug_ == 1) std::cout << "(*fun_).func_rom(" << (*a) << ") = " << (*fun_).func_rom((*a)) << std::endl;
//	//if (debug_ == 1) std::cout << "(*fun_).func_rom(" << (*b) << ") = " << (*fun_).func_rom((*b)) << std::endl;
//	//if (debug_ == 1) std::cout << "(*r_)(0,0) = " << (*r_)(0, 0) << std::endl;
//	//if (debug_ == 1) std::cin.get();
//
//	// compute R(n, 0)
//	for (int n = 1; n < max_steps_; n++) {
//		d_coeff = 0;							// set coefficient to 0 at each iteration cycle
//		for (int k = 1; k <= pow(2, (n - 1)); k++) {
//			d_coeff += (*fun_).func_rom((*a) + (2 * k - 1) * h_[n - 1]);
//		}
//		(*r_)(n, 0) = 0.5*(*r_)(n - 1, 0) + h_[n - 1] * d_coeff;
//	}
//
//	// compute R(n, m)
//	for (int n = 1; n < max_steps_; n++) {
//		for (int m = 1; m <= n; m++) {
//			(*r_)(n, m) = (pow(4, m)*(*r_)(n, m - 1) - (*r_)(n - 1, m - 1)) / (pow(4, m) - 1);
//		}
//		diff_ = fabs((*r_)(n, n) - ((*r_)(n - 1, n - 1)));
//		if (diff_< acc_ && n > min_steps_) {
//			if (debug_ == 1) std::cout << "achieved requested accuracy, diff = " << diff_ << std::endl << "iterations: " << n << std::endl;
//			return (*r_)(n, n);
//		}
//	}
//	if (debug_ == 1) {
//		(*r_).disp();
//		std::cout << "reached maximum number of iteration steps (" << max_steps_ << ")\n" << std::endl;
//	}
//	return (*r_)(max_steps_ - 1, max_steps_ - 1);
//}