#pragma once
/**
CS-11, matrix.hpp
shared code for various calling functions
Purpose:
matrix class and operations on matrices

@author Patrick Zimmermann, TU Graz
@year 2017
*/
#include<string>
#include <iostream>
#include<vector>
#include<cmath>
#include<ctime>
#include<omp.h>
#include<algorithm>
#include <cstddef>
#include<functional>
#include<fstream>
#include<sstream>
#include <iomanip>
#include <limits>
#include "../core/global.hpp"
using namespace std;

#ifndef matrix_h
#define matrix_h
#define index(i,j)((i)*dim_cols_ + (j))

// Klassendefinition
class matrix {
public:
	matrix(int dim_rows = 0, int dim_cols = 0); // Konstruktor
	matrix(int dim_rows, int dim_cols, double value); // Überladener Konstruktor
	matrix(const matrix& rhs); // Kopierkonstruktor
							   // Destruktor
	~matrix() {
		//std::cout << "dtor call: " << vals_ << std::endl;
		delete[] vals_; // Speicher freigeben
	};

	// Nötig für Kopierkonstruktor << nötig für Rechenoperationen und matrix-Klasse als Rückgabewert von Funktionen
	friend void swap(matrix& first, matrix& second);
	// nötig für Rechenoperationen und matrix-Klasse als Rückgabewert von Funktionen
	// A = B*C Speicherinhalte werden umgedreht!
	matrix& operator=(matrix other) // (1)
	{

		//if ((other).arrSize_ >= 100) {
		swap(*this, other); // (2)
	//}
	//else {
	//	(*this).dim_rows_ = other.dim_rows_;
	//	(*this).dim_cols_ = other.dim_cols_;
	//	(*this).arrSize_ = other.arrSize_;
	//	memcpy((*this).vals_, other.vals_, other.arrSize_ * sizeof(double));
	//}
		return *this;
	}

	matrix& operator=(matrix* other) // (1)
	{

		if (this != other) {
			(*this).dim_rows_ = (*other).dim_rows_;
			(*this).dim_cols_ = (*other).dim_cols_;

			if ((*this).arrSize_ != (*other).arrSize_) {
				delete[](*this).vals_;
				(*this).vals_ = new double[(*other).arrSize_];
			}
			memcpy((*this).vals_, (*other).vals_, (*other).arrSize_ * sizeof(double));
		}
		return *this;
	}

	void mat_input(int dim_rows, int dim_cols, double* inputarr);

	int rows()const { return dim_rows_; };
	int cols()const { return dim_cols_; };

	// Liefern Index kleinstem/größtem Wert in Zeile/Spalte
	int row_min(int row);
	int col_min(int col);
	int row_max(int row);
	int col_max(int col);

	// Liefern Matrizen mit allen Spaltenmaxima/-minima
	matrix max();
	matrix min();

	// Speicheradressen
	double* begin()const { return vals_; }
	double* end()const { return vals_ + dim_rows_ * dim_cols_; }
	double euclidMetric();

	void row(int row_qry, matrix* row_vec);
	void col(int col_qry, matrix* col_vec);
	void mat_transpose(matrix* transposed);
	void mat_inverse(matrix* Ainv);
	void mat_inverse(matrix* Ainv, matrix* fieldA);
	void swap_row(int a, int b);
	void swap_col(int a, int b);
	void writetofile(std::string file_name);
	void readfromfile(matrix* A, std::string file_name);
	void write_to_file(ofstream* writefile);

	friend matrix mat_add(matrix* A, matrix* B);
	friend matrix mat_sub(matrix* A, matrix* B);
	friend matrix mat_mult(matrix* A, matrix* B);
	friend matrix mat_eq(matrix* A, matrix* B);

	void disp();
	int* size();
	bool chkzero();

	double& operator()(int i, int j) { return vals_[index(i, j)]; };
	double* operator[](int i) { return vals_ + i * dim_rows_; };

	matrix operator+(matrix A) { return mat_add(this, &A); };
	matrix operator+(matrix* A) { return mat_add(this, A); };
	matrix operator+=(matrix A) { return (*this + mat_add(this, &A)); };
	matrix operator-(matrix A) { return mat_sub(this, &A); };
	matrix operator-(matrix* A) { return mat_sub(this, A); };
	matrix operator-=(matrix A) { return (*this - mat_sub(this, &A)); };
	matrix operator*(matrix A) { return mat_mult(this, &A); };
	matrix operator*(matrix* A) { return mat_mult(this, A); };
	matrix operator*=(matrix A) { return (*this + mat_mult(this, &A)); };
	matrix operator==(matrix A) { return (mat_eq(this, &A)); };

private:
	int dim_rows_;
	int dim_cols_;
	std::size_t arrSize_;
	int row_qry_;
	int col_qry_;
	double* vals_;
};

// Operatorendeklarationen
matrix operator+(matrix A, double b);
matrix operator+(double b, matrix A);
matrix operator-(matrix A, double b);
matrix operator-(double b, matrix A);
matrix operator*(matrix A, double b);
matrix operator*(double b, matrix A);
matrix operator/(matrix A, double b);
matrix operator/(double b, matrix A);
matrix mat_add1(matrix* A, double b);
matrix mat_add2(matrix* A, double b);
matrix mat_sub(double b, matrix* A);
matrix mat_add1(double b, matrix* A);
matrix mat_mult2(matrix* A, double b);
matrix div_mat(double b, matrix* A);
matrix exp(matrix* A);
matrix log(matrix* A);
matrix sqrt(matrix* A);
matrix pow(matrix* A, double d);
matrix elementprod(matrix* A, matrix* B);
void elementprod(matrix* A, matrix* B, matrix* res);
matrix elementdiv(matrix* A, matrix* B);
matrix mat_mult2(double b, matrix* A);
matrix el_inv(matrix* A);
matrix mat_add(matrix* A, matrix* B);
matrix mat_eq(matrix* A, matrix* B);
matrix mat_sub(matrix* A, matrix* B);
matrix mat_mult(matrix* A, matrix* B);
void mat_mult(matrix* A, matrix* B, matrix* Result);
void evalJacobian(double vec[], matrix* Jacobian, int n, int m, void(*fun)(double*, double*), double fval[]);
void LUPivot(matrix* A, matrix* L, matrix* U, int* piv);
void GCholesky(matrix* A, matrix* G, bool* errflag);
void forback(matrix* L, matrix* U, matrix* vec, matrix* x);
double determinant(matrix* A);
matrix mldivide(matrix* A, matrix* b);


#endif
