#include "matrix.hpp"

// Konstruktor: Allokation von benötigtem Speicherplatz
matrix::matrix(int dim_rows, int dim_cols) :
	dim_rows_(dim_rows),
	dim_cols_(dim_cols),
	arrSize_(dim_rows* dim_cols) {
	if (dim_rows_ >= 1 && dim_cols_ >= 1) {
		vals_ = new double[arrSize_];	// NICHT MIT 0 INITIALISIERT!!!
		for (auto k = begin(); k != end(); k++)
			*k = 0.0;						// Jetzt schon!!!
	}
	//std::cout << "ctor call: " << vals_ << std::endl;
};

// Überladener Konstruktor, Besetzung mit gegebenen double:
matrix::matrix(int dim_rows, int dim_cols, double value) :
	dim_rows_(dim_rows),
	dim_cols_(dim_cols),
	arrSize_(dim_rows* dim_cols) {
	if (dim_rows_ >= 1 && dim_cols_ >= 1) {
		vals_ = arrSize_ ? new double[arrSize_] : nullptr;	// NICHT MIT 0 INITIALISIERT!!!
		for (auto k = begin(); k != end(); k++)
			*k = value;						// Jetzt schon!!!
	}
	//std::cout << "ctor call: " << vals_ << std::endl;
};

// Kopierkonstruktor
matrix::matrix(const matrix& rhs)
	:dim_rows_(rhs.dim_rows_),
	dim_cols_(rhs.dim_cols_),
	arrSize_(rhs.arrSize_),
	vals_(arrSize_ ? new double[arrSize_] : nullptr) {
	memcpy(vals_, rhs.vals_, arrSize_ * sizeof(double)); // std::copy will hier nicht laufen! 
	//
	//std::cout << "cctor call: " << vals_ << " = " << rhs.vals_ << std::endl;
	//std::cin.get();
	//
};

// Nötig für Kopierkonstruktor << A = B*C: Speicherinhalte werden umgedreht!
void swap(matrix& first, matrix& second) // nothrow
{
	// enable ADL (not necessary in our case, but good practice)
	using std::swap;
	swap(first.dim_rows_, second.dim_rows_);
	swap(first.dim_cols_, second.dim_cols_);
	swap(first.arrSize_, second.arrSize_);
	swap(first.vals_, second.vals_);
};

double matrix::euclidMetric() {
	double res = 0;
	if (dim_cols_ == 1) {
		for (int i = 0; i < dim_rows_; i++)
			res += vals_[index(i, 0)] * vals_[index(i, 0)];
	}
	else {
		std::cout << "Wrong dimensions to calculate euclidean metric." << std::endl;
		res = -1;
	}
	return res;
};

// Zeilenminimum
int matrix::row_min(int a) {
	int k = 0;
	double temp_min;
	temp_min = vals_[index(a, k)];
	for (int i = 1; i < dim_cols_; i++) {
		if (vals_[index(a, i)] < temp_min) {
			temp_min = vals_[index(a, i)];
			k = i;
		}
	}
	return k;
};

// Zeilenmaximum
int matrix::row_max(int a) {
	int k = 0;
	double temp_max;
	temp_max = vals_[index(a, k)];
	for (int i = 1; i < dim_cols_; i++) {
		if (vals_[index(a, i)] > temp_max) {
			temp_max = vals_[index(a, i)];
			k = i;
		}
	}
	return k;
};

// Spaltenminimum
int matrix::col_min(int a) {
	int k = 0;
	double temp_min;
	temp_min = vals_[index(k, a)];
	for (int i = 1; i < dim_rows_; i++) {
		if (vals_[index(i, a)] < temp_min) {
			temp_min = vals_[index(i, a)];
			k = i;
		}
	}
	return k;
};

// Spaltenmaximum
int matrix::col_max(int a) {
	int k = 0;
	double temp_max;
	temp_max = vals_[index(k, a)];
	for (int i = 1; i < dim_rows_; i++) {
		if (vals_[index(i, a)] > temp_max) {
			temp_max = vals_[index(i, a)];
			k = i;
		}
	}
	return k;
};

// Vektor der Spaltenminima
matrix matrix::min() {
	matrix min_vec(1, dim_cols_);
	for (int i = 0; i < dim_cols_; i++) {
		min_vec(0, i) = vals_[index(col_min(i), i)];
	}
	return min_vec;
};

// Vektor der Spaltenmaxima
matrix  matrix::max() {
	matrix max_vec(1, dim_cols_);
	for (int i = 0; i < dim_rows_; i++) {
		max_vec(0, i) = vals_[index(col_max(i), i)];
	}
	return max_vec;
};


// Zeilen in einer Matrix tauschen
void matrix::swap_row(int a, int b) {
	double temp_val;
	if ((dim_rows_ > 1) && (a < dim_rows_) && (b < dim_rows_)) {
		for (int i = 0; i < dim_cols_; i++) {
			temp_val = vals_[index(a, i)];
			vals_[index(a, i)] = vals_[index(b, i)];
			vals_[index(b, i)] = temp_val;
		}
	}
	else {
		std::cout << "Fehler: Zeilentausch nicht m\x94 \bglich." << std::endl;
	}
};

// Spalten in einer Matrix tauschen
void matrix::swap_col(int a, int b) {
	double temp_val;
	if ((dim_cols_ > 1) && (a < dim_cols_) && (b < dim_cols_)) {
		for (int i = 0; i < dim_rows_; i++) {
			temp_val = vals_[index(i, a)];
			vals_[index(i, a)] = vals_[index(i, b)];
			vals_[index(i, b)] = temp_val;
		}
	}
	else {
		std::cout << "Fehler: Spaltentausch nicht m\x94 \bglich." << std::endl;
	}

};

// Zeilenvektor angegebener Spalte ausgeben 
void matrix::row(int row_qry, matrix* row_vec) {
	for (int i = 0; i < dim_cols_; i++) {
		(*row_vec)(0, i) = vals_[index(row_qry, i)];
	}
	return;
};

// Spaltenvektor angegebener Spalte ausgeben 
void matrix::col(int col_qry, matrix* col_vec) {
	for (int i = 0; i < dim_rows_; i++) {
		(*col_vec)(i, 0) = vals_[index(i, col_qry)];
	}
	return;
};

// Matrix ausgeben
void matrix::disp() {

	int i = -1;
	int j = -1;
	for (auto k = begin(); k != end(); k++) {
		j++;
		if (!(j % dim_cols_)) {
			i++;
			j = 0;
		}
		if ((j == 0) && (j == dim_cols_ - 1))
			std::cout << "(" << "\t" << *k << "\t" << ")" << std::endl;
		else if (j == 0)
			std::cout << "(" << "\t" << *k << "\t";
		else if (j == dim_cols_ - 1)
			std::cout << *k << "\t" << ")" << std::endl;
		else
			std::cout << *k << "\t";
	}
	/*
	std::cout << "k: " << k << ", *k: " << *k << ", int(*k): " << int(*k) << std::endl;
	std::cout << "Value: " << j % dim_cols_ << std::endl;
	std::cout << "if-Abfrage: "<< !(j % dim_cols_) << std::endl;
	std::cout << "(" << i << ", " << j << ")" << std::endl;
	for (int i = 0; i < dim_rows_; i++) {
	for (int j = 0; j < dim_cols_; j++) {
	if ((j == 0) && (j == dim_cols_ - 1))
	std::cout << "(" << "\t" << vals_[index(i, j)] << "\t" << ")" << std::endl;
	else if (j == 0)
	std::cout << "(" << "\t" << vals_[index(i, j)] << "\t";
	else if (j == dim_cols_ - 1)
	std::cout << vals_[index(i, j)] << "\t" << ")" << std::endl;
	else
	std::cout << vals_[index(i, j)] << "\t";
	}
	}
	*/
};

// Gib Zeiger auf Speicherplatz mit Größenangaben zurück
int* matrix::size() {
	int* size[2];
	size[0] = &dim_rows_;
	size[1] = &dim_cols_;
	return *size;
};

// Sind alle Elemente 0?
bool matrix::chkzero() {
	bool chkval = true;
	int i = 0;
	while ((chkval == true) && (i < dim_rows_ * dim_cols_)) {
		if (vals_[i] != 0)
			chkval = false;
		i++;
	}
	return chkval;
};

void matrix::mat_input(int dim_rows, int dim_cols, double* inputarr) {
	for (int i = 0; i < dim_rows; i++) {
		for (int j = 0; j < dim_cols; j++) {
			vals_[index(i, j)] = inputarr[index(i, j)];
		}
	}
};

void matrix::mat_transpose(matrix* Transposed) {
	for (int i = 0; i < dim_rows_; i++)
		for (int j = 0; j < dim_cols_; j++)
			(*Transposed)(j, i) = vals_[index(i, j)];
	return;
};

// Berechnung der inversen
void matrix::mat_inverse(matrix* Ainv) {
	int m = dim_rows_;
	int n = dim_cols_;
	int is;
	double tempval;
	double det = 1.0;

	matrix fieldA(n, 2 * n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			fieldA(i, j) = vals_[index(i, j)];
		for (int j = n; j < 2 * n; ++j)
			fieldA(i, j) = (i == j - n) ? 1.0 : 0.0;
	}

	for (int i = 0; i < n - 1; i++) {
		is = i;
		bool swap = 0;
		for (int j = is + 1; j < n; j++) {
			if (fabs(vals_[index(i, j)]) > fabs(vals_[index(i, j)])) {
				is = j;
				swap = 1;
			}
		}
		if (swap) {
			//(*A).swap_row(i, is);
			fieldA.swap_row(i, is);

			//tmp = piv[i];
			//piv[i] = piv[is];
			//piv[is] = tmp;
		}
		for (int j = i + 1; j < n; j++) {
			tempval = fieldA(j, i) / fieldA(i, i);
			for (int k = i; k < 2 * n; k++) {
				fieldA(j, k) -= tempval * fieldA(i, k);
			}
		}
	}
	for (int i = 0; i < n; i++)
		det += fieldA(i, i);

	if (det != 0.0) {
		for (int i = n - 1; i > 0; i--) {
			for (int j = i - 1; j >= 0; j--) {
				tempval = fieldA(j, i) / fieldA(i, i);
				for (int k = i; k < 2 * n; k++) {
					fieldA(j, k) -= tempval * fieldA(i, k);
				}
			}
		}
		for (int i = 0; i < n; i++) {
			tempval = fieldA(i, i);
			for (int j = 0; j < n; j++) {
				(*Ainv)(i, j) = fieldA(i, j + n) / tempval;
			}
		}
	}
	else {
		std::cout << "A nicht invertierbar";
	};
};

void matrix::mat_inverse(matrix* Ainv, matrix* fieldA) {
	int m = dim_rows_;
	int n = dim_cols_;
	int is;
	double tempval;
	double det = 1.0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			(*fieldA)(i, j) = vals_[index(i, j)];
		for (int j = n; j < 2 * n; ++j)
			(*fieldA)(i, j) = (i == j - n) ? 1.0 : 0.0;
	}

	for (int i = 0; i < n - 1; i++) {
		is = i;
		bool swap = 0;
		for (int j = is + 1; j < n; j++) {
			if (fabs(vals_[index(i, j)]) > fabs(vals_[index(i, j)])) {
				is = j;
				swap = 1;
			}
		}
		if (swap) {
			//(*A).swap_row(i, is);
			(*fieldA).swap_row(i, is);

			//tmp = piv[i];
			//piv[i] = piv[is];
			//piv[is] = tmp;
		}
		for (int j = i + 1; j < n; j++) {
			tempval = (*fieldA)(j, i) / (*fieldA)(i, i);
			for (int k = i; k < 2 * n; k++) {
				(*fieldA)(j, k) -= tempval * (*fieldA)(i, k);
			}
		}
	}
	for (int i = 0; i < n; i++)
		det += (*fieldA)(i, i);

	if (det != 0.0) {
		for (int i = n - 1; i > 0; i--) {
			for (int j = i - 1; j >= 0; j--) {
				tempval = (*fieldA)(j, i) / (*fieldA)(i, i);
				for (int k = i; k < 2 * n; k++) {
					(*fieldA)(j, k) -= tempval * (*fieldA)(i, k);
				}
			}
		}
		for (int i = 0; i < n; i++) {
			tempval = (*fieldA)(i, i);
			for (int j = 0; j < n; j++) {
				(*Ainv)(i, j) = (*fieldA)(i, j + n) / tempval;
			}
		}
	}
	else {
		std::cout << "A nicht invertierbar";
	};
};

void matrix::writetofile(std::string file_name) {
	ofstream file;
	file.open(file_name);
	int i = -1;
	int j = -1;
	for (auto k = begin(); k != end(); k++) {
		j++;
		if (!(j % dim_cols_)) {
			i++;
			j = 0;
		}
		if ((j == 0) && (j == dim_cols_ - 1))
			file << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << *k << std::endl;
		else if (j == 0)
			file << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << *k << "\t";
		else if (j == dim_cols_ - 1)
			file << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << *k << "\t" << std::endl;
		else
			file << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << *k << "\t";
	}
	file.close();
};

void matrix::readfromfile(matrix* A, std::string file_name) {

	// Open file
	ifstream file;
	file.open(file_name);
	if (!(file.is_open())) {
		std::cout << "Unable to open matrix-file!" << std::endl;
		//std::cin.ignore();
		return;
	}

	// Initialize ints and strings
	int n_lines = 0;
	int n_cols = 0;
	char token = '\t';
	std::string line, result;

	// Parse file to determine n_cols and n_rows
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		n_lines++;
		if (n_lines == 1) {
			while (std::getline(iss, result, token)) {
				n_cols++;
			}
		}
	}

	// Create matrix
	matrix B(n_lines, n_cols);

	// Init indices
	int i = 0;	// row index
	int j = 0;	// column index

	// Reopen file
	file.close();
	file.open(file_name);

	// Write values
	while (std::getline(file, line)) {
		std::istringstream iss2(line);
		j = 0;
		while (std::getline(iss2, result, token)) {
			B(i, j) = std::stod(result);
			j++;
		}
		i++;
	}

	// Copy to input matrix
	(*A) = B;

	file.close();
};

void matrix::write_to_file(ofstream* writefile) {

	int i = -1;
	int j = -1;
	for (auto k = begin(); k != end(); k++) {
		j++;
		if (!(j % dim_cols_)) {
			i++;
			j = 0;
		}
		if ((j == 0) && (j == dim_cols_ - 1))
			*writefile << "(" << "\t" << *k << "\t" << ")" << std::endl;
		else if (j == 0)
			*writefile << "(" << "\t" << *k << "\t";
		else if (j == dim_cols_ - 1)
			*writefile << *k << "\t" << ")" << std::endl;
		else
			*writefile << *k << "\t";
	}
}

// Definiere Matrix-Skalar Operatoren
matrix operator+(matrix A, double b) {
	//matrix mat_add1(matrix* A, double b);
	return mat_add1(&A, b);
};
matrix operator+(double b, matrix A) {
	//matrix mat_add1(double b, matrix* A);
	return mat_add1(b, &A);
};
matrix operator-(matrix A, double b) {
	//matrix mat_add1(matrix* A, double b);
	return mat_add1(&A, -b);
};
matrix operator-(double b, matrix A) {
	//matrix mat_add1(matrix* A, double b);
	return mat_sub(b, &A);
};
matrix operator*(matrix A, double b) {
	//matrix mat_mult2(matrix* A, double b);
	return mat_mult2(&A, b);
};
matrix operator*(double b, matrix A) {
	//matrix mat_mult2(double b, matrix* A);
	return mat_mult2(b, &A);
};
matrix operator/(matrix A, double b) {
	//matrix mat_mult2(matrix* A, double b);
	return mat_mult2(&A, 1 / b);
};

matrix operator/(double b, matrix A) {
	return div_mat(b, &A);
};

// Addition: Matrix + Skalar
matrix mat_add1(matrix* A, double b) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Addition
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = (*A)(i, j) + b;
		}
	}
	return result;
};

// Subtraktion skalar - matrix
matrix mat_sub(double b, matrix* A) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Addition
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = b - (*A)(i, j);
		}
	}
	return result;
};

// Addition: Matrix + Skalar mit iterator ... aus irgendeinem Grund langsamer als mat_add1
matrix mat_add2(matrix* A, double b) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Addition
	int i = -1;
	int j = -1;
	for (auto k = (*A).begin(); k != (*A).end(); k++) {
		j++;
		if (!(j % (*A).cols())) {
			i++;
			j = 0;
		}
		(result)(i, j) = *k + b;
	}
	return result;
};

// Elementweise Multiplikation
matrix elementprod(matrix* A, matrix* B) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Addition
	int i = -1;
	int j = -1;
	for (auto k = (*A).begin(), l = (*B).begin(); k != (*A).end(); k++, l++) {
		j++;
		if (!(j % (*A).cols())) {
			i++;
			j = 0;
		}
		(result)(i, j) = *k * *l;
	}
	return result;
}

// Elementweise Division
matrix elementdiv(matrix* A, matrix* B) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Addition
	int i = -1;
	int j = -1;
	for (auto k = (*A).begin(), l = (*B).begin(); k != (*A).end(); k++, l++) {
		j++;
		if (!(j % (*A).cols())) {
			i++;
			j = 0;
		}
		(result)(i, j) = *k / *l;
	}
	return result;
}

// Addition: Skalar + Matrix
matrix mat_add1(double b, matrix* A) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Addition
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = (*A)(i, j) + b;
		}
	}
	return result;
};

// Multiplikation: Matrix * Skalar
matrix mat_mult2(matrix* A, double b) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Multiplikation
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = (*A)(i, j) * b;
		}
	}
	return result;
};

// Multiplikation: Skalar * Matrix
matrix mat_mult2(double b, matrix* A) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Multiplikation
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = (*A)(i, j) * b;
		}
	}
	return result;
};

// Division eines Skalars durch alle Matrix-Elemente
matrix div_mat(double b, matrix* A) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Addition
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = b / (*A)(i, j);
		}
	}
	return result;
};

// Elementweise Exponentialfunktion
matrix exp(matrix* A) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Multiplikation
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = exp((*A)(i, j));
		}
	}
	return result;
};

// Elementweiser Logarithmus
matrix log(matrix* A) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Multiplikation
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = log((*A)(i, j));
		}
	}
	return result;
};

// Elementweise Wurzel
matrix sqrt(matrix* A) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Multiplikation
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = sqrt((*A)(i, j));
		}
	}
	return result;
};

// Elementweise Potenzieren
matrix pow(matrix* A, double d) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Multiplikation
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = pow((*A)(i, j), d);
		}
	}
	return result;
};

// Elementweise invertieren
matrix el_inv(matrix* A) {
	int* size1[2];
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	matrix result(*(*size1), (*(*size1 + 1)));
	// Die eigentliche Multiplikation
	for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size1 + 1); j++) {
			result(i, j) = 1 / ((*A)(i, j));
		}
	}
	return result;
};

// Matrix-Addition zweier übergebener Matrizen: output = input1 * input2 
matrix mat_add(matrix* A, matrix* B) {

	// Check auf Kompatibilität
	int* size1[2];
	int* size2[2];
	int checkval = 0;
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	*size2 = (*B).size();

	matrix result(*(*size1), (*(*size1 + 1)));

	if ((*(*size1 + 1) == *(*size2 + 1)) && (*(*size1) == *(*size2))) {
		// Alles OK
		checkval = 1;
	}
	else {
		std::cout << "Matrizen inkompatibel für Matrizenaddition:" << std::endl
			<< "m(A) = " << *(*size1) << ", n(A) = " << *(*size1 + 1) << std::endl
			<< "m(B) = " << *(*size2) << ", n(B) = " << *(*size2 + 1) << std::endl;
	}

	// Die eigentliche Addition
	if (checkval) {
		/*for (int i = 0; i < *(*size1); i++) {
		for (int j = 0; j < *(*size2 + 1); j++) {
		result(i, j) = (*A)(i, j) + (*B)(i, j);
		}
		}*/
		auto k = (*A).begin();
		auto l = (*B).begin();
		auto m = result.begin();
		for (; k != (*A).end(); k++, l++, m++) {
			*m = (*k + *l);
		}
	}
	//Result.disp();
	return result;
};

matrix mat_eq(matrix* A, matrix* B) {

	// Check auf Kompatibilität
	int* size1[2];
	int* size2[2];
	int checkval = 0;
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	*size2 = (*B).size();

	matrix result(*(*size1), (*(*size1 + 1)));

	if ((*(*size1 + 1) == *(*size2 + 1)) && (*(*size1) == *(*size2))) {
		// Alles OK
		checkval = 1;
	}
	else {
		std::cout << "Matrizen inkompatibel für Matrizenvergleich:" << std::endl
			<< "m(A) = " << *(*size1) << ", n(A) = " << *(*size1 + 1) << std::endl
			<< "m(B) = " << *(*size2) << ", n(B) = " << *(*size2 + 1) << std::endl;
	}

	// Die eigentliche Addition
	if (checkval) {
		auto k = (*A).begin();
		auto l = (*B).begin();
		auto m = result.begin();
		for (; k != (*A).end(); k++, l++, m++) {
			*m = (fabs(*k - *l) < DBL_EPSILON);
		}
	}
	//Result.disp();
	return result;
};

// Matrix-Subtraktion zweier übergebener Matrizen: output = input1 * input2 
matrix mat_sub(matrix* A, matrix* B) {

	// Check auf Kompatibilität
	int* size1[2];
	int* size2[2];
	int checkval = 0;
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	*size2 = (*B).size();

	matrix result(*(*size1), (*(*size1 + 1)));

	if ((*(*size1 + 1) == *(*size2 + 1)) && (*(*size1) == *(*size2))) {
		// Alles OK
		checkval = 1;
	}
	else {
		std::cout << "Matrizen inkompatibel für Matrizenaddition:" << std::endl
			<< "m(A) = " << *(*size1) << ", n(A) = " << *(*size1 + 1) << std::endl
			<< "m(B) = " << *(*size2) << ", n(B) = " << *(*size2 + 1) << std::endl;
	}

	// Die eigentliche Addition
	if (checkval) {
		for (int i = 0; i < *(*size1); i++) {
			for (int j = 0; j < *(*size2 + 1); j++) {
				result(i, j) = (*A)(i, j) - (*B)(i, j);
			}
		}
	}
	//Result.disp();
	return result;
};

// Matrix-Multiplikation zweier übergebener Matrizen: output = input1 * input2 
matrix mat_mult(matrix* A, matrix* B) {

	// Check auf Kompatibilität
	int* size1[2];
	int* size2[2];
	int checkval = 0;
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	*size2 = (*B).size();
	matrix Result(*(*size1), *(*size2 + 1));
	matrix tempB(*(*size2 + 1), *(*size2));
	double res = 0;

	if (*(*size1 + 1) == *(*size2)) {
		// Alles OK
		checkval = 1;
	}
	else {
		std::cout << "Matrizen inkompatibel für Matrixmultiplikation: n(A) = " << *(*size1 + 1) << ", m(B) = " << *(*size2) << std::endl;
	}

	// Die eigentliche Multiplikation
	if (checkval) {
		(*B).mat_transpose(&tempB);

		for (int i = 0; i < *(*size1); i++) {
			for (int j = 0; j < *(*size2 + 1); j++) {
				res = 0;
				for (int k = 0; k < *(*size2); k++) {
					res += (*A)(i, k) * tempB(j, k);
				}
				Result(i, j) = res;
			}
		}
	}
	//std::cout << Result.cols() << std::endl;
	//Result.disp();
	return Result;
};

// Matrix-Multiplikation zweier übergebener Matrizen: Result = input1 * input2 
void mat_mult(matrix* A, matrix* B, matrix* Result) {

	// Check auf Kompatibilität
	int* size1[2];
	int* size2[2];
	int checkval = 0;
	// Alloziere Ergebnismatrix
	*size1 = (*A).size();
	*size2 = (*B).size();
	//matrix Result(*(*size1), *(*size2 + 1));
	matrix tempB(*(*size2 + 1), *(*size2));
	double res = 0;

	if (*(*size1 + 1) == *(*size2)) {
		// Alles OK
		checkval = 1;
	}
	else {
		std::cout << "Matrizen inkompatibel für Matrixmultiplikation: n(A) = " << *(*size1 + 1) << ", m(B) = " << *(*size2) << std::endl;
	}

	// Die eigentliche Multiplikation
	if (checkval) {
		(*B).mat_transpose(&tempB);


		for (int i = 0; i < *(*size1); i++) {
			for (int j = 0; j < *(*size2 + 1); j++) {
				res = 0;
				for (int k = 0; k < *(*size2); k++) {
					res += (*A)(i, k) * tempB(j, k);
				}
				(*Result)(i, j) = res;
			}
		}
	}
	//Result.disp();
	return;
};

// Berechnung einer Jacobi-Matrix
void evalJacobian(double vec[], matrix* Jacobian, int n, int m, void(*fun)(double*, double*), double fval[])
{
	double dx = 1e-6;

	double* temp_vec = new double[n]; // (double*)calloc(n, sizeof(double)); das is C
	double* temp_fval = new double[n]; //(double*)calloc(n, sizeof(double));

	for (int j = 0; j < m; j++) temp_vec[j] = vec[j];

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			temp_vec[j] = vec[j] + dx;
			(*fun)(temp_vec, temp_fval);
			(*fun)(vec, fval);
			(*Jacobian)(i, j) = (temp_fval[i] - fval[i]) / dx;
			temp_vec[j] = vec[j];
			//std::cout << Jacobian(i, j) << "\t";
		}
		//std::cout << std::endl;
	}
	delete[] temp_vec;
	delete[] temp_fval;
	return;
};

/* Standard-Matrixmultiplikation und ikj Algorithmus --> Doppelt so langsam wie obige versionen.
Aufrufe
matrix Cr(n, m);
tstart = std::clock();
mat_mult2(&A, &B, &Cr);
duration = (std::clock() - tstart) / (double)CLOCKS_PER_SEC;
std::cout << "Elapsed time: " << duration << " s" << std::endl;
cout << Cr(1, 1) << ", " << Cr(25, 25) << ", " << Cr(99, 99) << std::endl;
std::cin.get();

matrix Cs(n, m);
tstart = std::clock();
mat_mult3(&A, &B, &Cs);
duration = (std::clock() - tstart) / (double)CLOCKS_PER_SEC;
std::cout << "Elapsed time: " << duration << " s" << std::endl;
cout << Cs(1, 1) << ", " << Cs(25, 25) << ", " << Cs(99, 99) << std::endl;
std::cin.get();
void mat_mult2(matrix* A, matrix* B, matrix* Result) {

// Check auf Kompatibilität
int* size1[2];
int* size2[2];
int checkval = 0;
// Alloziere Ergebnismatrix
*size1 = (*A).size();
*size2 = (*B).size();
//matrix Result(*(*size1), *(*size2 + 1));
matrix tempB(*(*size2 + 1), *(*size2));
double res = 0;
matrix temprowvec(1, *(*size1 + 1));
matrix tempcolvec(*(*size2), 1);

if (*(*size1 + 1) == *(*size2)) {
// Alles OK
checkval = 1;

}
else {
std::cout << "Matrizen inkompatibel für Matrixmultiplikation: n(A) = " << *(*size1 + 1) << ", m(B) = " << *(*size2) << std::endl;
}

// Die eigentliche Multiplikation
if (checkval) {
//(*B).mat_transpose(&tempB);
for (int i = 0; i < *(*size1); i++) {
//(*A).row(i, &temprowvec);
for (int k = 0; k < *(*size2); k++) {
//res = 0;
//(*B).col(j, &tempcolvec);
for (int j = 0; j < *(*size2 + 1); j++) {
(*Result)(i, j) += (*A)(i, k)*(*B)(k, j);
}
}
}
}
//Result.disp();
return;
}

void mat_mult3(matrix* A, matrix* B, matrix* Result) {

// Check auf Kompatibilität
int* size1[2];
int* size2[2];
int checkval = 0;
// Alloziere Ergebnismatrix
*size1 = (*A).size();
*size2 = (*B).size();
//matrix Result(*(*size1), *(*size2 + 1));
matrix tempB(*(*size2 + 1), *(*size2));
double res = 0;
matrix temprowvec(1, *(*size1 + 1));
matrix tempcolvec(*(*size2), 1);

if (*(*size1 + 1) == *(*size2)) {
// Alles OK
checkval = 1;

}
else {
std::cout << "Matrizen inkompatibel für Matrixmultiplikation: n(A) = " << *(*size1 + 1) << ", m(B) = " << *(*size2) << std::endl;
}

// Die eigentliche Multiplikation
if (checkval) {
//(*B).mat_transpose(&tempB);
for (int i = 0; i < *(*size1); i++) {
//(*A).row(i, &temprowvec);
for (int j = 0; j < *(*size2 + 1); j++) {
//res = 0;
//(*B).col(j, &tempcolvec);
for (int k = 0; k < *(*size2); k++) {
(*Result)(i, j) += (*A)(i, k)*(*B)(k, j);
}
}
}
}
//Result.disp();
return;
}*/

// Berechne LR-Zerlegung
void LUPivot(matrix* A, matrix* L, matrix* U, int* piv) {
	int m = (*A).rows();
	int n = (*A).cols();
	int tmp;
	int is;
	(*U) = (*A);
	for (int i = 0; i < n; i++) {
		piv[i] = i;
	}

	for (int i = 0; i < n - 1; i++) {
		is = i;
		bool swap = 0;
		for (int j = is + 1; j < n; j++) {
			if (fabs((*U)(j, i)) > fabs((*U)(is, i))) {
				is = j;
				swap = 1;
			}
		}
		if (swap) {
			(*U).swap_row(i, is);
			(*L).swap_row(i, is);
			tmp = piv[i];
			piv[i] = piv[is];
			piv[is] = tmp;
		}
		for (int j = i + 1; j < n; j++) {
			(*L)(j, i) = (*U)(j, i) / (*U)(i, i);
			for (int k = i; k < n; k++)
				(*U)(j, k) -= (*L)(j, i) * (*U)(i, k);
		}
	}
	for (int i = 0; i < n; i++) {
		(*L)(i, i) = 1;
	}
};

// Berechne Choleskyzerlegung
void GCholesky(matrix* A, matrix* G, bool* errflag) {
	int n = (*A).rows();
	//(*G) = (*A);
	double sum = 0;
	for (int i = 0; i < n; i++) {
		//std::cout << "Iter i = " << i << std::endl;
		for (int j = 0; j <= i; j++) {
			sum = (*A)(i, j);
			//std::cout << "Iter j = " << j << std::endl;
			for (int k = 0; k < j; k++) {
				sum -= (*G)(i, k) * (*G)(j, k);
				//std::cout << "Iter k = " << k <<std::endl;
			}
			//std::cout << "SUMME = " << sum << std::endl;
			if (i > j)
				(*G)(i, j) = sum / (*G)(j, j);
			else if (sum > 0.0)
				(*G)(i, i) = sqrt(sum);
			else {
				std::cout << "Cholesky-ERROR" << std::endl;
				*errflag = 1;
			}
		}
	}
	/*for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	if (i < j)
	(*G)(i, j) = 0;
	}
	}*/
};

// Vorwärts-Rückwärtseinsetzen zum Lösen linearer Gleichungssysteme
void forback(matrix* L, matrix* U, matrix* vec, matrix* x) {
	int n = (*vec).rows();
	double tmp;
	matrix vec_hat(n, 1);
	// Forward elimination: L*vec_hat = vec
	vec_hat(0, 0) = (*vec)(0, 0) / (*L)(0, 0);
	for (int i = 1; i < n; i++) {
		tmp = 0;
		for (int j = 0; j < i; j++) {
			tmp += (*L)(i, j) * vec_hat(j, 0);
		}
		vec_hat(i, 0) = (1 / ((*L)(i, i))) * ((*vec)(i, 0) - tmp);
	}
	// Back substitution: U*x = vec_hat
	(*x)(n - 1, 0) = vec_hat(n - 1, 0) / (*U)(n - 1, n - 1);
	for (int i = n - 2; i >= 0; i--) {
		tmp = 0;
		for (int j = i + 1; j < n; j++) {
			tmp += (*U)(i, j) * (*x)(j, 0);
		}
		(*x)(i, 0) = (1 / ((*U)(i, i))) * (vec_hat(i, 0) - tmp);
	}
	return;
};

// Berechne Determinante auf Basis der LU Zerlegung
double determinant(matrix* A) {
	double det = 1;
	int n = (*A).rows();
	matrix L(n, n);
	matrix U(n, n);
	std::cout << "n = " << n << std::endl;
	std::cin.get();
	int* piv;
	piv = new int[n];
	LUPivot(A, &L, &U, piv);
	U.disp();
	for (int i = 0; i < n; i++)
		det *= U(i, i);
	delete piv;
	return det;
};

// Löse lineares Gleichungssystem A*x = b, teste zuerst Cholesky, dann LR
matrix mldivide(matrix* A, matrix* b) {
	int n = (*b).rows();
	matrix x(n, 1);
	bool err = 0;
	matrix G(n, n);
	matrix GT(n, n);
	GCholesky(A, &G, &err);
	if (!err) {
		G.mat_transpose(&GT);
		forback(&G, &GT, b, &x);
	}
	else {
		matrix L(n, n);
		matrix U(n, n);
		int* piv;
		int* tmppiv;
		piv = new int[n];
		tmppiv = new int[n];
		LUPivot(A, &L, &U, piv);
		int tmpval = 0;
		for (int i = 0; i < n; i++)
			tmppiv[i] = i;
		for (int i = 0; i < n; i++) {
			if (tmppiv[i] == piv[i]) {}
			else {
				int k = 0;
				for (int j = i + 1; j < n; j++) {
					if (tmppiv[j] == piv[i]) {
						k = j;
					}
				}
				tmpval = tmppiv[i];
				tmppiv[i] = tmppiv[k];
				tmppiv[k] = tmpval;
				(*b).swap_row(i, k);
			}
		}
		forback(&L, &U, b, &x);
		delete[] piv;
		delete[] tmppiv;
	}

	return x;
};