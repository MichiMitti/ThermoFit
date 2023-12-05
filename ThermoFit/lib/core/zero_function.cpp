#include "zero_function.hpp"

std::vector<double> transfer_struct::sigmoid_forward_Inner(std::vector<double>* vec) {

	std::vector<double> vec_for(vec->size(), 0.0);

	for (int i = 0; i < vec->size(); i++) {
		// if both boundaries are "0.0" no transformation is performed
		if ((boundary_inner[i][1] == 0.0) && (boundary_inner[i][0] == 0.0)) {		// == dangerous! double is not able to store certain values exactly "-" even more dangerous
			vec_for[i] = (*vec)[i];
		}
		// transformation is performed (0...lower bound, 1...higher bound)
		else {
			vec_for[i] = -std::log(((boundary_inner[i][1] - boundary_inner[i][0]) / ((*vec)[i] - boundary_inner[i][0])) - 1.0);
		}
	}
	return vec_for;
}

std::vector<double> transfer_struct::sigmoid_forward_Outer(std::vector<double>* vec) {

	std::vector<double> vec_for(vec->size(), 0.0);

	for (int i = 0; i < vec->size(); i++) {
		// if both boundaries are "0.0" no transformation is performed!
		if ((boundary_outer[i][1] == 0.0) && (boundary_outer[i][0] == 0.0)) {		// == dangerous! double is not able to store certain values exactly "-" even more dangerous
			vec_for[i] = (*vec)[i];
		}
		// transformation is performed (0...lower bound, 1...higher bound)
		else {
			vec_for[i] = -std::log(((boundary_outer[i][1] - boundary_outer[i][0]) / ((*vec)[i] - boundary_outer[i][0])) - 1.0);
		}
	}
	return vec_for;
}

std::vector<double> transfer_struct::sigmoid_backward_Inner(std::vector<double>* vec) {

	std::vector<double> vec_back(vec->size(), 0.0);

	for (int i = 0; i < vec->size(); i++) {
		// if both boundaries are "0.0" no transformation is performed
		if ((boundary_inner[i][1] == 0.0) && (boundary_inner[i][0] == 0.0)) {		// == dangerous! double is not able to store certain values exactly "-" even more dangerous
			vec_back[i] = (*vec)[i];
		}
		// transformation is performed (0...lower bound, 1...higher bound)
		else {
			vec_back[i] = boundary_inner[i][0] + (boundary_inner[i][1] - boundary_inner[i][0]) / (1.0 + std::exp(-(*vec)[i]));
		}
	}
	return vec_back;
}

std::vector<double> transfer_struct::sigmoid_backward_Outer(std::vector<double>* vec) {

	std::vector<double> vec_back(vec->size(), 0.0);

	for (int i = 0; i < vec->size(); i++) {
		// if both boundaries are "0.0" no transformation is performed
		if ((boundary_outer[i][1] == 0.0) && (boundary_outer[i][0] == 0.0)) {		// == dangerous! double is not able to store certain values exactly "-" even more dangerous
			vec_back[i] = (*vec)[i];
		}
		// transformation is performed (0...lower bound, 1...higher bound)
		else {
			vec_back[i] = boundary_outer[i][0] + (boundary_outer[i][1] - boundary_outer[i][0]) / (1.0 + std::exp(-(*vec)[i]));
		}
	}
	return vec_back;
}