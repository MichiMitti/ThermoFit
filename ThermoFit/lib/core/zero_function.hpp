
#ifndef ZERO_FUNCTION_HPP

#define ZERO_FUNCTION_HPP

#pragma once
#include <vector>
#include <cmath>
#include "global.hpp"

class CallClass;	 // forward declaration of CallClass so that the ptr_ in Zero_Function can be of type CallClass*

// Struct for passing solutions between Inner and Outer Problem as vectors of doubles
struct transfer_struct {
	std::vector<std::vector<double>> boundary_inner;	// vector with boundaries for Inner problem
	std::vector<std::vector<double>> boundary_outer;	// vector with boundaries for Outer problem
	std::vector<double> outer;							// Stores solution of the Outer problem
	std::vector<double> inner;							// Stores solution of the Inner problem
	std::vector<double> special_variables;				// Storage for special cases only use if necesarry

	//  forward transformation with the sigmoid function whithin the Inner boundaries
	std::vector<double> sigmoid_forward_Inner(std::vector<double>* vec);
	//  forward transformation with the sigmoid function whithin the Outer boundaries
	std::vector<double> sigmoid_forward_Outer(std::vector<double>* vec);
	//  backward transformation with the sigmoid function whithin the Inner boundaries
	std::vector<double> sigmoid_backward_Inner(std::vector<double>* vec);
	//  backward transformation with the sigmoid function whithin the Outer boundaries
	std::vector<double> sigmoid_backward_Outer(std::vector<double>* vec);

};

// Hawedere nice to meet you!
// I am the legendary Zero_Function, I will help you that your dreams come true and all your problems will be solved!
// My class hold the following properties:
// zero_function ... virtual void function which is overwritten by the problem the function takes two input vectors, the result is written in the second vector
// ptr_ ... void pointer to the CallClass, this is used to call the if necessary the Inner problem
// transfer_ ... is a pointer to a transfer_struct for passing values between Inner and Outer problem
class Zero_Function {
	// Class to solve Problem
public:
	virtual void zero_function(std::vector<double>* vec, std::vector<double>* fval) {}
	virtual void zero_function(double(*vec), double(*fval)) {}
	CallClass* call_class_ptr;						// Pointer to CallClass
	transfer_struct* transfer_;		// Transfer struct pointer to transfer values between Inner and Outer problem
};
#endif // !ZERO_FUNCTION_HPP


