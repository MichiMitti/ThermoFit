/*
Initialize CallClass function file

@author: Michael Mitterlindner, TU Graz
@email: mitterlindner@gmx.at
@github: MichiMitti
@year: 2023
*/

#pragma once
#include <vector>
#include <cmath>
#include "../lib/core/global.hpp"

// Includes call function 
#include "../lib/call_function.hpp"

// Include used problem files:
#include "test/test_functions.hpp"

// Initializes the problem(s) and pases them to the call function
// 0....Information
// 1....Single problem selected
// 2....Two problems selected
// 11...Benchmark test functions
// 21...Only inital guess for Inner problem is searched randomly
// 22...Only inital guess for Outer problem is searched randomly
// 23...Inital guess searched randomly for Inner and Outer probem
int  init_call_function(int case_selection);
