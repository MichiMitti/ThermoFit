#     ThermoFit
<!--
  __________________________________________________________________
 |         ________                              __________         |
 |        /_  __/ /_  ___  _________ ___  ____  / ___ /_/ /_        |
 |         / / / __ \/ _ \/ ___/ __ `__ \/ __ \/ /_  / / __/        |
 |        / / / / / /  __/ /  / / / / / / /_/ / __/ / / /_          |
 |       /_/ /_/ /_/\___/_/  /_/ /_/ /_/\____/_/   /_/\__/          |
 |          The Tool for Optimization in Thermodynamic              |
 |__________________________________________________________________|
-->


## Important Functions

please include these two headers in your problem headers for useing the CallClass!

**"lib/call_function.hpp"** and **"lib/zero_function.hpp"**

for all methods where you can pass a string with "inner" this method is also possible to set it for "outer"!

### Initialization function

This functions shows how the syntax for solveing your problem(s) looks like.
Here you must include your problem header for the test and benchmark functions it is:
**#include "test/call_function_test.hpp"** with the int passed to the function you can select your problem, pass "0" to see the information.

**int  init_call_function(int case_selection)**

### The CallClass constructor

This are the constructors for your problem if you only have a single problem only pass one "Zero_Function*", for two nested problems
please pass first the Inner problem and then the Outer one. You can also pass your boundaries here if you only have one for the Inner
and Outer problem this vector will then be splitted into two vectors for Inner and Outer. If only a one dimensional vector is provided 
this boundaries will be used for all dimensions. 

**CallClass(Zero_Function\* Fun_Inner, std::vector<std::vector<double>>\* boundaries = nullptr)**

**CallClass(Zero_Function\* Fun_Inner, Zero_Function\* Fun_Outer, std::vector<std::vector<double>>\* boundaries = nullptr)**


### The specification method on how you want to solve your problem

Please specify your problems after the constructor with the following to Functions:

Declare method as string to solve the INNER problem options are: "pso", ("powell" coming soon) , "simplex", "ga"
If no boundaries are given the boundaries will be used from the class call or default values will be used this values can be modified by 
changing these two values: shift_for_std_boundary and multiplier_for_std_boundary.
If initGuess_Inner is empty or a nullptr the USER is asked to specify the number of dimensions and the initial guess will be searched randomly
but other methods can be choosen with the set_initalGuess_options() function.
If no options are specified the standart options are used.

**void solve_Inner_problem_with(std::string method_Inner, std::vector<std::vector<double>>\* boundaries = nullptr, std::vector<double>\* initGuess_Inner = nullptr, void\* opts_Inner = nullptr)**


Declare method as string to solve the Outer problem options are: "pso", ("powell" coming soon) , "simplex", "ga"
If no boundaries are given the boundaries will be used from the class call or default values will be used this values can be modified by 
changing these two values: shift_for_std_boundary and multiplier_for_std_boundary.
If initGuess_Inner is empty or a nullptr the USER is asked to specify the number of dimensions and the initial guess will be searched randomly
but other methods can be choosen with the set_initalGuess_options() function. If both guesses were not given the guess will be searched randomly for both problems!
If no options are specified the standart options are used.

**void solve_Outer_problem_with(std::string method_Outer, std::vector<std::vector<double>>\* boundaries = nullptr, std::vector<double>\* initGuess_Outer = nullptr, void\* opts_Outer = nullptr)**


### Solveing methods

Solves the declared problem(s)

**solve()**

Function should be called from Outer Function to solve Inner Function.
The syntax is: **"call_class_ptr->call_inner_problem(call_class_ptr);"** the result is set in the transfer_struct
for this the following two header files must be included in the Inner problem and Outer problem: 

**"lib/call_function.hpp"** and **"lib/zero_function.hpp"**

**void call_inner_problem(Call_Class\* ptr_to_class)**

If you only want to calculate a inital guess randomly you can use these two methods but be aware that some criterions (dimensions, boundaries, transfer_struct_)
must be specified first please look at the init_call function to see how this works.

This method searches an inital guess randomly for a Inner or Outer problem 
For finding Inner guess please set all dependencies most important the inital guess for the outer problem 
set_Outer_initial_guess_for_single_problem(&initGuess_Outer); (the guess will be passed to the transfer_struct_.outer)
For finding Outer guess the dependencies must be set with the following methods:
set_initalGuess_options("pso");	 				(if you want to search it with pso instead of random)
set_Boundary_conditions(&boundaries_inner, "inner");
set_Inner_initial_guess_for_single_problem(&initGuess_Inner);
set_inital_guess_options("inner","simplex"); 			(if you want a diffrent method for solveing the Inner problem)

**void random_guess_search(std::vector<std::vector<double>>\* boundaries, std::string problem = "inner")**

This method searches an inital guess randomly for two specified problems within boudaries randomly (the base can be set with the next method)
(the Inner zero_function will not be solved!)

**void random_guess_search(std::vector<std::vector<double>>\* boundaries_inner, std::vector<std::vector<double>>\* boundaries_outer)**

Sets base for how many guesses for random searched will be used **guesses = base^dimensions**

**void set_base_for_random_guess_search(int base)**


### Available transfer methods

Please use Inner and Outer method in you Problem because therefore the right boundaries will be used!

**std::vector<double> sigmoid_forward_Inner(std::vector<double> vec)**

**std::vector<double> sigmoid_forward_Outer(std::vector<double> vec)**

**std::vector<double> sigmoid_backward_Inner(std::vector<double> vec)**

**std::vector<double> sigmoid_backward_Outer(std::vector<double> vec)**


### Available set-methods

**void set_transfer_struct_special_variables(std::vector<double> special_variables)**

**void set_Boundary_conditions(std::vector<std::vector<double>>\* boundaries = nullptr, std::string problem = "inner")**

**void set_Outer_initial_guess(std::vector<double>\* initGuess_Outer)**

**void set_Inner_initial_guess(std::vector<double>\* initGuess_Inner)**

**void set_inital_guess_options(std::string method = "random", std::string options_Init_Guess = "pso", void\* options_initGuess = nullptr)**

**void set_solveing_method(std::string method = "", std::string problem = "inner")**


### Display methods

Displays results for Inner and Outer problem

**void display_results()**

Displays results for inital guess calculation

**void display_inital_guess_results()**
