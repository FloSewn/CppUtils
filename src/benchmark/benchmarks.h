/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>
#include <string>

/*********************************************************************
* The main benchmark function
*********************************************************************/
int run_benchmarks(const std::string& library);

/*********************************************************************
* Benchmark functions
*********************************************************************/
void run_benchmarks_QuadTree();
void run_benchmarks_Container();
void run_benchmarks_RTreeND();

