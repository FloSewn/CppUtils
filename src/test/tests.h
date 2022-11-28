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
* The main test function
*********************************************************************/
int run_tests(const std::string& library);

/*********************************************************************
* Test functions
*********************************************************************/
void run_tests_MathUtility();
void run_tests_StringOps(); 
void run_tests_Vec2(); 
void run_tests_Geometry(); 
void run_tests_QuadTree();
void run_tests_Container();
void run_tests_ParaReader();
void run_tests_VtkIO();
void run_tests_Log();
void run_tests_Matrix();
void run_tests_BTree();
void run_tests_RTree();
void run_tests_VecND();
void run_tests_BBoxND();


