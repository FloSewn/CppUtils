/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <cassert>

#include "tests.h"

#include "Matrix.h"
#include "Testing.h"

namespace MatrixTests 
{
using namespace CppUtils;

/*********************************************************************
* Test Matrix class
*********************************************************************/
void general()
{
  // Create a 10x10 matrix
  Matrix<int> m(10,10);

  CHECK( m.width() == 10 );
  CHECK( m.height() == 10 );

  // Access matrix
  m[0][0] = 1;
  m[1][1] = 1;
  m[2][2] = 1;
  m[3][3] = 1;

  CHECK( m[0][0] == 1 );
  CHECK( m[1][1] == 1 );
  CHECK( m[2][2] == 1 );
  CHECK( m[3][3] == 1 );

  // Copy Matrix
  Matrix n = m;

  CHECK( n[0][0] == 1 );
  CHECK( n[1][1] == 1 );
  CHECK( n[2][2] == 1 );
  CHECK( n[3][3] == 1 );

  // Move Matrix
  Matrix c = std::move(m);

  CHECK( c[0][0] == 1 );
  CHECK( c[1][1] == 1 );
  CHECK( c[2][2] == 1 );
  CHECK( c[3][3] == 1 );

  c[0][0] = 10;
  c[1][1] = 10;
  c[2][2] = 10;
  c[3][3] = 10;

  // Copy back 
  m = c;

  CHECK( m[0][0] == 10 );
  CHECK( m[1][1] == 10 );
  CHECK( m[2][2] == 10 );
  CHECK( m[3][3] == 10 );

  CHECK( c[0][0] == 10 );
  CHECK( c[1][1] == 10 );
  CHECK( c[2][2] == 10 );
  CHECK( c[3][3] == 10 );

  // Swap
  c.swap( n );

  CHECK( c[0][0] == 1 );
  CHECK( c[1][1] == 1 );
  CHECK( c[2][2] == 1 );
  CHECK( c[3][3] == 1 );

  CHECK( n[0][0] == 10 );
  CHECK( n[1][1] == 10 );
  CHECK( n[2][2] == 10 );
  CHECK( n[3][3] == 10 );

} // general()

} // namespace MatrixTests


/*********************************************************************
* Run tests for: Matrix.h
*********************************************************************/
void run_tests_Matrix()
{
  MatrixTests::general();

} // run_tests_Matrix()
