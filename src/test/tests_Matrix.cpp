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

  CHECK( m.rows() == 10 );
  CHECK( m.columns() == 10 );

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

/*********************************************************************
* Test Matrix class
*********************************************************************/
void elongated()
{
  // Create a 10x2 matrix
  Matrix<int> m(4,2);

  CHECK( m.rows() == 4 );
  CHECK( m.columns() == 2 );

  // Access matrix
  m[0][0] = 1;
  m[0][1] = 2;
  m[1][0] = 3;
  m[1][1] = 4;
  m[2][0] = 5;
  m[2][1] = 6;
  m[3][0] = 7;
  m[3][1] = 8;

  CHECK( m[0][0] == 1 );
  CHECK( m[0][1] == 2 );
  CHECK( m[1][0] == 3 );
  CHECK( m[1][1] == 4 );
  CHECK( m[2][0] == 5 );
  CHECK( m[2][1] == 6 );
  CHECK( m[3][0] == 7 );
  CHECK( m[3][1] == 8 );

  // Copy Matrix
  Matrix n = m;

  CHECK( n[0][0] == 1 );
  CHECK( n[0][1] == 2 );
  CHECK( n[1][0] == 3 );
  CHECK( n[1][1] == 4 );
  CHECK( n[2][0] == 5 );
  CHECK( n[2][1] == 6 );
  CHECK( n[3][0] == 7 );
  CHECK( n[3][1] == 8 );

  // Move Matrix
  Matrix c = std::move(m);

  CHECK( c[0][0] == 1 );
  CHECK( c[0][1] == 2 );
  CHECK( c[1][0] == 3 );
  CHECK( c[1][1] == 4 );
  CHECK( c[2][0] == 5 );
  CHECK( c[2][1] == 6 );
  CHECK( c[3][0] == 7 );
  CHECK( c[3][1] == 8 );

  c[0][1] = 10;
  c[1][1] = 10;
  c[2][1] = 10;
  c[3][1] = 10;

  // Copy back 
  m = c;

  CHECK( m[0][0] == 1 );
  CHECK( m[0][1] == 10);
  CHECK( m[1][0] == 3 );
  CHECK( m[1][1] == 10);
  CHECK( m[2][0] == 5 );
  CHECK( m[2][1] == 10);
  CHECK( m[3][0] == 7 );
  CHECK( m[3][1] == 10);

  // Swap
  c.swap( n );

  CHECK( n[0][0] == 1 );
  CHECK( n[0][1] == 10);
  CHECK( n[1][0] == 3 );
  CHECK( n[1][1] == 10);
  CHECK( n[2][0] == 5 );
  CHECK( n[2][1] == 10);
  CHECK( n[3][0] == 7 );
  CHECK( n[3][1] == 10);

  CHECK( c[0][0] == 1 );
  CHECK( c[0][1] == 2 );
  CHECK( c[1][0] == 3 );
  CHECK( c[1][1] == 4 );
  CHECK( c[2][0] == 5 );
  CHECK( c[2][1] == 6 );
  CHECK( c[3][0] == 7 );
  CHECK( c[3][1] == 8 );


} // elongated()

/*********************************************************************
* Test Matrix class
*********************************************************************/
void resize()
{
  Matrix<int> m;

  CHECK( m.size() == 0 );
  CHECK( m.rows() == 0 );
  CHECK( m.columns() == 0 );

  m.resize(4, 2);

  CHECK( m.rows() == 4 );
  CHECK( m.columns() == 2 );

  // Access matrix
  m[0][0] = 1;
  m[0][1] = 2;
  m[1][0] = 3;
  m[1][1] = 4;
  m[2][0] = 5;
  m[2][1] = 6;
  m[3][0] = 7;
  m[3][1] = 8;

  CHECK( m[0][0] == 1 );
  CHECK( m[0][1] == 2 );
  CHECK( m[1][0] == 3 );
  CHECK( m[1][1] == 4 );
  CHECK( m[2][0] == 5 );
  CHECK( m[2][1] == 6 );
  CHECK( m[3][0] == 7 );
  CHECK( m[3][1] == 8 );

  m.resize(2, 2);

  CHECK( m.rows() == 2 );
  CHECK( m.columns() == 2 );

  CHECK( m[0][0] == 1 );
  CHECK( m[0][1] == 2 );
  CHECK( m[1][0] == 3 );
  CHECK( m[1][1] == 4 );


  // Resize matrix


} // resize()

} // namespace MatrixTests


/*********************************************************************
* Run tests for: Matrix.h
*********************************************************************/
void run_tests_Matrix()
{
  MatrixTests::general();
  MatrixTests::elongated();
  MatrixTests::resize();

} // run_tests_Matrix()
