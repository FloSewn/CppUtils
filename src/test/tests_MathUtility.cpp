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

#include "MathUtility.h"

namespace MathUtilityTests 
{

using namespace CppUtils;

/*--------------------------------------------------------------------
| Test abs
--------------------------------------------------------------------*/
static void abs()
{
  // Integers
  int a = -5;
  CHECK( (ABS(a) == 5) );
  CHECK( (ABS(a) != a) );

  // Float
  float b = -5.0f;
  CHECK( (ABS(b) == 5.0f) );
  CHECK( (ABS(b) != b) );

  // Double
  double c = -5.0;
  CHECK( (ABS(c) == 5.0) );
  CHECK( (ABS(c) != c) );
}

/*--------------------------------------------------------------------
| Test eq
--------------------------------------------------------------------*/
static void eq()
{
  // Integers
  int a = -5;
  CHECK(  EQ(a, a) );
  CHECK( !EQ(a, 5) );

  // Float
  float b = -5.0f;
  CHECK(  EQ(b, b) );
  CHECK( !EQ(b, 5.0f) );

  // Double
  double c = -5.0;
  CHECK(  EQ(c, c) );
  CHECK( !EQ(c, 5.0) );
}

/*--------------------------------------------------------------------
| Test eq0
--------------------------------------------------------------------*/
static void eq0()
{
  // Integers
  CHECK(  EQ0(1-1) );
  CHECK( !EQ0(1) );

  // Float
  CHECK(  EQ0(1.0f-1.0f) );
  CHECK( !EQ0(1.0f) );

  // Double
  CHECK(  EQ0(1.0-1.0) );
  CHECK( !EQ0(1.0) );
}

/*--------------------------------------------------------------------
| Test min
--------------------------------------------------------------------*/
static void min()
{
  // Integers
  CHECK( ( MIN(1,-1) == -1 ) );

  // Float
  CHECK( EQ( MIN(1.0f,-1.0f), -1.0f ) );

  // Double
  CHECK( EQ( MIN(1.0,-1.0), -1.0 ) );
}

/*--------------------------------------------------------------------
| Test max
--------------------------------------------------------------------*/
static void max()
{
  // Integers
  CHECK( ( MAX(1,-1) == 1 ) );

  // Float
  CHECK( EQ( MAX(1.0f,-1.0f), 1.0f ) );

  // Double
  CHECK( EQ( MAX(1.0,-1.0), 1.0 ) );
}

/*--------------------------------------------------------------------
| Test mod
--------------------------------------------------------------------*/
static void mod()
{
  CHECK( MOD(4, 3) == 1 );
}

} // namespace MathUtilityTests 


/*********************************************************************
* Run tests for: MathUtility.h
*********************************************************************/
void run_tests_MathUtility()
{
  MathUtilityTests::abs();
  MathUtilityTests::eq();
  MathUtilityTests::eq0();
  MathUtilityTests::min();
  MathUtilityTests::max();
  MathUtilityTests::mod();

} // run_tests_MathUtility()
