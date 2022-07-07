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

#include "Vec2.h"
#include "Testing.h"


namespace Vec2Tests 
{

using namespace CppUtils;

/*--------------------------------------------------------------------
| Test addition
--------------------------------------------------------------------*/
static void addition() 
{
  const Vec2i a { 1, 1};
  const Vec2i b {-1, 1};
  const Vec2i c = a + b;

  CHECK( (c.x == 0) );
  CHECK( (c.y == 2) );

  Vec2i d = { 3, -3 };

  d += a;

  CHECK( (d.x == 4) );
  CHECK( (d.y == -2) );

} // addition()

/*--------------------------------------------------------------------
| Test scalar product
--------------------------------------------------------------------*/
static void scalar_product() 
{
  const Vec2i a { 2, 3};
  const Vec2i b {-4, 5};

  CHECK( dot(a,b) == 7 );

} // scalar_prodcut()

/*--------------------------------------------------------------------
| Test cross product
--------------------------------------------------------------------*/
static void cross_product() 
{
  const Vec2i a { 2, 3};
  const Vec2i b {-4, 5};

  CHECK( cross(a,b) == 22 );
  
} // cross_product()


} // namespace Vec2Tests


/*********************************************************************
* Run tests for: Vec2.h
*********************************************************************/
void run_tests_Vec2()
{
  Vec2Tests::addition();
  Vec2Tests::scalar_product();
  Vec2Tests::cross_product();

} // run_tests_Vec2()
