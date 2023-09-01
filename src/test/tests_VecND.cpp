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

#include "VecND.h"
#include "Testing.h"
#include "Log.h"
#include "MathUtility.h"

namespace VecNDTests 
{

using namespace CppUtils;

/*--------------------------------------------------------------------
| Test addition
--------------------------------------------------------------------*/
static void addition()
{
  VecND<int,5> a { 1, 2, 3, 4, 5};
  VecND<int,5> b { a };
  VecND<int,5> c = -b;

  VecND<int,5> d = a;
  d += b;

  VecND<int,5> e = a;
  e -= b;

  VecND<int,5> f = a;
  f *= b;

  VecND<int,5> g = a;
  g /= b;

  VecND<int,5> h = a;
  h += 2;

  VecND<int,5> l = h;
  l -= 2;

  VecND<int,5> m = a;
  m *= 2;

  VecND<int,5> n = m;
  n /= 2;

  for ( size_t i = 1; i <= 5; ++i )
  {
    CHECK( a[i-1] == i );
    CHECK( b[i-1] == i );
    CHECK( c[i-1] == -i );
    CHECK( d[i-1] == 2*i );
    CHECK( e[i-1] == 0 );
    CHECK( f[i-1] == i*i );
    CHECK( g[i-1] == 1 );
    CHECK( h[i-1] == i+2 );
    CHECK( l[i-1] == i );
    CHECK( m[i-1] == 2*i );
    CHECK( n[i-1] == i );
  }

  int norm_sqr = a.norm_sqr();
  CHECK( norm_sqr == 55 );
  int norm = a.norm();
  CHECK( norm == 7 );

  CHECK( e.is_zero() );

  CHECK( a == b );

  CHECK( a+e == b );

  CHECK( a/g == b );

  CHECK( a*g == b );

  CHECK( a+2 == h );
  CHECK( 2+a == h );

  CHECK( h-2 == a );
  CHECK( 2-h == -a );

  CHECK( a*2 == d );
  CHECK( 2*a == d );

  CHECK( d/2 == a );

  CHECK( 2/g == 2*g );


  CHECK( a.dot(b) == 55 );
  CHECK( dot(a,b) == 55 );
  CHECK( a.cross(b) == 0 );
  CHECK( cross(a,b) == 0 );


  VecND<int,2> a2 { 1, 2 };
  VecND<int,2> b2 = 2 + a2;

  CHECK( a2.cross(b2) == -2 );
  CHECK( cross(a2, b2) == -2 );
  // In 2D: z should point to x
  CHECK( a2.z == 1 );
  a2.z = 2;
  CHECK( a2.x == 2 );

  // In 1D: z & y should point to x
  VecND<int,1> a1 { 1 };
  CHECK( a1.x == 1 );
  CHECK( a1.y == 1 );
  CHECK( a1.z == 1 );
  a1.y = 2;
  CHECK( a1.x == 2 );
  CHECK( a1.y == 2 );
  CHECK( a1.z == 2 );
  a1.z = 3;
  CHECK( a1.x == 3 );
  CHECK( a1.y == 3 );
  CHECK( a1.z == 3 );


  VecND<int,3> a3 { 1, 2, 3 };
  VecND<int,3> b3 = 2 + a3;

  VecND<int,3> c3 = a3.cross(b3);
  CHECK( c3.x == -2 );
  CHECK( c3.y ==  4 );
  CHECK( c3.z == -2 );

  VecND<int,3> d3 = cross(a3,b3);
  CHECK( d3.x == -2 );
  CHECK( d3.y ==  4 );
  CHECK( d3.z == -2 );


  VecND<int,5>::os_width = 2;
  VecND<int,5>::os_precision = 0;
  LOG(INFO) << a;

  Vec2<int> ai {1,2};



  Vec3d ang_1 { 0.1, 0.4, 0.7 };
  Vec3d ang_2 { 0.23, 0.56, 0.89};

  CHECK( ABS(ang_1.angle(ang_2) - 0.1010) < 1.0E-4 );
  CHECK( ABS(angle(ang_1,ang_2) - 0.1010) < 1.0E-4 );



  CHECK( EQ(ang_1.min(), 0.1) );
  CHECK( EQ(ang_1.max(), 0.7) );

} // addition()

} // namespace VecNDTests

/*********************************************************************
* Run tests for: VecND.h
*********************************************************************/
void run_tests_VecND()
{
  VecNDTests::addition();

} // run_tests_VecND()
