/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <memory>
#include <vector>
#include <cassert>

#include <CppUtilsConfig.h>

#include "VecND.h"
#include "Testing.h"
#include "BBoxND.h"

namespace BBoxNDTests 
{
using namespace CppUtils;

/*--------------------------------------------------------------------
| Test constructor
--------------------------------------------------------------------*/
void constructor()
{
  BBoxND<int, 3> a3 { {0,0,0}, {2,2,2} };
  BBoxND<int, 3> b3 { {0,0,0}, {1,1,1} };

  CHECK( b3.bbox_intersection(a3) == a3.bbox_intersection(b3) );
  CHECK( b3.bbox_intersection(a3) == 1 );
  CHECK( b3.bbox_cover(a3) == a3 );

  BBoxND<int, 2> a2 { {0,0}, {2,2} };
  BBoxND<int, 2> b2 { {0,0}, {1,1} };

  CHECK( b2.bbox_intersection(a2) == a2.bbox_intersection(b2) );
  CHECK( b2.bbox_intersection(a2) == 1 );
  CHECK( b2.bbox_cover(a2) == a2 );

  BBoxND<int, 1> a1 { {0}, {2} };
  BBoxND<int, 1> b1 { {0}, {1} };

  CHECK( b1.bbox_intersection(a1) == a1.bbox_intersection(b1) );
  CHECK( b1.bbox_intersection(a1) == 1 );
  CHECK( b1.bbox_cover(a1) == a1 );

  LOG(INFO) << b1;
  LOG(INFO) << b2;
  LOG(INFO) << b3;

} // constructor()

/*--------------------------------------------------------------------
| Test intersection
--------------------------------------------------------------------*/
void intersection()
{
  // d = 2
  BBoxND<int, 2> a1_d2 { {1,1}, {4,4} };
  BBoxND<int, 2> b1_d2 { {4,1}, {6,4} };
  BBoxND<int, 2> b2_d2 { {3,3}, {8,6} };
  VecND<int, 2>  p1_d2 { 5, 2 };
  VecND<int, 2>  p2_d2 { 4, 2 };
  VecND<int, 2>  p3_d2 { 3, 2 };

  // a1_d2 and b1_d2 are adjacent (they touch, but don't intersect)
  CHECK( !a1_d2.bbox_intersect(b1_d2) ); 
  CHECK( a1_d2.bbox_intersect_touch(b1_d2) ); 
  CHECK( a1_d2.bbox_touch(b1_d2) ); 
  CHECK( a1_d2.bbox_intersection(b1_d2) == 0 ); 

  // a1_d2 and b2_d2 intersect
  CHECK( a1_d2.bbox_intersect(b2_d2) ); 
  CHECK( a1_d2.bbox_intersection(b2_d2) == 1 ); 
  CHECK( !a1_d2.bbox_touch(b2_d2) ); 

  // Points
  CHECK( !a1_d2.point_inside( p1_d2 ) );
  CHECK( !a1_d2.point_inside_touch( p1_d2 ) );
  CHECK( !a1_d2.point_touch( p1_d2 ) );

  CHECK( !a1_d2.point_inside( p2_d2 ) );
  CHECK(  a1_d2.point_inside_touch( p2_d2 ) );
  CHECK(  a1_d2.point_touch( p2_d2 ) );

  CHECK(  a1_d2.point_inside( p3_d2 ) );
  CHECK(  a1_d2.point_inside_touch( p3_d2 ) );
  CHECK(  !a1_d2.point_touch( p3_d2 ) );



  // d = 3 
  BBoxND<float, 3> a1_d3 { {1.0f,1.0f,1.0f}, {4.0f,4.0f,4.0f} };
  BBoxND<float, 3> b1_d3 { {4.0f,1.0f,1.0f}, {6.0f,4.0f,2.0f} };
  BBoxND<float, 3> b2_d3 { {3.0f,3.0f,2.0f}, {8.0f,6.0f,3.0f} };
  BBoxND<float, 3> b3_d3 { {0.0f,4.0f,4.0f}, {1.0f,5.0f,5.0f} };

  // a1_d3 and b1_d3 are adjacent (they touch, but don't intersect)
  CHECK( !a1_d3.bbox_intersect(b1_d3) ); 
  CHECK( a1_d3.bbox_intersect_touch(b1_d3) ); 
  CHECK( a1_d3.bbox_touch(b1_d3) ); 
  CHECK( a1_d3.bbox_intersection(b1_d3) == 0 ); 

  // a1_d3 and b2_d3 intersect
  CHECK( a1_d3.bbox_intersect(b2_d3) ); 
  CHECK( a1_d3.bbox_intersection(b2_d3) == 1 ); 
  CHECK( !a1_d3.bbox_touch(b2_d3) ); 

  // a1_d3 and b3_d3 touch at one vertex
  CHECK( !a1_d3.bbox_intersect(b3_d3) ); 
  CHECK( a1_d3.bbox_intersect_touch(b3_d3) ); 
  CHECK( !a1_d3.bbox_intersection(b3_d3) == 1 ); 
  CHECK( a1_d3.bbox_touch(b3_d3) ); 



} // intersection()

} // namespace BBoxNDTests


/*********************************************************************
* Run tests for: BBoxND.h
*********************************************************************/
void run_tests_BBoxND()
{
  BBoxNDTests::constructor();
  BBoxNDTests::intersection();

} // run_tests_BBoxND()