#include <iostream>
#include <cassert>

#include "tests.h"

#include "Geometry.h"


namespace GeometryTests 
{

using namespace CppUtils;

/*--------------------------------------------------------------------
| Test orientation
--------------------------------------------------------------------*/
static void orientation() 
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { -0.1, -0.1 };
  const Vec2d q1 { -0.2, -0.3 };
  const Vec2d r1 {  0.5,  0.2 };

  CHECK( (orientation(p1,q1,r1) == Orientation::CCW) );
  CHECK( (orientation(p1,r1,q1) == Orientation::CW)  );
  CHECK( (orientation(p1,r1,r1) == Orientation::CL)  );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 { -0.1f, -0.1f };
  const Vec2f q2 { -0.2f, -0.3f };
  const Vec2f r2 {  0.5f,  0.2f };

  CHECK( (orientation(p2,q2,r2) == Orientation::CCW) );
  CHECK( (orientation(p2,r2,q2) == Orientation::CW)  );
  CHECK( (orientation(p2,r2,r2) == Orientation::CL)  );

  /*------------------------------------------------------------------
  | Int
  ------------------------------------------------------------------*/
  const Vec2i p3 { -1, -1 };
  const Vec2i q3 { -2, -3 };
  const Vec2i r3 {  5,  2 };

  CHECK( (orientation(p3,q3,r3) == Orientation::CCW) );
  CHECK( (orientation(p3,r3,q3) == Orientation::CW)  );
  CHECK( (orientation(p3,r3,r3) == Orientation::CL)  );

} // addition()



} // namespace GeometryTests


/*********************************************************************
* Run tests for: Geometry.h
*********************************************************************/
void run_tests_Geometry()
{
  GeometryTests::orientation();

} // run_tests_Geometry()
