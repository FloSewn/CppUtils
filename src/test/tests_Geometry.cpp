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

} // orientation()

/*********************************************************************
* Test is_left()
*********************************************************************/
void is_left()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { -0.1, -0.1 };
  const Vec2d q1 {  0.1,  0.1 };
  const Vec2d r1 { -0.1,  0.1 };
  const Vec2d s1 {  0.1, -0.1 };
  const Vec2d t1 {  0.0, -0.0 };

  CHECK( ( is_left(p1,q1,r1) ) );
  CHECK( !(is_left(p1,q1,s1) ) );
  CHECK( !(is_left(p1,q1,t1) ) );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 { -0.1f, -0.1f };
  const Vec2f q2 {  0.1f,  0.1f };
  const Vec2f r2 { -0.1f,  0.1f };
  const Vec2f s2 {  0.1f, -0.1f };
  const Vec2f t2 {  0.0f, -0.0f };

  CHECK( ( is_left(p2,q2,r2) ) );
  CHECK( !(is_left(p2,q2,s2) ) );
  CHECK( !(is_left(p2,q2,t2) ) );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 { -1, -1 };
  const Vec2i q3 {  1,  1 };
  const Vec2i r3 { -1,  1 };
  const Vec2i s3 {  1, -1 };
  const Vec2i t3 {  0, -0 };

  CHECK( ( is_left(p3,q3,r3) ) );
  CHECK( !(is_left(p3,q3,s3) ) );
  CHECK( !(is_left(p3,q3,t3) ) );

} // is_left()

/*********************************************************************
* Test is_lefton()
*********************************************************************/
void is_lefton()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { -0.1, -0.1 };
  const Vec2d q1 {  0.1,  0.1 };
  const Vec2d r1 { -0.1,  0.1 };
  const Vec2d s1 {  0.1, -0.1 };
  const Vec2d t1 {  0.0,  0.0 };

  CHECK(  (is_lefton(p1,q1,r1) ) );
  CHECK( !(is_lefton(p1,q1,s1) ) );
  CHECK(  (is_lefton(p1,q1,t1) ) );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2d p2 { -0.1f, -0.1f };
  const Vec2d q2 {  0.1f,  0.1f };
  const Vec2d r2 { -0.1f,  0.1f };
  const Vec2d s2 {  0.1f, -0.1f };
  const Vec2d t2 {  0.0f,  0.0f };

  CHECK(  (is_lefton(p2,q2,r2) ) );
  CHECK( !(is_lefton(p2,q2,s2) ) );
  CHECK(  (is_lefton(p2,q2,t2) ) );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2d p3 { -1, -1 };
  const Vec2d q3 {  1,  1 };
  const Vec2d r3 { -1,  1 };
  const Vec2d s3 {  1, -1 };
  const Vec2d t3 {  0,  0 };

  CHECK(  (is_lefton(p3,q3,r3) ) );
  CHECK( !(is_lefton(p3,q3,s3) ) );
  CHECK(  (is_lefton(p3,q3,t3) ) );

} // is_lefton()

/*********************************************************************
* Test in_segment()
*********************************************************************/
void in_segment()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { -0.1, -0.1 };
  const Vec2d q1 {  0.1,  0.1 };
  const Vec2d r1 { -0.1,  0.1 };
  const Vec2d s1 {  0.1, -0.1 };
  const Vec2d t1 {  0.2,  0.2 };
  const Vec2d u1 { -0.2, -0.2 };
  const Vec2d v1 {  0.0,  0.0 };

  CHECK( !(in_segment(p1,q1,r1) ) );
  CHECK( !(in_segment(p1,q1,s1) ) );
  CHECK( !(in_segment(p1,q1,t1) ) );
  CHECK( !(in_segment(p1,q1,u1) ) );
  CHECK(  (in_segment(p1,q1,v1) ) );
  CHECK( !(in_segment(p1,q1,p1) ) );
  CHECK( !(in_segment(p1,q1,q1) ) );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 { -0.1f, -0.1f };
  const Vec2f q2 {  0.1f,  0.1f };
  const Vec2f r2 { -0.1f,  0.1f };
  const Vec2f s2 {  0.1f, -0.1f };
  const Vec2f t2 {  0.2f,  0.2f };
  const Vec2f u2 { -0.2f, -0.2f };
  const Vec2f v2 {  0.0f,  0.0f };

  CHECK( !(in_segment(p2,q2,r2) ) );
  CHECK( !(in_segment(p2,q2,s2) ) );
  CHECK( !(in_segment(p2,q2,t2) ) );
  CHECK( !(in_segment(p2,q2,u2) ) );
  CHECK(  (in_segment(p2,q2,v2) ) );
  CHECK( !(in_segment(p2,q2,p2) ) );
  CHECK( !(in_segment(p2,q2,q2) ) );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 { -1, -1 };
  const Vec2i q3 {  1,  1 };
  const Vec2i r3 { -1,  1 };
  const Vec2i s3 {  1, -1 };
  const Vec2i t3 {  2,  2 };
  const Vec2i u3 { -2, -2 };
  const Vec2i v3 {  0,  0 };

  CHECK( !(in_segment(p3,q3,r3) ) );
  CHECK( !(in_segment(p3,q3,s3) ) );
  CHECK( !(in_segment(p3,q3,t3) ) );
  CHECK( !(in_segment(p3,q3,u3) ) );
  CHECK(  (in_segment(p3,q3,v3) ) );
  CHECK( !(in_segment(p3,q3,p3) ) );
  CHECK( !(in_segment(p3,q3,q3) ) );

} // in_segment()

/*********************************************************************
* Test in_on_segment()
*********************************************************************/
void in_on_segment()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { -0.1, -0.1 };
  const Vec2d q1 {  0.1,  0.1 };
  const Vec2d r1 { -0.1,  0.1 };
  const Vec2d s1 {  0.1, -0.1 };
  const Vec2d t1 {  0.2,  0.2 };
  const Vec2d u1 { -0.2, -0.2 };
  const Vec2d v1 {  0.0,  0.0 };

  CHECK( !(in_on_segment(p1,q1,r1) ) );
  CHECK( !(in_on_segment(p1,q1,s1) ) );
  CHECK( !(in_on_segment(p1,q1,t1) ) );
  CHECK( !(in_on_segment(p1,q1,u1) ) );
  CHECK(  (in_on_segment(p1,q1,v1) ) );
  CHECK(  (in_on_segment(p1,q1,p1) ) );
  CHECK(  (in_on_segment(p1,q1,q1) ) );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 { -0.1f, -0.1f };
  const Vec2f q2 {  0.1f,  0.1f };
  const Vec2f r2 { -0.1f,  0.1f };
  const Vec2f s2 {  0.1f, -0.1f };
  const Vec2f t2 {  0.2f,  0.2f };
  const Vec2f u2 { -0.2f, -0.2f };
  const Vec2f v2 {  0.0f,  0.0f };

  CHECK( !(in_on_segment(p2,q2,r2) ) );
  CHECK( !(in_on_segment(p2,q2,s2) ) );
  CHECK( !(in_on_segment(p2,q2,t2) ) );
  CHECK( !(in_on_segment(p2,q2,u2) ) );
  CHECK(  (in_on_segment(p2,q2,v2) ) );
  CHECK(  (in_on_segment(p2,q2,p2) ) );
  CHECK(  (in_on_segment(p2,q2,q2) ) );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 { -1, -1 };
  const Vec2i q3 {  1,  1 };
  const Vec2i r3 { -1,  1 };
  const Vec2i s3 {  1, -1 };
  const Vec2i t3 {  2,  2 };
  const Vec2i u3 { -2, -2 };
  const Vec2i v3 {  0,  0 };

  CHECK( !(in_on_segment(p3,q3,r3) ) );
  CHECK( !(in_on_segment(p3,q3,s3) ) );
  CHECK( !(in_on_segment(p3,q3,t3) ) );
  CHECK( !(in_on_segment(p3,q3,u3) ) );
  CHECK(  (in_on_segment(p3,q3,v3) ) );
  CHECK(  (in_on_segment(p3,q3,p3) ) );
  CHECK(  (in_on_segment(p3,q3,q3) ) );

} // in_on_segment()

/*********************************************************************
* Test line_line_intersection
*********************************************************************/
void line_line_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 {-0.1, -0.1};
  const Vec2d q1 { 0.1,  0.1};

  const Vec2d r1 {-0.1,  0.1};
  const Vec2d s1 { 0.1, -0.1};

  const Vec2d t1 { 0.0,  0.0};

  const Vec2d u1 {-0.1,  0.0};
  const Vec2d v1 { 0.1,  0.2};

  CHECK( !( line_line_intersection(p1,q1,p1,q1) ) );
  CHECK(  ( line_line_intersection(p1,q1,r1,s1) ) );
  CHECK(  ( line_line_intersection(p1,q1,t1,q1) ) );
  CHECK( !( line_line_intersection(p1,q1,u1,v1) ) );
  CHECK( !( line_line_intersection(p1,q1,q1,v1) ) );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 {-0.1f, -0.1f};
  const Vec2f q2 { 0.1f,  0.1f};

  const Vec2f r2 {-0.1f,  0.1f};
  const Vec2f s2 { 0.1f, -0.1f};

  const Vec2f t2 { 0.0f,  0.0f};

  const Vec2f u2 {-0.1f,  0.0f};
  const Vec2f v2 { 0.1f,  0.2f};

  CHECK( !( line_line_intersection(p2,q2,p2,q2) ) );
  CHECK(  ( line_line_intersection(p2,q2,r2,s2) ) );
  CHECK(  ( line_line_intersection(p2,q2,t2,q2) ) );
  CHECK( !( line_line_intersection(p2,q2,u2,v2) ) );
  CHECK( !( line_line_intersection(p2,q2,q2,v2) ) );


  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 {-1, -1};
  const Vec2i q3 { 1,  1};

  const Vec2i r3 {-1,  1};
  const Vec2i s3 { 1, -1};

  const Vec2i t3 { 0,  0};

  const Vec2i u3 {-1,  0};
  const Vec2i v3 { 1,  2};


  CHECK( !( line_line_intersection(p3,q3,p3,q3) ) );
  CHECK(  ( line_line_intersection(p3,q3,r3,s3) ) );
  CHECK(  ( line_line_intersection(p3,q3,t3,q3) ) );
  CHECK( !( line_line_intersection(p3,q3,u3,v3) ) );
  CHECK( !( line_line_intersection(p3,q3,q3,v3) ) );


} // line_line_intersection() */

/*********************************************************************
* Test line_tri_intersection
*********************************************************************/
void line_tri_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 {-0.1, -0.1};
  const Vec2d q1 { 0.1,  0.1};

  const Vec2d r1 {-0.1,  0.1};
  const Vec2d s1 { 0.1, -0.1};
  const Vec2d t1 { 0.2,  0.2};

  const Vec2d u1 { 0.2, -0.1};
  const Vec2d v1 { 0.3,  0.1};
  const Vec2d w1 { 0.0,  0.0};

  CHECK(  ( line_tri_intersection(p1,q1, r1,s1,t1) ) );
  CHECK(  ( line_tri_intersection(p1,q1, r1,t1,s1) ) );
  CHECK( !( line_tri_intersection(u1,v1, r1,t1,s1) ) );
  CHECK( !( line_tri_intersection(r1,s1, r1,t1,s1) ) );
  CHECK( !( line_tri_intersection(t1,v1, r1,t1,s1) ) );
  CHECK(  ( line_tri_intersection(t1,p1, r1,t1,s1) ) );
  CHECK(  ( line_tri_intersection(r1,w1, r1,t1,s1) ) );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 {-0.1f, -0.1f};
  const Vec2f q2 { 0.1f,  0.1f};

  const Vec2f r2 {-0.1f,  0.1f};
  const Vec2f s2 { 0.1f, -0.1f};
  const Vec2f t2 { 0.2f,  0.2f};

  const Vec2f u2 { 0.2f, -0.1f};
  const Vec2f v2 { 0.3f,  0.1f};
  const Vec2f w2 { 0.0f,  0.0f};

  CHECK(  ( line_tri_intersection(p2,q2, r2,s2,t2) ) );
  CHECK(  ( line_tri_intersection(p2,q2, r2,t2,s2) ) );
  CHECK( !( line_tri_intersection(u2,v2, r2,t2,s2) ) );
  CHECK( !( line_tri_intersection(r2,s2, r2,t2,s2) ) );
  CHECK( !( line_tri_intersection(t2,v2, r2,t2,s2) ) );
  CHECK(  ( line_tri_intersection(t2,p2, r2,t2,s2) ) );
  CHECK(  ( line_tri_intersection(r2,w2, r2,t2,s2) ) );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 {-1, -1};
  const Vec2i q3 { 1,  1};

  const Vec2i r3 {-1,  1};
  const Vec2i s3 { 1, -1};
  const Vec2i t3 { 2,  2};

  const Vec2i u3 { 2, -1};
  const Vec2i v3 { 3,  1};
  const Vec2i w3 { 0,  0};

  CHECK(  ( line_tri_intersection(p3,q3, r3,s3,t3) ) );
  CHECK(  ( line_tri_intersection(p3,q3, r3,t3,s3) ) );
  CHECK( !( line_tri_intersection(u3,v3, r3,t3,s3) ) );
  CHECK( !( line_tri_intersection(r3,s3, r3,t3,s3) ) );
  CHECK( !( line_tri_intersection(t3,v3, r3,t3,s3) ) );
  CHECK(  ( line_tri_intersection(t3,p3, r3,t3,s3) ) );
  CHECK(  ( line_tri_intersection(r3,w3, r3,t3,s3) ) );

} // line_tri_intersection()

/*********************************************************************
* Test line_quad_intersection
*********************************************************************/
void line_quad_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 {-0.1, -0.1};
  const Vec2d q1 { 0.1, -0.1};
  const Vec2d r1 { 0.1,  0.1};
  const Vec2d s1 {-0.1,  0.1};

  const Vec2d t1 { 0.0,  0.0};
  const Vec2d u1 { 0.2,  0.1};
  const Vec2d v1 { 0.2, -0.1};

  CHECK(  ( line_quad_intersection(t1,u1, p1,q1,r1,s1) ) );
  CHECK(  ( line_quad_intersection(t1,u1, s1,r1,q1,p1) ) );
  CHECK( !( line_quad_intersection(u1,v1, p1,q1,r1,s1) ) );
  CHECK( !( line_quad_intersection(r1,u1, p1,q1,r1,s1) ) );
  CHECK( !( line_quad_intersection(r1,s1, p1,q1,r1,s1) ) );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 {-0.1f, -0.1f};
  const Vec2f q2 { 0.1f, -0.1f};
  const Vec2f r2 { 0.1f,  0.1f};
  const Vec2f s2 {-0.1f,  0.1f};

  const Vec2f t2 { 0.0f,  0.0f};
  const Vec2f u2 { 0.2f,  0.1f};
  const Vec2f v2 { 0.2f, -0.1f};

  CHECK(  ( line_quad_intersection(t2,u2, p2,q2,r2,s2) ) );
  CHECK(  ( line_quad_intersection(t2,u2, s2,r2,q2,p2) ) );
  CHECK( !( line_quad_intersection(u2,v2, p2,q2,r2,s2) ) );
  CHECK( !( line_quad_intersection(r2,u2, p2,q2,r2,s2) ) );
  CHECK( !( line_quad_intersection(r2,s2, p2,q2,r2,s2) ) );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 {-1, -1};
  const Vec2i q3 { 1, -1};
  const Vec2i r3 { 1,  1};
  const Vec2i s3 {-1,  1};

  const Vec2i t3 { 0,  0};
  const Vec2i u3 { 2,  1};
  const Vec2i v3 { 2, -1};

  CHECK(  ( line_quad_intersection(t3,u3, p3,q3,r3,s3) ) );
  CHECK(  ( line_quad_intersection(t3,u3, s3,r3,q3,p3) ) );
  CHECK( !( line_quad_intersection(u3,v3, p3,q3,r3,s3) ) );
  CHECK( !( line_quad_intersection(r3,u3, p3,q3,r3,s3) ) );
  CHECK( !( line_quad_intersection(r3,s3, p3,q3,r3,s3) ) );

} // line_quad_intersection()

/*********************************************************************
* Test tri_tri_intersection
*********************************************************************/
void tri_tri_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { 3.0, 2.0 };
  const Vec2d q1 { 9.0, 2.0 };
  const Vec2d r1 { 6.0, 6.0 };

  const Vec2d s1 { 3.0, 2.0 };
  const Vec2d t1 { 6.0, 2.0 };
  const Vec2d u1 { 5.0,-1.0 };

  const Vec2d v1 { 6.0, 4.0 };
  const Vec2d w1 {11.0, 4.0 };
  const Vec2d x1 {10.0, 8.0 };

  const Vec2d a1 { 7.0, 2.0 };
  const Vec2d b1 { 7.0,-1.0 };
  const Vec2d c1 { 9.0,-1.0 };

  const Vec2d d1 {-4.0, 2.0 };
  const Vec2d e1 {-2.0, 2.0 };
  const Vec2d f1 {-3.0, 4.0 };

  const Vec2d g1 { 3.0, 2.0 };
  const Vec2d h1 { 6.0, 6.0 };
  const Vec2d i1 { 2.0, 6.0 };

  CHECK( (tri_tri_intersection( p1,q1,r1, s1,t1,u1 )) );
  CHECK( (tri_tri_intersection( p1,q1,r1, v1,w1,x1 )) );
  CHECK( (tri_tri_intersection( p1,q1,r1, a1,b1,c1 )) );
  CHECK(!(tri_tri_intersection( p1,q1,r1, d1,e1,f1 )) );
  CHECK(!(tri_tri_intersection( p1,q1,r1, g1,h1,i1 )) );
  CHECK(!(tri_tri_intersection( p1,q1,r1, p1,q1,r1 )) );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 { 3.0f, 2.0f };
  const Vec2f q2 { 9.0f, 2.0f };
  const Vec2f r2 { 6.0f, 6.0f };

  const Vec2f s2 { 3.0f, 2.0f };
  const Vec2f t2 { 6.0f, 2.0f };
  const Vec2f u2 { 5.0f,-1.0f };

  const Vec2f v2 { 6.0f, 4.0f };
  const Vec2f w2 {11.0f, 4.0f };
  const Vec2f x2 {10.0f, 8.0f };

  const Vec2f a2 { 7.0f, 2.0f };
  const Vec2f b2 { 7.0f,-1.0f };
  const Vec2f c2 { 9.0f,-1.0f };

  const Vec2f d2 {-4.0f, 2.0f };
  const Vec2f e2 {-2.0f, 2.0f };
  const Vec2f f2 {-3.0f, 4.0f };

  const Vec2f g2 { 3.0f, 2.0f };
  const Vec2f h2 { 6.0f, 6.0f };
  const Vec2f i2 { 2.0f, 6.0f };

  CHECK( (tri_tri_intersection( p2,q2,r2, s2,t2,u2 )) );
  CHECK( (tri_tri_intersection( p2,q2,r2, v2,w2,x2 )) );
  CHECK( (tri_tri_intersection( p2,q2,r2, a2,b2,c2 )) );
  CHECK(!(tri_tri_intersection( p2,q2,r2, d2,e2,f2 )) );
  CHECK(!(tri_tri_intersection( p2,q2,r2, g2,h2,i2 )) );
  CHECK(!(tri_tri_intersection( p2,q2,r2, p2,q2,r2 )) );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 { 3, 2 };
  const Vec2i q3 { 9, 2 };
  const Vec2i r3 { 6, 6 };

  const Vec2i s3 { 3, 2 };
  const Vec2i t3 { 6, 2 };
  const Vec2i u3 { 5,-1 };

  const Vec2i v3 { 6, 4 };
  const Vec2i w3 {11, 4 };
  const Vec2i x3 {10, 8 };

  const Vec2i a3 { 7, 2 };
  const Vec2i b3 { 7,-1 };
  const Vec2i c3 { 9,-1 };

  const Vec2i d3 {-4, 2 };
  const Vec2i e3 {-2, 2 };
  const Vec2i f3 {-3, 4 };

  const Vec2i g3 { 3, 2 };
  const Vec2i h3 { 6, 6 };
  const Vec2i i3 { 2, 6 };

  CHECK( (tri_tri_intersection( p3,q3,r3, s3,t3,u3 )) );
  CHECK( (tri_tri_intersection( p3,q3,r3, v3,w3,x3 )) );
  CHECK( (tri_tri_intersection( p3,q3,r3, a3,b3,c3 )) );
  CHECK(!(tri_tri_intersection( p3,q3,r3, d3,e3,f3 )) );
  CHECK(!(tri_tri_intersection( p3,q3,r3, g3,h3,i3 )) );
  CHECK(!(tri_tri_intersection( p3,q3,r3, p3,q3,r3 )) );

} // tri_tri_intersection()

/*********************************************************************
* Test quad_quad_intersection
*********************************************************************/
void quad_quad_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1d {-0.5,-0.5 };
  const Vec2d q1d { 0.5,-0.5 };
  const Vec2d r1d { 0.5, 0.5 };
  const Vec2d s1d {-0.5, 0.5 };

  // Regular intersecion
  const Vec2d p2d { 0.0,-1.0 };
  const Vec2d q2d { 1.0,-1.0 };
  const Vec2d r2d { 1.0, 0.0 };
  const Vec2d s2d { 0.0, 0.0 };

  // No intersection
  const Vec2d p3d {-2.5, 0.5 };
  const Vec2d q3d {-2.0, 0.0 };
  const Vec2d r3d {-1.5, 0.0 };
  const Vec2d s3d {-1.5, 1.0 };

  // Adjacent edge --> no intersection
  const Vec2d p4d {-0.5, 0.5 };
  const Vec2d q4d { 0.5, 0.5 };
  const Vec2d r4d { 0.5, 1.5 };
  const Vec2d s4d {-0.5, 1.5 };

  // Adjacent vertex --> no intersection
  const Vec2d p5d { 0.5, 0.5 };
  const Vec2d q5d { 1.5, 0.5 };
  const Vec2d r5d { 1.5, 1.5 };
  const Vec2d s5d { 0.5, 1.5 };

  // Half-Edge intersection
  const Vec2d p6d { 0.5,-0.5 };
  const Vec2d q6d { 1.0,-0.5 };
  const Vec2d r6d { 1.0, 0.0 };
  const Vec2d s6d { 0.5, 0.0 };

  CHECK( (quad_quad_intersection(p1d, q1d, r1d, s1d, p2d, q2d, r2d, s2d)) );
  CHECK(!(quad_quad_intersection(p1d, q1d, r1d, s1d, p3d, q3d, r3d, s3d)) );
  CHECK(!(quad_quad_intersection(p1d, q1d, r1d, s1d, p4d, q4d, r4d, s4d)) );
  CHECK(!(quad_quad_intersection(p1d, q1d, r1d, s1d, p5d, q5d, r5d, s5d)) );
  CHECK( (quad_quad_intersection(p1d, q1d, r1d, s1d, p6d, q6d, r6d, s6d)) );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p1f {-0.5f,-0.5f };
  const Vec2f q1f { 0.5f,-0.5f };
  const Vec2f r1f { 0.5f, 0.5f };
  const Vec2f s1f {-0.5f, 0.5f };

  // Regular intersecion
  const Vec2f p2f { 0.0f,-1.0f };
  const Vec2f q2f { 1.0f,-1.0f };
  const Vec2f r2f { 1.0f, 0.0f };
  const Vec2f s2f { 0.0f, 0.0f };

  // No intersection
  const Vec2f p3f {-2.5f, 0.5f };
  const Vec2f q3f {-2.0f, 0.0f };
  const Vec2f r3f {-1.5f, 0.0f };
  const Vec2f s3f {-1.5f, 1.0f };

  // Adjacent edge --> no intersection
  const Vec2f p4f {-0.5f, 0.5f };
  const Vec2f q4f { 0.5f, 0.5f };
  const Vec2f r4f { 0.5f, 1.5f };
  const Vec2f s4f {-0.5f, 1.5f };

  // Adjacent vertex --> no intersection
  const Vec2f p5f { 0.5f, 0.5f };
  const Vec2f q5f { 1.5f, 0.5f };
  const Vec2f r5f { 1.5f, 1.5f };
  const Vec2f s5f { 0.5f, 1.5f };

  // Half-Edge intersection
  const Vec2f p6f { 0.5f,-0.5f };
  const Vec2f q6f { 1.0f,-0.5f };
  const Vec2f r6f { 1.0f, 0.0f };
  const Vec2f s6f { 0.5f, 0.0f };

  CHECK( (quad_quad_intersection(p1f, q1f, r1f, s1f, p2f, q2f, r2f, s2f)) );
  CHECK(!(quad_quad_intersection(p1f, q1f, r1f, s1f, p3f, q3f, r3f, s3f)) );
  CHECK(!(quad_quad_intersection(p1f, q1f, r1f, s1f, p4f, q4f, r4f, s4f)) );
  CHECK(!(quad_quad_intersection(p1f, q1f, r1f, s1f, p5f, q5f, r5f, s5f)) );
  CHECK( (quad_quad_intersection(p1f, q1f, r1f, s1f, p6f, q6f, r6f, s6f)) );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p1i {-5,-5 };
  const Vec2i q1i { 5,-5 };
  const Vec2i r1i { 5, 5 };
  const Vec2i s1i {-5, 5 };

  // Regular intersecion
  const Vec2i p2i { 0,-10 };
  const Vec2i q2i { 10,-10 };
  const Vec2i r2i { 10, 0 };
  const Vec2i s2i { 0, 0 };

  // No intersection
  const Vec2i p3i {-25, 5 };
  const Vec2i q3i {-20, 0 };
  const Vec2i r3i {-15, 0 };
  const Vec2i s3i {-15, 10 };

  // Adjacent edge --> no intersection
  const Vec2i p4i {-5, 5 };
  const Vec2i q4i { 5, 5 };
  const Vec2i r4i { 5, 15 };
  const Vec2i s4i {-5, 15 };

  // Adjacent vertex --> no intersection
  const Vec2i p5i { 5, 5 };
  const Vec2i q5i { 15, 5 };
  const Vec2i r5i { 15, 15 };
  const Vec2i s5i { 5, 15 };

  // Half-Edge intersection
  const Vec2i p6i { 5,-5 };
  const Vec2i q6i { 10,-5 };
  const Vec2i r6i { 10, 0 };
  const Vec2i s6i { 5, 0 };

  CHECK( (quad_quad_intersection(p1i, q1i, r1i, s1i, p2i, q2i, r2i, s2i)) );
  CHECK(!(quad_quad_intersection(p1i, q1i, r1i, s1i, p3i, q3i, r3i, s3i)) );
  CHECK(!(quad_quad_intersection(p1i, q1i, r1i, s1i, p4i, q4i, r4i, s4i)) );
  CHECK(!(quad_quad_intersection(p1i, q1i, r1i, s1i, p5i, q5i, r5i, s5i)) );
  CHECK( (quad_quad_intersection(p1i, q1i, r1i, s1i, p6i, q6i, r6i, s6i)) );


} // quad_quad_intersection()

/*********************************************************************
* Test tri_quad_intersection
*********************************************************************/
void tri_quad_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d pd {-1.0,-0.5 };
  const Vec2d qd { 0.5,-0.5 };
  const Vec2d rd { 0.5, 1.5 };
  const Vec2d sd {-1.0, 1.5 };

  // Regular intersection
  const Vec2d m1d {-2.0, 0.5 };
  const Vec2d n1d {-1.5,-1.5 };
  const Vec2d o1d {-0.5, 0.0 };

  // No intersection
  const Vec2d m2d { 1.5,-1.0 };
  const Vec2d n2d { 2.5,-1.0 };
  const Vec2d o2d { 2.0,-0.5 };

  // Adjacent edges --> no intersection
  const Vec2d m3d { 0.5,-0.5 };
  const Vec2d n3d { 1.5, 0.5 };
  const Vec2d o3d { 0.5,-1.0 };

  // Adjacent vertex --> no intersection
  const Vec2d m4d { 0.5, 1.5 };
  const Vec2d n4d { 1.5, 1.0 };
  const Vec2d o4d { 1.5, 1.5 };

  // Half edge intersection
  const Vec2d m5d { 0.5,-0.5 };
  const Vec2d n5d { 1.5,-0.5 };
  const Vec2d o5d { 0.5, 0.5 };

  CHECK( (tri_quad_intersection(m1d, n1d, o1d, pd, qd, rd, sd)) );
  CHECK(!(tri_quad_intersection(m2d, n2d, o2d, pd, qd, rd, sd)) );
  CHECK(!(tri_quad_intersection(m3d, n3d, o3d, pd, qd, rd, sd)) );
  CHECK(!(tri_quad_intersection(m4d, n4d, o4d, pd, qd, rd, sd)) );
  CHECK( (tri_quad_intersection(m5d, n5d, o5d, pd, qd, rd, sd)) );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f pf {-1.0f,-0.5f };
  const Vec2f qf { 0.5f,-0.5f };
  const Vec2f rf { 0.5f, 1.5f };
  const Vec2f sf {-1.0f, 1.5f };

  // Regular intersection
  const Vec2f m1f {-2.0f, 0.5f };
  const Vec2f n1f {-1.5f,-1.5f };
  const Vec2f o1f {-0.5f, 0.0f };

  // No intersection
  const Vec2f m2f { 1.5f,-1.0f };
  const Vec2f n2f { 2.5f,-1.0f };
  const Vec2f o2f { 2.0f,-0.5f };

  // Adjacent edges --> no intersection
  const Vec2f m3f { 0.5f,-0.5f };
  const Vec2f n3f { 1.5f, 0.5f };
  const Vec2f o3f { 0.5f,-1.0f };

  // Adjacent vertex --> no intersection
  const Vec2f m4f { 0.5f, 1.5f };
  const Vec2f n4f { 1.5f, 1.0f };
  const Vec2f o4f { 1.5f, 1.5f };

  // Half edge intersection
  const Vec2f m5f { 0.5f,-0.5f };
  const Vec2f n5f { 1.5f,-0.5f };
  const Vec2f o5f { 0.5f, 0.5f };

  CHECK( (tri_quad_intersection(m1f, n1f, o1f, pf, qf, rf, sf)) );
  CHECK(!(tri_quad_intersection(m2f, n2f, o2f, pf, qf, rf, sf)) );
  CHECK(!(tri_quad_intersection(m3f, n3f, o3f, pf, qf, rf, sf)) );
  CHECK(!(tri_quad_intersection(m4f, n4f, o4f, pf, qf, rf, sf)) );
  CHECK( (tri_quad_intersection(m5f, n5f, o5f, pf, qf, rf, sf)) );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i pi {-10,-5 };
  const Vec2i qi { 5,-5 };
  const Vec2i ri { 5, 15 };
  const Vec2i si {-10, 15 };

  // Regular intersection
  const Vec2i m1i {-20, 5 };
  const Vec2i n1i {-15,-15 };
  const Vec2i o1i {-5, 0 };

  // No intersection
  const Vec2i m2i { 15,-10 };
  const Vec2i n2i { 25,-10 };
  const Vec2i o2i { 20,-5 };

  // Adjacent edges --> no intersection
  const Vec2i m3i { 5,-5 };
  const Vec2i n3i { 15, 5 };
  const Vec2i o3i { 5,-10 };

  // Adjacent vertex --> no intersection
  const Vec2i m4i { 5, 15 };
  const Vec2i n4i { 15, 10 };
  const Vec2i o4i { 15, 15 };

  // Half edge intersection
  const Vec2i m5i { 5,-5 };
  const Vec2i n5i { 15,-5 };
  const Vec2i o5i { 5, 5 };

  CHECK( (tri_quad_intersection(m1i, n1i, o1i, pi, qi, ri, si)) );
  CHECK(!(tri_quad_intersection(m2i, n2i, o2i, pi, qi, ri, si)) );
  CHECK(!(tri_quad_intersection(m3i, n3i, o3i, pi, qi, ri, si)) );
  CHECK(!(tri_quad_intersection(m4i, n4i, o4i, pi, qi, ri, si)) );
  CHECK( (tri_quad_intersection(m5i, n5i, o5i, pi, qi, ri, si)) );

} // tri_quad_intersection()


/*********************************************************************
* Test rect_overlap()
*********************************************************************/
void rect_overlap()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d al_d { -0.2, -0.2 };
  const Vec2d au_d {  0.2,  0.2 };

  const Vec2d bl_d {  0.0,  0.0 };
  const Vec2d bu_d {  0.2,  0.2 };

  const Vec2d cl_d {  0.1,  0.1 };
  const Vec2d cu_d {  0.1,  0.4 };

  const Vec2d dl_d {  0.2, -0.2 };
  const Vec2d du_d {  0.4,  0.3 };

  const Vec2d el_d {  0.2,  0.2 };
  const Vec2d eu_d {  0.4,  0.4 };

  const Vec2d fl_d {  0.3,  0.3 };
  const Vec2d fu_d {  0.4,  0.4 };

  CHECK( (rect_overlap(al_d, au_d, bl_d, bu_d)) );
  CHECK( (rect_overlap(al_d, au_d, cl_d, cu_d)) );
  CHECK( (rect_overlap(al_d, au_d, dl_d, du_d)) );
  CHECK( (rect_overlap(al_d, au_d, el_d, eu_d)) );
  CHECK(!(rect_overlap(al_d, au_d, fl_d, fu_d)) );

} // rect_overlap()


/*********************************************************************
* Test in_on_rect()
*********************************************************************/
void in_on_rect()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d al_d {-0.2, -0.2};
  const Vec2d au_d { 0.2,  0.2};
  
  const Vec2d p_d { 0.0,  0.0};
  const Vec2d q_d { 0.2,  0.0};
  const Vec2d r_d { 0.2,  0.2};
  const Vec2d s_d { 0.3,  0.3};


  CHECK( (in_on_rect(p_d, al_d, au_d)) );
  CHECK( (in_on_rect(q_d, al_d, au_d)) );
  CHECK( (in_on_rect(r_d, al_d, au_d)) );
  CHECK(!(in_on_rect(s_d, al_d, au_d)) );


} // in_on_rect()


/*********************************************************************
* Test in_rect()
*********************************************************************/
void in_rect()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d al_d {-0.2, -0.2};
  const Vec2d au_d { 0.2,  0.2};
  
  const Vec2d p_d { 0.0,  0.0};
  const Vec2d q_d { 0.2,  0.0};
  const Vec2d r_d { 0.2,  0.2};
  const Vec2d s_d { 0.3,  0.3};


  CHECK( (in_rect(p_d, al_d, au_d)) );
  CHECK(!(in_rect(q_d, al_d, au_d)) );
  CHECK(!(in_rect(r_d, al_d, au_d)) );
  CHECK(!(in_rect(s_d, al_d, au_d)) );


} // in_rect()

/*********************************************************************
* Test in_on_triangle()
*********************************************************************/
void in_on_triangle()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d t1_d {-0.2, -0.2};
  const Vec2d t2_d { 0.2, -0.2};
  const Vec2d t3_d { 0.0,  0.2};

  const Vec2d a_d { 0.0 , 0.0};
  const Vec2d b_d { 0.0 ,-0.2};
  const Vec2d c_d {-0.2 ,-0.2};
  const Vec2d d_d {-0.2 , 0.2};

  CHECK( (in_on_triangle(a_d, t1_d, t2_d, t3_d)) );
  CHECK( (in_on_triangle(b_d, t1_d, t2_d, t3_d)) );
  CHECK( (in_on_triangle(c_d, t1_d, t2_d, t3_d)) );
  CHECK(!(in_on_triangle(d_d, t1_d, t2_d, t3_d)) );


} // in_on_triangle()

/*********************************************************************
* Test in_triangle()
*********************************************************************/
void in_triangle()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d t1_d {-0.2, -0.2};
  const Vec2d t2_d { 0.2, -0.2};
  const Vec2d t3_d { 0.0,  0.2};

  const Vec2d a_d { 0.0 , 0.0};
  const Vec2d b_d { 0.0 ,-0.2};
  const Vec2d c_d {-0.2 ,-0.2};
  const Vec2d d_d {-0.2 , 0.2};

  CHECK( (in_triangle(a_d, t1_d, t2_d, t3_d)) );
  CHECK(!(in_triangle(b_d, t1_d, t2_d, t3_d)) );
  CHECK(!(in_triangle(c_d, t1_d, t2_d, t3_d)) );
  CHECK(!(in_triangle(d_d, t1_d, t2_d, t3_d)) );

} // in_triangle()

/*********************************************************************
* Test in_on_quad()
*********************************************************************/
void in_on_quad()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d q1_d { -0.2, -0.2 };
  const Vec2d q2_d {  0.2, -0.2 };
  const Vec2d q3_d {  0.2,  0.2 };
  const Vec2d q4_d { -0.2,  0.2 };

  const Vec2d a_d { 0.0, 0.0 };
  const Vec2d b_d {-0.2, 0.0 };
  const Vec2d c_d {-0.2, 0.0 };
  const Vec2d d_d {-0.3, 0.3 };

  CHECK( (in_on_quad(a_d, q1_d, q2_d, q3_d, q4_d)) );
  CHECK( (in_on_quad(b_d, q1_d, q2_d, q3_d, q4_d)) );
  CHECK( (in_on_quad(c_d, q1_d, q2_d, q3_d, q4_d)) );
  CHECK(!(in_on_quad(d_d, q1_d, q2_d, q3_d, q4_d)) );

} // in_on_quad()

/*********************************************************************
* Test in_quad()
*********************************************************************/
void in_quad()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d q1_d { -0.2, -0.2 };
  const Vec2d q2_d {  0.2, -0.2 };
  const Vec2d q3_d {  0.2,  0.2 };
  const Vec2d q4_d { -0.2,  0.2 };

  const Vec2d a_d { 0.0, 0.0 };
  const Vec2d b_d {-0.2, 0.0 };
  const Vec2d c_d {-0.2, 0.0 };
  const Vec2d d_d {-0.3, 0.3 };

  CHECK( (in_quad(a_d, q1_d, q2_d, q3_d, q4_d)) );
  CHECK(!(in_quad(b_d, q1_d, q2_d, q3_d, q4_d)) );
  CHECK(!(in_quad(c_d, q1_d, q2_d, q3_d, q4_d)) );
  CHECK(!(in_quad(d_d, q1_d, q2_d, q3_d, q4_d)) );

} // in_quad()

} // namespace GeometryTests


/*********************************************************************
* Run tests for: Geometry.h
*********************************************************************/
void run_tests_Geometry()
{
  GeometryTests::orientation();
  GeometryTests::is_left();
  GeometryTests::is_lefton();
  GeometryTests::in_segment();
  GeometryTests::in_on_segment();
  GeometryTests::line_line_intersection();
  GeometryTests::line_tri_intersection();
  GeometryTests::line_quad_intersection();
  GeometryTests::tri_tri_intersection();
  GeometryTests::quad_quad_intersection();
  GeometryTests::tri_quad_intersection();
  GeometryTests::rect_overlap();
  GeometryTests::in_on_rect();
  GeometryTests::in_rect();
  GeometryTests::in_on_triangle();
  GeometryTests::in_triangle();
  GeometryTests::in_on_quad();
  GeometryTests::in_quad();

} // run_tests_Geometry()
