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

#include "tests.h"

#include "Vec2.h"
#include "Testing.h"
#include "RTree.h"

namespace RTreeTests 
{
using namespace CppUtils;

/*--------------------------------------------------------------------
| A simple vertex class that is stored in the QuadTree structure
--------------------------------------------------------------------*/
template<typename T>
class VertexType
{
public:
  VertexType(T x, T y) 
  : xy_ {x, y} 
  {}

  const Vec2<T>& xy() const { return xy_; }
  const Vec2<T>& lowleft() const { return xy_; }
  const Vec2<T>& upright() const { return xy_; }

private:
  Vec2<T> xy_;
}; 

/*--------------------------------------------------------------------
| A simple edge class that is stored in the QuadTree structure
--------------------------------------------------------------------*/
template<typename T>
class EdgeType
{
public:
  EdgeType(VertexType<T>& v1, VertexType<T>& v2)
  : v1_ { &v1 }
  , v2_ { &v2 }
  {
    xy_ = ( v1.xy() + v2.xy() ) / 2.0;

    double x_min = MIN( v1.xy().x, v2.xy().x );
    double x_max = MAX( v1.xy().x, v2.xy().x );

    double y_min = MIN( v1.xy().y, v2.xy().y );
    double y_max = MAX( v1.xy().y, v2.xy().y );

    lowleft_ = { x_min, y_min };
    upright_ = { x_max, y_max };
  }

  const Vec2<T>& xy() const { return xy_; }
  const Vec2<T>& lowleft() const { return lowleft_; }
  const Vec2<T>& upright() const { return upright_; }

  VertexType<T>& v1() { return *v1_; };
  const VertexType<T>& v1() const { return *v1_; };

  VertexType<T>& v2() { return *v2_; };
  const VertexType<T>& v2() const { return *v2_; };

private:
  VertexType<T>* v1_ { nullptr };
  VertexType<T>* v2_ { nullptr };
  Vec2<T>        xy_ { 0, 0 };
  Vec2<T>        lowleft_ { 0, 0 };
  Vec2<T>        upright_ { 0, 0 };
}; 


/*--------------------------------------------------------------------
| Query function for nearest point to an edge
--------------------------------------------------------------------*/
template <typename T, typename V>
static inline bool nearest_edge_fun(T* edge,
                                    const Vec2<V>& query,
                                    V& min_dist_sqr)
{
  const VertexType<V>& v1 = edge->v1();
  const VertexType<V>& v2 = edge->v2();

  const Vec2<V>& xy_1 = v1.xy();
  const Vec2<V>& xy_2 = v2.xy();

  const V d_sqr = vertex_edge_dist_sqr(query, xy_1, xy_2);

  if (d_sqr < min_dist_sqr)
  {
    min_dist_sqr = d_sqr;
    return true;
  }

  return false;
}


/*--------------------------------------------------------------------
| 
--------------------------------------------------------------------*/
using Vertex = VertexType<double>;
using Edge   = EdgeType<double>;
using std::unique_ptr;
using std::vector;
using std::make_unique;


/*--------------------------------------------------------------------
| Test constructor
--------------------------------------------------------------------*/
void constructor()
{

  Vec2d lowleft = {  0.0,  0.0 };
  Vec2d upright = { 10.0, 10.0 };

  RTree<Vertex,2> rtree { lowleft, upright };


} // constructor()

} // namespace RTreeTests


/*********************************************************************
* Run tests for: RTree.h
*********************************************************************/
void run_tests_RTree()
{
  RTreeTests::constructor();

} // run_tests_RTree()
