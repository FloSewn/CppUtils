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

#include "tests.h"

#include "VecND.h"
#include "Testing.h"
#include "RTreeND.h"

namespace RTreeNDTests 
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
using BBox2d = BBoxND<double, 2>;


/*--------------------------------------------------------------------
| Rectangular class for testing of the RTreeND
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
class TestObject
{
  using Vec  = VecND<T,N>;
  using BBox = BBoxND<T,N>;

public:
  //TestObject() {}

  TestObject(const Vec& ll, const Vec& ur)
  : bbox_ { ll, ur }
  {
    xy_ = 0.5 * (ll + ur);
  }

  Vec& xy() { return xy_; }
  const Vec& xy() const { return xy_; }

  BBox& bbox() { return bbox_; }
  const BBox& bbox() const { return bbox_; }

private:
  BBox bbox_ {};
  Vec  xy_   {};

}; // TestObject

using TestRect = TestObject<double, 2>;

/*--------------------------------------------------------------------
| Test constructor
--------------------------------------------------------------------*/
void constructor()
{
  RTreeND<TestRect,3,double,2> tree {};


  std::vector<TestRect> rectangles;

  rectangles.push_back( { {0.0,0.0}, {1.0,1.0} } );
  rectangles.push_back( { {0.5,1.0}, {1.5,2.0} } );
  rectangles.push_back( { {2.0,0.5}, {3.0,1.5} } );
  rectangles.push_back( { {2.0,2.0}, {3.0,3.0} } );
  rectangles.push_back( { {4.0,4.0}, {5.0,5.0} } );
  rectangles.push_back( { {4.0,1.0}, {4.5,1.5} } );
  rectangles.push_back( { {3.0,4.5}, {5.0,5.5} } );

  for (TestRect& r : rectangles)
    tree.insert( r );

  //CHECK( tree.root().n_entries() == 2 );


  tree.print(std::cout);

  /*

  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/test_data/RTree_constructor.txt" };

  tree.write_to_file( file_name );
  */


} // constructor()

/*--------------------------------------------------------------------
| Test bulk insertion of many objects
--------------------------------------------------------------------*/
void bulk_insertion()
{
  /*
  RTreeND<TestRect, 3> tree {};

  std::vector<TestRect> rectangles;

  rectangles.push_back( { {3.0,4.5}, {5.0,5.5} } );
  rectangles.push_back( { {0.0,0.0}, {1.0,1.0} } );
  rectangles.push_back( { {4.0,4.0}, {5.0,5.0} } );
  rectangles.push_back( { {2.0,0.5}, {3.0,1.5} } );
  rectangles.push_back( { {0.5,1.0}, {1.5,2.0} } );
  rectangles.push_back( { {2.0,2.0}, {3.0,3.0} } );
  rectangles.push_back( { {4.0,1.0}, {4.5,1.5} } );

  tree.insert( rectangles );


  tree.print(std::cout);


  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };
  std::string file_name 
  { source_dir + "/auxiliary/test_data/RTree_bulk_insertion.txt" };

  tree.write_to_file( file_name );
  */


} // bulk_insertion()

} // namespace RTreeNDTests


/*********************************************************************
* Run tests for: RTreeND.h
*********************************************************************/
void run_tests_RTreeND()
{
  RTreeNDTests::constructor();

  RTreeNDTests::bulk_insertion();

} // run_tests_RTreeND()
