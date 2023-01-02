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
#include "QuadTree.h"

namespace QuadTreeTests 
{
using namespace CppUtils;

/*--------------------------------------------------------------------
| A simple vertex class that is stored in the QuadTree structure
--------------------------------------------------------------------*/
template<typename T>
class VertexType
{
public:
  VertexType(T x, T y) : xy_ {x, y} {}

  const Vec2<T>& xy() const { return xy_; }

private:
  Vec2<T> xy_ {};
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
  }

  const Vec2<T>& xy() const { return xy_; }

  VertexType<T>& v1() { return *v1_; };
  const VertexType<T>& v1() const { return *v1_; };

  VertexType<T>& v2() { return *v2_; };
  const VertexType<T>& v2() const { return *v2_; };

private:
  VertexType<T>* v1_ { nullptr };
  VertexType<T>* v2_ { nullptr };
  Vec2<T>        xy_ { 0.0, 0.0 };
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

  const V d_sqr = distance_point_edge_sqr(query, xy_1, xy_2);

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
using std::list;
using std::make_unique;


/*--------------------------------------------------------------------
| Test constructor
--------------------------------------------------------------------*/
void constructor()
{
  const Vec2d center { 0.0, 0.0 };
  double scale       { 1.0 };
  size_t max_item    { 10 };
  size_t max_depth   { 3 };

  QuadTree<Vertex,double> 
    quadtree { scale, max_item, max_depth, center };

  CHECK( ( EQ(quadtree.lowleft().x, -0.5) ) );
  CHECK( ( EQ(quadtree.lowleft().y, -0.5) ) );

  CHECK( ( EQ(quadtree.upright().x,  0.5) ) );
  CHECK( ( EQ(quadtree.upright().y,  0.5) ) );

  CHECK( ( quadtree.size() == 0 ) );

  CHECK(!( quadtree.children()[0] ) );
  CHECK(!( quadtree.children()[1] ) );
  CHECK(!( quadtree.children()[2] ) );
  CHECK(!( quadtree.children()[3] ) );

} // constructor()

/*--------------------------------------------------------------------
| Test add
--------------------------------------------------------------------*/
void add()
{
  /*------------------------------------------------------------------
  | Initialize structure and vertex container
  ------------------------------------------------------------------*/
  const Vec2d center { 0.0, 0.0 };
  double scale       { 1.0 };
  size_t max_item    { 3 };
  size_t max_depth   { 5 };

  QuadTree<Vertex,double> 
    quadtree { scale, max_item, max_depth, center };

  // Vertex-container
  vector<unique_ptr<Vertex>> vertices;

  vertices.push_back( make_unique<Vertex>( 0.1, 0.1 ) );
  vertices.push_back( make_unique<Vertex>( 0.2, 0.2 ) );
  vertices.push_back( make_unique<Vertex>( 0.4, 0.2 ) );
  vertices.push_back( make_unique<Vertex>( 0.4, 0.3 ) );

  /*------------------------------------------------------------------
  | Add vertices
  ------------------------------------------------------------------*/
  quadtree.add( vertices[0].get() );
  quadtree.add( vertices[1].get() );
  quadtree.add( vertices[2].get() );

  CHECK(!( quadtree.split() ) );

  /*------------------------------------------------------------------
  | quadtree should split 
  ------------------------------------------------------------------*/
  quadtree.add( vertices[3].get() );

  CHECK( ( quadtree.split() ) );
  CHECK( ( quadtree.size() == 4 ) );

  CHECK( ( quadtree.children()[0]->split() ) );

  CHECK( ( quadtree.children()[0]->children()[0]->size() == 1 ) );
  CHECK( ( quadtree.children()[0]->children()[1]->size() == 0 ) );
  CHECK( ( quadtree.children()[0]->children()[2]->size() == 2 ) );
  CHECK( ( quadtree.children()[0]->children()[3]->size() == 1 ) );

  CHECK( ( quadtree.n_leafs() == 7 ) );

  /*------------------------------------------------------------------
  | Check splitting behaviour for collision 
  ------------------------------------------------------------------*/
  vertices.push_back( make_unique<Vertex>( 0.5,-0.5 ) );
  vertices.push_back( make_unique<Vertex>( 0.5,-0.5 ) );
  vertices.push_back( make_unique<Vertex>( 0.5,-0.5 ) );
  vertices.push_back( make_unique<Vertex>( 0.5,-0.5 ) );
  vertices.push_back( make_unique<Vertex>( 0.5,-0.5 ) );

  quadtree.add( vertices[4].get() );
  quadtree.add( vertices[5].get() );
  quadtree.add( vertices[6].get() );
  quadtree.add( vertices[7].get() );
  quadtree.add( vertices[8].get() );

  // Container in lower right corner should contain more elements
  // than given as maximum number
  CHECK( ( quadtree.children()[3]->children()[3]->children()[3]->
            children()[3]->children()[3]->size() == 5 ) );

  CHECK( ( quadtree.size() == 9 ) );

} // add()

/*--------------------------------------------------------------------
| Test remove
--------------------------------------------------------------------*/
void remove()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  const Vec2d center { 0.0, 0.0 };
  double scale       { 1.0 };
  size_t max_item    { 3 };
  size_t max_depth   { 5 };

  QuadTree<Vertex,double> 
    quadtree { scale, max_item, max_depth, center };

  // Vertex-container
  vector<unique_ptr<Vertex>> vertices;

  vertices.push_back( make_unique<Vertex>( 0.1, 0.1 ) );
  quadtree.add( vertices[0].get() );

  vertices.push_back( make_unique<Vertex>( 0.2, 0.2 ) );
  quadtree.add( vertices[1].get() );

  vertices.push_back( make_unique<Vertex>( 0.4, 0.2 ) );
  quadtree.add( vertices[2].get() );

  vertices.push_back( make_unique<Vertex>( 0.4, 0.3 ) );
  quadtree.add( vertices[3].get() );

  /*------------------------------------------------------------------
  | Remove / add / remove vertices
  ------------------------------------------------------------------*/
  quadtree.remove( vertices[2].get() );
  quadtree.add( vertices[2].get() );
  quadtree.remove( vertices[2].get() );

  CHECK( ( quadtree.size() == 3 ) );
  CHECK(!( quadtree.split() ) );
  CHECK( ( quadtree.n_leafs() == 1 ) );

  quadtree.remove( vertices[1].get() );
  quadtree.add( vertices[1].get() );
  quadtree.remove( vertices[1].get() );

  quadtree.remove( vertices[3].get() );
  quadtree.add( vertices[3].get() );
  quadtree.remove( vertices[3].get() );

  quadtree.remove( vertices[0].get() );
  quadtree.add( vertices[0].get() );
  quadtree.remove( vertices[0].get() );

  CHECK( ( quadtree.size() == 0 ) );

} // remove()

/*--------------------------------------------------------------------
| Test get_items()
--------------------------------------------------------------------*/
void get_items()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  const Vec2d center { 0.0, 0.0 };
  double scale       { 1.0 };
  size_t max_item    { 3 };
  size_t max_depth   { 5 };

  QuadTree<Vertex,double> 
    quadtree { scale, max_item, max_depth, center };

  // Vertex-container
  vector<unique_ptr<Vertex>> vertices;

  vertices.push_back( make_unique<Vertex>( 0.1, 0.1 ) );
  quadtree.add( vertices[0].get() );

  vertices.push_back( make_unique<Vertex>( 0.2, 0.2 ) );
  quadtree.add( vertices[1].get() );

  vertices.push_back( make_unique<Vertex>( 0.4, 0.2 ) );
  quadtree.add( vertices[2].get() );

  vertices.push_back( make_unique<Vertex>( 0.4, 0.3 ) );
  quadtree.add( vertices[3].get() );

  vertices.push_back( make_unique<Vertex>(-0.2,-0.2 ) );
  quadtree.add( vertices[4].get() );

  vertices.push_back( make_unique<Vertex>(-0.4,-0.4 ) );
  quadtree.add( vertices[5].get() );

  vertices.push_back( make_unique<Vertex>( 0.0,-0.3 ) );
  quadtree.add( vertices[6].get() );

  vertices.push_back( make_unique<Vertex>(-0.1,-0.2 ) );
  quadtree.add( vertices[7].get() );

  vertices.push_back( make_unique<Vertex>( 0.1,-0.1 ) );
  quadtree.add( vertices[8].get() );


  /*------------------------------------------------------------------
  | Get items in rectangle
  ------------------------------------------------------------------*/
  vector<Vertex*> r_found {};
  
  size_t n_r = quadtree.get_items( {-0.25,-0.25}, {0.25,0.25}, r_found );

  CHECK( ( n_r == 5 ) );
  CHECK( ( n_r == r_found.size() ) );

  /*------------------------------------------------------------------
  | Get items in circle
  ------------------------------------------------------------------*/
  vector<Vertex*> c_found {};
  
  size_t n_c = quadtree.get_items( {0.25,0.25}, 0.25, c_found );

  CHECK( ( n_c == 4 ) );
  CHECK( ( n_c == c_found.size() ) );

} // get_items()

/*--------------------------------------------------------------------
| Test add() / remove()
--------------------------------------------------------------------*/
void add_remove()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  const Vec2d center { 0.0, 0.0 };
  double scale       { 2.0 };
  size_t max_item    { 3 };
  size_t max_depth   { 5 };

  QuadTree<Vertex,double> 
    quadtree { scale, max_item, max_depth, center };

  // Vertex-container
  vector<unique_ptr<Vertex>> vertices;

  vertices.push_back( make_unique<Vertex>( 0.1, 0.1 ) );
  quadtree.add( vertices[0].get() );

  vertices.push_back( make_unique<Vertex>( 0.2, 0.2 ) );
  quadtree.add( vertices[1].get() );

  vertices.push_back( make_unique<Vertex>( 0.3, 0.3 ) );
  quadtree.add( vertices[2].get() );

  vertices.push_back( make_unique<Vertex>( 0.3, 0.3 ) );
  quadtree.add( vertices[3].get() );

  CHECK( ( quadtree.split() ) );

  quadtree.remove( vertices[3].get() );
  CHECK( !( quadtree.split() ) );

  quadtree.remove( vertices[2].get() );
  quadtree.remove( vertices[1].get() );

  quadtree.add( vertices[3].get() );
  quadtree.add( vertices[2].get() );
  quadtree.add( vertices[1].get() );
  CHECK( ( quadtree.split() ) );

} // add_remove()


/*--------------------------------------------------------------------
| Test get_nearest()
--------------------------------------------------------------------*/
void get_nearest()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  const Vec2d center { 2.0, 2.0 };
  double scale       { 4.0 };
  size_t max_item    { 2 };
  size_t max_depth   { 5 };

  QuadTree<Vertex,double> 
    quadtree { scale, max_item, max_depth, center };

  // Vertex-container
  vector<unique_ptr<Vertex>> vertices;

  vertices.push_back( make_unique<Vertex>( 1.0, 1.0 ) );
  quadtree.add( vertices[0].get() );

  CHECK( (vertices[0].get() == quadtree.get_nearest( {2.1, 2.2} ) ) );

  vertices.push_back( make_unique<Vertex>( 2.5, 0.5 ) );
  quadtree.add( vertices[1].get() );

  CHECK( ( !quadtree.split() ) );

  vertices.push_back( make_unique<Vertex>( 3.0, 3.0 ) );
  quadtree.add( vertices[2].get() );

  CHECK( ( quadtree.split() ) );

  vertices.push_back( make_unique<Vertex>( 3.5, 0.5 ) );
  quadtree.add( vertices[3].get() );

  vertices.push_back( make_unique<Vertex>( 3.5, 1.5 ) );
  quadtree.add( vertices[4].get() );

  CHECK( (vertices[4].get() == quadtree.get_nearest( {3.5, 2.25} ) ) );
  CHECK( (vertices[2].get() == quadtree.get_nearest( {1.0, 3.00} ) ) );

} // get_nearest() 


/*--------------------------------------------------------------------
| Test get_leaf()
--------------------------------------------------------------------*/
void get_leaf()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  const Vec2d center { 2.0, 2.0 };
  double scale       { 4.0 };
  size_t max_item    { 2 };
  size_t max_depth   { 5 };

  QuadTree<Vertex,double> 
    quadtree { scale, max_item, max_depth, center };

  // Vertex-container
  vector<unique_ptr<Vertex>> vertices;

  vertices.push_back( make_unique<Vertex>( 1.0, 1.0 ) );
  vertices.push_back( make_unique<Vertex>( 2.5, 0.5 ) );
  vertices.push_back( make_unique<Vertex>( 3.0, 3.0 ) );
  vertices.push_back( make_unique<Vertex>( 3.5, 0.5 ) );
  vertices.push_back( make_unique<Vertex>( 3.5, 1.5 ) );

  quadtree.add( vertices[0].get() );
  quadtree.add( vertices[1].get() );
  quadtree.add( vertices[2].get() );
  quadtree.add( vertices[3].get() );
  quadtree.add( vertices[4].get() );


  const QuadTree<Vertex,double>* leaf = nullptr;

  leaf = quadtree.get_leaf( {3.4, 1.4} );
  CHECK( leaf );
  CHECK( (!leaf->split()) );
  CHECK( (leaf->size() == 1) );
  CHECK( (leaf->items().front() == vertices[4].get()) );

  leaf = quadtree.get_leaf( {0.9, 3.4} );
  CHECK( leaf );
  CHECK( (!leaf->split()) );  
  CHECK( (leaf->size() == 0) );

  leaf = quadtree.get_leaf( {2.2, 2.2} );
  CHECK( leaf );
  CHECK( (!leaf->split()) );  
  CHECK( (leaf->size() == 1) );
  CHECK( (leaf->items().front() == vertices[2].get()) );


} // get_leaf()


/*--------------------------------------------------------------------
| Test get_neaerst() for edges
--------------------------------------------------------------------*/
void get_nearest_edge()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  const Vec2d center { 4.0, 4.0 };
  double scale       { 8.0 };
  size_t max_item    { 2 };
  size_t max_depth   { 5 };

  QuadTree<Edge,double> 
    quadtree { scale, max_item, max_depth, center };

  /*------------------------------------------------------------------
  | Initialize vertices
  ------------------------------------------------------------------*/
  vector<unique_ptr<Vertex>> vertices;
  vertices.push_back( make_unique<Vertex>( 0.5, 0.5 ) ); // 0
  vertices.push_back( make_unique<Vertex>( 1.5, 1.5 ) ); // 1
  vertices.push_back( make_unique<Vertex>( 7.0, 5.5 ) ); // 2
  vertices.push_back( make_unique<Vertex>( 7.0, 6.5 ) ); // 3
  vertices.push_back( make_unique<Vertex>( 6.0, 2.0 ) ); // 4
  vertices.push_back( make_unique<Vertex>( 7.0, 1.0 ) ); // 5
  vertices.push_back( make_unique<Vertex>( 4.5, 4.5 ) ); // 6
  vertices.push_back( make_unique<Vertex>( 5.5, 5.0 ) ); // 7
  vertices.push_back( make_unique<Vertex>( 7.5, 4.5 ) ); // 8

  /*------------------------------------------------------------------
  | Initialize edges
  ------------------------------------------------------------------*/
  vector<unique_ptr<Edge>> edges;
  edges.push_back( make_unique<Edge>( *vertices[0], *vertices[1] ) ); // 0
  edges.push_back( make_unique<Edge>( *vertices[2], *vertices[3] ) ); // 1
  edges.push_back( make_unique<Edge>( *vertices[4], *vertices[5] ) ); // 2

  /*------------------------------------------------------------------
  | Add edges to quadtree
  ------------------------------------------------------------------*/
  quadtree.add( edges[0].get() ); 
  quadtree.add( edges[1].get() ); 

  CHECK( !quadtree.split() );
  quadtree.add( edges[2].get() ); 
  CHECK(  quadtree.split() );


  Edge* e_nearest = nullptr;

  e_nearest = quadtree.get_nearest( {2.0,1.0}, nearest_edge_fun );
  CHECK( e_nearest == edges[0].get() );

  e_nearest = quadtree.get_nearest( {1.0,4.5}, nearest_edge_fun );
  CHECK( e_nearest == edges[0].get() );


  /*------------------------------------------------------------------
  | More edges
  ------------------------------------------------------------------*/
  edges.push_back( make_unique<Edge>( *vertices[6], *vertices[7] ) ); // 3
  edges.push_back( make_unique<Edge>( *vertices[7], *vertices[8] ) ); // 4

  quadtree.add( edges[3].get() ); 
  quadtree.add( edges[4].get() ); 

  e_nearest = quadtree.get_nearest( {6.0,5.5}, nearest_edge_fun );
  CHECK( e_nearest == edges[4].get() );


} // get_nearest_edge()

} // namespace QuadTreeTests


/*********************************************************************
* Run tests for: QuadTree.h
*********************************************************************/
void run_tests_QuadTree()
{
  QuadTreeTests::constructor();
  QuadTreeTests::add();
  QuadTreeTests::remove();
  QuadTreeTests::get_items();
  QuadTreeTests::add_remove();
  QuadTreeTests::get_nearest();
  QuadTreeTests::get_leaf();
  QuadTreeTests::get_nearest_edge();

} // run_tests_QuadTree()
