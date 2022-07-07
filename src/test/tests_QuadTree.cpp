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

  const Vec2<T> xy() const { return xy_; }

private:
  Vec2<T>                   xy_;
}; 

using Vertex = VertexType<double>;
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

} // run_tests_QuadTree()
