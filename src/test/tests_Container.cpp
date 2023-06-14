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
#include "Container.h"
#include "Testing.h"

#include <list>           // std::list
#include <memory>           // std::list

namespace ContainerTests 
{
using namespace CppUtils;

class Edge;

/*********************************************************************
* A simple vertex class that is stored in the Container
*********************************************************************/
class Vertex : public ContainerEntry<Vertex>
{
public:

  using EdgeList = std::list<Edge*>;

  /*------------------------------------------------------------------
  | Constructor 
  ------------------------------------------------------------------*/
  Vertex(double x, double y) : ContainerEntry<Vertex>(x,y) {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const EdgeList& edges() const { return edges_;}
  const Edge& edges(size_t i) const 
  { 
    ASSERT(!( i < 0 || i >= edges_.size() ),
            "Invalid access to vertex edge list." );
    auto iter = edges_.begin();
    std::advance( iter, i );
    ASSERT( *iter, "Edge vertex is nullptr.");
    return *(*iter);
  }

  void add_edge(Edge* e) { edges_.push_back(e); }
  void remove_edge(Edge* e) { edges_.remove(e); }

private:

  /*------------------------------------------------------------------
  | Vertex attributes 
  ------------------------------------------------------------------*/
  EdgeList edges_ {};

}; 


/*********************************************************************
* A simple edge class that is stored in the Container
*********************************************************************/
class Edge : public ContainerEntry<Edge>
{
public:
  //friend Container<Edge>;

  /*------------------------------------------------------------------
  | Constructor 
  ------------------------------------------------------------------*/
  Edge(Vertex& v1, Vertex& v2) 
  : ContainerEntry<Edge>(0.5 * (v1.xy() + v2.xy()))
  , v1_ {&v1}
  , v2_ {&v2}
  {
    ASSERT((&v1_ && &v2_),
        "Failed to create edge structure due to given nullptr." );
      
    const Vec2d d_xy = v2_->xy() - v1_->xy();

    xy_     = 0.5 * ( v1_->xy() + v2_->xy() );
    length_ = d_xy.norm();
    tang_   = d_xy / length_;

    norm_.x = -tang_.y;
    norm_.y =  tang_.x;

    v1_->add_edge( this );
    v2_->add_edge( this );
  }

  /*------------------------------------------------------------------
  | Destructor 
  ------------------------------------------------------------------*/
  ~Edge() {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const Vertex& v1() const { return *v1_; };
  const Vertex& v2() const { return *v2_; };
  Vertex& v1() { return *v1_; };
  Vertex& v2() { return *v2_; };

  double length() const { return length_; }
  const Vec2d& normal() const { return norm_;}
  const Vec2d& tangent() const { return tang_;}

  /*------------------------------------------------------------------
  | Container functions 
  ------------------------------------------------------------------*/
  void container_destructor() override
  { 
    if (v1_) v1_->remove_edge( this );
    if (v2_) v2_->remove_edge( this );

    v1_ = nullptr;
    v2_ = nullptr;
  }


private:
  /*------------------------------------------------------------------
  | Edge attributes 
  ------------------------------------------------------------------*/
  Vertex*     v1_;
  Vertex*     v2_;

  Vec2d       xy_;
  Vec2d       tang_;
  Vec2d       norm_;
  double      length_;
}; 

/*********************************************************************
* Vertex ostream operator overload
*********************************************************************/
static std::ostream& operator<<(std::ostream& os, 
                                const Vertex& v)
{ return os << v.xy(); }

/*********************************************************************
* Edge ostream operator overload
*********************************************************************/
static std::ostream& operator<<(std::ostream& os, 
                                const Edge& e)
{ return os << e.v1() << " --> "  << e.v2(); }

/***********************************************************
* Vertex equality operator 
***********************************************************/
static bool operator==(const Vertex& v1, const Vertex& v2)
{ return &v1 == &v2; }
static bool operator!=(const Vertex& v1, const Vertex& v2)
{ return !(v1 == v2); }

/***********************************************************
* Edge equality operator 
***********************************************************/
static bool operator==(const Edge& e1, const Edge& e2)
{ return &e1 == &e2; }
static bool operator!=(const Edge& e1, const Edge& e2)
{ return !(e1 == e2); }


/*********************************************************************
* Test Container::push_back()
*********************************************************************/
void push_back()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  Container<Vertex> vertices {};

  vertices.push_back( 0.0, 0.0 );
  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  vertices.push_back( 2.0, 1.0 );
  vertices.push_back( 3.0, 0.0 );
  vertices.push_back( 2.0, 0.0 );

  CHECK(v1.in_container());

  CHECK( (vertices.size() == 5) );

  CHECK( (vertices[1] == v1) );

  // Check Container::back() function
  CHECK( (vertices.back() == vertices[4]) );

  // Check Container::front() function
  CHECK( (vertices.front() == vertices[0]) );


} // push_back()

/*********************************************************************
* Test Container::insert()
*********************************************************************/
void insert()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  Container<Vertex> vertices {};

  vertices.push_back( 0.0, 0.0 );

  Vertex& v3 = vertices.push_back( 3.0, 1.0 );

  // Insert before vertices[1]
  Vertex& v1 = vertices.insert( vertices[1].pos(), 1.0, 1.0 );
  Vertex& v2 = vertices.insert( v3.pos(), 2.0, 1.0 );

  vertices.push_back( 2.0, 0.0 );

  CHECK( (vertices.size() == 5) );

  CHECK( (vertices[1] == v1) );
  CHECK( (vertices[1].xy().x == 1.0 && vertices[1].xy().y == 1.0 ) );

  CHECK( (vertices[2] == v2) );
  CHECK( (vertices[2].xy().x == 2.0 && vertices[2].xy().y == 1.0 ) );

} // insert()

/*********************************************************************
* Test Container::remove()
*********************************************************************/
void remove()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  Container<Vertex> vertices {};

  vertices.push_back( 0.0, 0.0 );
  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  vertices.push_back( 2.0, 1.0 );
  vertices.push_back( 3.0, 0.0 );
  Vertex& v4 = vertices.push_back( 2.0, 0.0 );


  vertices.remove( vertices[3] );

  CHECK( (vertices.size() == 4) );

  CHECK( (vertices.waste().size() == 1) );

  vertices.remove( v1 );

  CHECK( (vertices.size() == 3) );

  CHECK( (vertices.waste().size() == 2) );

  vertices.remove( v4 );

  CHECK( (vertices.size() == 2) );

  vertices.remove( vertices[0] );
  vertices.remove( vertices[0] );

  CHECK( (vertices.size() == 0) );

  // Free everything in the garbage collector
  vertices.clear_waste();
  CHECK( (vertices.waste().size() == 0) );

} // remove()

/*********************************************************************
* Test Container of Vertices and Edges
*********************************************************************/
void vertices_and_edges()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  Container<Vertex> vertices {};
  Container<Edge>   edges {};

  vertices.push_back( 0.0, 0.0 );
  vertices.push_back( 1.0, 1.0 );
  vertices.push_back( 2.0, 1.0 );
  vertices.push_back( 3.0, 1.0 );
  vertices.push_back( 2.0, 0.0 );

  for (size_t i = 0; i< vertices.size(); i++)
    edges.push_back( vertices[i], vertices[(i+1)%vertices.size()] );

  CHECK( (edges.size() == 5) );
  CHECK( (vertices[0].edges().size() == 2) );
  CHECK( (vertices[0].edges(0) == edges[0] ) );
  CHECK( (vertices[1].edges(1) == edges[1] ) );
  CHECK( (vertices[3].edges(0) == edges[2] ) );

  edges.remove( edges[3] );

  CHECK( (edges.size() == 4) );
  CHECK( (vertices[0].edges().size() == 2) );
  CHECK( (vertices[3].edges().size() == 1) );
  CHECK( (vertices[4].edges().size() == 1) );

  CHECK( (vertices[3].edges(0) == edges[2] ) );
  CHECK( (vertices[4].edges(0) == edges[3] ) );

  // Check that pointer equality works
  Edge& e = edges.push_back( vertices[1], vertices[3] );
  CHECK( (e == edges[4]) );
  CHECK( (&e == &edges[4]) );


} // vertices_and_edges()

/*********************************************************************
* Test Container::sort()
*********************************************************************/
void sort()
{
  Container<Vertex> vertices {};

  vertices.push_back( 0.0, 0.0 );
  vertices.push_back( 1.0, 1.0 );
  vertices.push_back( 2.0, 1.0 );
  vertices.push_back( 3.0, 1.0 );
  vertices.push_back( 2.0, 0.0 );

  // Sort vertices by ascending distance to point "x"
  const Vec2d x {1.0,1.0};

  vertices.sort(
  [x] ( Container<Vertex>::value_type& a, 
        Container<Vertex>::value_type& b )
  {
    const double l1 = ((*a).xy()-x).norm_sqr();
    const double l2 = ((*b).xy()-x).norm_sqr();
    return ( l1 < l2 );
  });

  double l_old = 0.0;
  for ( const auto& v_ptr : vertices )
  {
    const double l = ((*v_ptr).xy()-x).norm();
    CHECK( (l_old <= l) );
    l_old = l;
  }

} // sort()

/*********************************************************************
* Test inheritance of generic container class
*********************************************************************/
// A simple class similar to the container class
template <typename T>
class FastList
{
public:

  using List = std::list<std::unique_ptr<T>>;
  using iterator = typename List::iterator;
  using const_iterator = typename List::const_iterator;

  // constructor()
  FastList() {}

  // insert()
  template <typename... Args>
  T& insert( const_iterator pos, Args&&... args )
  {
    std::unique_ptr<T> u_ptr = std::make_unique<T>(args...);
    T* ptr = u_ptr.get();
    iterator iter = items_.insert( pos, std::move(u_ptr) );
    ptr->pos_          = iter;
    ptr->in_container_ = true;
    
    return *ptr;
  }

  // push_back()
  template <typename... Args>
  T& push_back( Args&&... args )
  { return insert( items_.end(), args... ); }

  // Attributes
  List items_;

}; 

// The base class from which we will derive classes that are stored
// in the FastList
template<typename Derived>
class FastListEntry
{
public:
  using List = typename FastList<Derived>::List;
  using Iterator = typename List::iterator;

  FastListEntry(double x, double y) : xy_ {x, y} {}

  const Vec2d xy() const { return xy_; }
  const Iterator& pos() const { return pos_; }
  bool in_container() const { return in_container_; }

  Vec2d    xy_            {};
  Iterator pos_           {};
  bool     in_container_  {false};
}; 

// A dummy class that is inherited from FastListEntry
class InheritedEntry : public FastListEntry<InheritedEntry>
{
public:

  InheritedEntry(double x, double y) : FastListEntry<InheritedEntry>(x,y) 
  {}

}; 

/*********************************************************************
* The actual test function for the inheritance class
*********************************************************************/
void container_entry_inheritance()
{
  // Generic Test
  FastList<InheritedEntry> test_inheritance {};

  auto& t = test_inheritance.push_back( 1.0, 1.0 );

  CHECK( t.in_container_ );

} // container_entry_inheritance()


} // namespace ContainerTests


/*********************************************************************
* Run tests for: Container.h
*********************************************************************/
void run_tests_Container()
{
  ContainerTests::push_back();
  ContainerTests::insert();
  ContainerTests::remove();
  ContainerTests::vertices_and_edges();
  ContainerTests::sort();
  ContainerTests::container_entry_inheritance();

} // run_tests_Container()
