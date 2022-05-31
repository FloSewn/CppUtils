/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <cassert>

#include "benchmarks.h"

#include "Vec2.h"
#include "Container.h"
#include "Timer.h"
#include "Helpers.h"

namespace ContainerBenchmarks 
{
using namespace CppUtils;

class Edge;

/*********************************************************************
* A spiral function for the generation of vertices in the QuadTree
* benchmark
*********************************************************************/
class SpiralFunction
{
public:
  SpiralFunction(double a, double b, double c)
  : a_ { a }, b_ { b }, c_ { c } {}

  inline Vec2d eval( const double t ) const
  {
    double x = (a_ + b_ * t) * cos(t) + c_ * sin(40.*t);
    double y = (a_ + b_ * t) * sin(t) + c_ * cos(40.*t);
    return { x,y };
  }

  inline void eval( const double t, double& x, double& y)
  {
    x = (a_ + b_ * t) * cos(t) + c_ * sin(40.*t);
    y = (a_ + b_ * t) * sin(t) + c_ * cos(40.*t);
  }

private:
  double a_  { 0.0 };
  double b_  { 0.0 };
  double c_  { 0.0 };

}; 

/*********************************************************************
* A simple vertex class that is stored in the Container
*********************************************************************/
class Vertex
{
public:

  friend Container<Vertex>;
  using List = Container<Vertex>::List; 
  using Iterator = List::iterator;
  using EdgeList = std::list<Edge*>;

  /*------------------------------------------------------------------
  | Constructor 
  ------------------------------------------------------------------*/
  Vertex(double x, double y) : xy_ {x, y} {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const Vec2d xy() const { return xy_; }
  const Iterator& pos() const { return pos_; }
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

  /*------------------------------------------------------------------
  | Container functions 
  ------------------------------------------------------------------*/
  bool in_container() const { return in_container_; }

private:
  /*------------------------------------------------------------------
  | Container functions 
  ------------------------------------------------------------------*/
  void container_destructor() {}

  /*------------------------------------------------------------------
  | Vertex attributes 
  ------------------------------------------------------------------*/
  Vec2d         xy_;
  EdgeList      edges_ {};

  Iterator      pos_ {nullptr};
  bool          in_container_;

}; 


/*********************************************************************
* Benchmark for Container
*********************************************************************/
void benchmark(size_t n) 
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  Container<Vertex> vertices { 1000.0 };
  Timer timer {};
  SpiralFunction fun( 2.5, -3.4, 6.5);

  constexpr double t0 =  0.0;
  constexpr double t1 =  5.0 * M_PI;
  const double dt = (t1-t0) / static_cast<double>(n);

  double x, y;

  /*------------------------------------------------------------------
  | Push back new vertices to the container
  ------------------------------------------------------------------*/
  timer.count("Container - push_back()");

  for (size_t i = 0; i < n; ++i)
  {
    const double t = t0 + i*dt;
    fun.eval(t, x, y);
    vertices.push_back( x, y );
  }

  /*------------------------------------------------------------------
  | Remove vertices from the container
  ------------------------------------------------------------------*/
  timer.count("Container - remove()");

  for (size_t i = 0; i < n; ++i)
  {
    Vertex& v_last = vertices.back();
    vertices.remove( v_last );
  }

  /*------------------------------------------------------------------
  | Finalize benchmark - output times to user
  ------------------------------------------------------------------*/
  timer.count();

  // User output
  CppUtils::SimpleLogger MSG(std::cout, "  ");

  MSG << std::endl;
  MSG << "----------------------------------" << std::endl;
  MSG << "n=" << n << std::endl;
  MSG << "Container vertex insertion          : " 
      << timer.delta(0) << "s" << std::endl;
  MSG << "Container vertex removal            : " 
      << timer.delta(1) << "s" << std::endl;

} // benchmark()

} // namespace ContainerBenchmarks 

/*********************************************************************
* Run benchmarks for Container.h
*********************************************************************/
void run_benchmarks_Container()
{
  ContainerBenchmarks::benchmark( 1000 );
  ContainerBenchmarks::benchmark( 10000 );
  ContainerBenchmarks::benchmark( 100000 );
  ContainerBenchmarks::benchmark( 1000000 );
  
} // run_benchmarks_Container()

