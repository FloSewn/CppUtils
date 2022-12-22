/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <cassert>
#include <fstream>

#include "benchmarks.h"

#include "CppUtilsConfig.h"

#include "VecND.h"
#include "QuadTree.h"
#include "Timer.h"
#include "Helpers.h"
#include "Log.h"

namespace QuadTreeBenchmarks 
{
using namespace CppUtils;

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
* A simple vertex class that is stored in the QuadTree structure
*********************************************************************/
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

/*********************************************************************
* A simple benchmark for the QuadTree structure
*********************************************************************/
void benchmark(int n, int r, size_t imax, size_t dmax, 
               bool brute_force, bool export_qtree=false)
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  list<unique_ptr<Vertex>> vertices;
  Timer timer {};
  SpiralFunction fun( 2.5, -3.4, 6.5);

  const Vec2d center { 0.0, 0.0 };
  double scale       { 200.0 };
  size_t max_item    { imax };
  size_t max_depth   { dmax };
  QuadTree<Vertex,double> 
    qtree { scale, max_item, max_depth, center };

  constexpr double t0 =  0.0;
  constexpr double t1 =  5.0 * M_PI;
  const double dt = (t1-t0) / static_cast<double>(n);

  double x, y;
  int item_count_qt, item_count_bf;

  // Rectangle size for vertex search (scales with number of vertices)
  double dx = scale / static_cast<double>(n) ;
  double dy = scale / static_cast<double>(n) ;

  
  /*------------------------------------------------------------------
  | Create vertices
  ------------------------------------------------------------------*/
  timer.count("Vertex initialization");

  for (int i = 0; i < n; ++i)
  {
    const double t = t0 + i*dt;
    fun.eval(t, x, y);
    vertices.push_back( make_unique<Vertex>( x, y ) );
  }


  /*------------------------------------------------------------------
  | Add vertices to qtree
  ------------------------------------------------------------------*/
  timer.count("QuadTree initialization");

  for ( const auto& v_ptr : vertices )
    qtree.add( v_ptr.get() );

  /*------------------------------------------------------------------
  | Find objects within rectangle using QuadTree
  ------------------------------------------------------------------*/
  timer.count("QuadTree search");

  item_count_qt = 0;

  for (int i = 0; i < n; ++i)
  {
    const double t = t0 + i*dt;
    fun.eval( t, x, y );

    vector<Vertex*> items {};
    qtree.get_items( {x-dx,y-dy}, {x+dx,y+dy}, items );

    item_count_qt += items.size();
  }

  /*------------------------------------------------------------------
  | Find nearest neighbors via QuadTree
  ------------------------------------------------------------------*/
  timer.count("QuadTree nearest neigbhor search");

  for (int i = 0; i < n; ++i)
  {
    const double t = t0 + i*dt;
    fun.eval( t, x, y );

    qtree.get_nearest({x,y});
  }


  /*------------------------------------------------------------------
  | Find objects within rectangle using brute force
  ------------------------------------------------------------------*/
  timer.count("Brute force search");

  if (brute_force)
  {
    item_count_bf = 0;

    for (int i = 0; i < n; ++i)
    {
      const double t = t0 + i*dt;
      fun.eval( t, x, y );

      for ( const auto& v_ptr : vertices )
      {
        bool in_rect = in_on_rect( 
            (*v_ptr).xy(), {x-dx,y-dy}, {x+dx,y+dy} );

        if ( in_rect )
          ++item_count_bf;
      }
    }

    // Assert that qtree found all items
    ASSERT( (item_count_qt == item_count_bf),
            "QuadTree::get_items() failed.");
  }



  /*------------------------------------------------------------------
  | Find nearest neighbors via brute force
  ------------------------------------------------------------------*/
  timer.count("Brute force nearest neigbhor search");

  if (brute_force)
  {
    for (int i = 0; i < n; ++i)
    {
      const double t = t0 + i*dt;
      fun.eval( t, x, y );

      const Vec2d xy {x,y};

      double best_dist = 10.0 * scale;
      Vertex* winner = nullptr;

      for ( const auto& v_ptr : vertices )
      {
        double dist = (xy-v_ptr->xy()).norm_sqr();

        if (dist < best_dist)
        {
          best_dist = dist;
          winner = v_ptr.get();
        }
      }
    }
  }

  /*------------------------------------------------------------------
  | Randomly remove vertices 
  ------------------------------------------------------------------*/
  timer.count("Remove random vertices");

  if ( r > 0)
  {
    std::srand(123);

    for (int i = 0; i < r; ++i)
    {
      auto it = vertices.begin();
      int pick = (std::rand() % (n-i));
      std::advance(it, pick);

      qtree.remove( it->get() );
      vertices.erase( it );
    }

    ASSERT( (qtree.size() == n-r), 
            "QuadTree::remove() failed.");
  }


  /*------------------------------------------------------------------
  | Repeat search benchmark, to verify that removal 
  | of objects succeeded
  | AGAIN: Find objects within rectangle using qtree 
  ------------------------------------------------------------------*/
  timer.count("QuadTree search after removal");
  
  if ( r > 0)
  {
    item_count_qt = 0;

    for (int i = 0; i < n; ++i)
    {
      const double t = t0 + i*dt;
      fun.eval( t, x, y );

      std::vector<Vertex*> items {};
      qtree.get_items( {x-dx,y-dy}, {x+dx,y+dy}, items );

      item_count_qt += items.size();
    }
  }


  /*------------------------------------------------------------------
  | AGAIN: Find objects within rectangle using brute force
  ------------------------------------------------------------------*/
  timer.count("Brute force search after removal");

  if (brute_force && r > 0)
  {
    item_count_bf = 0;

    for (int i = 0; i < n; ++i)
    {
      const double t = t0 + i*dt;
      fun.eval( t, x, y );

      for ( const auto& v_ptr : vertices )
      {
        bool in_rect = in_on_rect( 
            (*v_ptr).xy(), {x-dx,y-dy}, {x+dx,y+dy} );

        if ( in_rect )
          ++item_count_bf;
      }
    }
    // Assert that qtree found all items
    ASSERT( (item_count_qt == item_count_bf),
            "QuadTree::get_items() failed.");
  }


  /*------------------------------------------------------------------
  | Finalize benchmark - output times to user
  ------------------------------------------------------------------*/
  timer.count("Finalize");

  // User output
  // User output
  LOG(INFO) << "\n----------------------------------";
  LOG(INFO) << "n=" << n << ", r=" << r 
      << ", max. items= " << imax << ", max. depth=" << dmax;
  LOG(INFO) << std::setprecision(7) << std::fixed 
      << "dx=" << dx << ", dy=" << dy;
  LOG(INFO) << "Vertex initialization            : " 
            << timer.delta(0) << "s";
  LOG(INFO) << "QuadTree initialization          : " 
            << timer.delta(1) << "s";
  LOG(INFO) << "QuadTree search                  : " 
            << timer.delta(2) << "s";
  LOG(INFO) << "QuadTree nearest search          : " 
            << timer.delta(3) << "s";

  if (brute_force)
  {
    LOG(INFO) << "BruteForce search                : " 
              << timer.delta(4) << "s";
    LOG(INFO) << "BruteForce nearest search        : " 
              << timer.delta(5) << "s";
    LOG(INFO) << "QuadTree speedup                 : "
              << timer.delta(4) / timer.delta(2);
    LOG(INFO) << "QuadTree speedup (nearest)       : "
              << timer.delta(5) / timer.delta(3);
  }
  if ( r > 0 )
  {
    LOG(INFO) << "Data removal                     : "
              << timer.delta(6);
    LOG(INFO) << "QuadTree search after removal    : " 
              << timer.delta(7) << "s";
  }
  if ( brute_force && r > 0 )
  {
    LOG(INFO) << "BruteForce search after removal  : " 
              << timer.delta(8) << "s";
    LOG(INFO) << "QuadTree speedup after removal   : "
              << timer.delta(8) / timer.delta(7);
  }


  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  if (export_qtree)
  {
    LOG(INFO) << "EXPORT QTREE STRUCTURE";

    std::cout << "VERTICES " << vertices.size() << std::endl;
    for ( const auto& v_ptr : vertices )
      std::cout << std::setprecision(5) << std::fixed 
                << (*v_ptr).xy().x << "," 
                << (*v_ptr).xy().y << std::endl;

    std::cout << "QTREE-LEAFS " << qtree.n_leafs() << std::endl;
    std::cout << qtree;
  }

} // benchmark()

} // namespace QuadTreeBenchmarks



/*********************************************************************
* Run benchmarks for QuadTree.h
*********************************************************************/
void run_benchmarks_QuadTree()
{
  int n = 100;
  size_t e = 11;
  for (size_t i = 0; i < e; ++i)
  {
    QuadTreeBenchmarks::benchmark(n, n/2, 100, 25, true, false);
    n *= 2;
  }
  
} // run_benchmarks_Container()
