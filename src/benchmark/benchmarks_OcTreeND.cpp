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
#include <algorithm>
#include <random>

#include "benchmarks.h"

#include "CppUtilsConfig.h"

#include "VecND.h"
#include "BBoxND.h"
#include "OcTreeND.h"

#include "Timer.h"
#include "Helpers.h"
#include "Log.h"

namespace OcTreeNDBenchmarks 
{
using namespace CppUtils;

/*********************************************************************
* A 2D spiral function for the generation of vertices 
*********************************************************************/
class SpiralFunction
{
public:
  SpiralFunction(double a, double b, double c)
  : a_ { a }, b_ { b }, c_ { c } {}

  inline void eval( const double t, double& x, double& y) const
  {
    x = (a_ + b_ * t) * cos(t) + c_ * sin(40.*t);
    y = (a_ + b_ * t) * sin(t) + c_ * cos(40.*t);
  }

  inline Vec2d eval( const double t ) const
  {
    Vec2d v {};
    this->eval(t, v.x, v.y);
    return v;
  }

private:
  double a_  { 0.0 };
  double b_  { 0.0 };
  double c_  { 0.0 };
}; 

/*********************************************************************
* A simple vertex class 
*********************************************************************/
template<typename CoordType, std::size_t DIM>
class VertexType
{
public:
  VertexType(const VecND<CoordType,DIM>& pos) 
  : pos_ { pos } 
  , bbox_ { pos, pos }
  {}

  template<
    class... TT,
    class E = std::enable_if_t<(std::is_same_v<TT, CoordType> && ...)>
  >
  VertexType(TT... tt) 
  : pos_ { tt... } 
  {
    bbox_  = { pos_, pos_ };
  }

  const VecND<CoordType,DIM>& pos() const { return pos_; }
  const BBoxND<CoordType,DIM>& bbox() const { return bbox_; }

private:
  VecND<CoordType,DIM>  pos_;
  BBoxND<CoordType,DIM> bbox_;

}; // VertexType

/*********************************************************************
* Typedefs 
*********************************************************************/
using Vertex2d = VertexType<double,2>;


/*********************************************************************
* The actual benchmark  
*********************************************************************/
template<std::size_t TREE_M>
void benchmark(std::size_t n_samples, 
               std::size_t n_query, 
               bool object_removal=false,
               bool brute_force=false)
{
  /*------------------------------------------------------------------
  | Initialize 
  ------------------------------------------------------------------*/
  std::vector<Vertex2d> vertices_sorted;
  std::vector<Vertex2d> vertices;

  Timer timer {};
  SpiralFunction fun( 2.5, -3.4, 6.5);

  // Rectangle size for vertex search (scales with number of vertices)
  double scale       { 200.0 };
  double dx = scale / static_cast<double>(n_query) ;
  double dy = scale / static_cast<double>(n_query) ;


  OcTreeND<Vertex2d,TREE_M,double,2> tree 
  { {{-scale,-scale},{scale,scale}} };

  constexpr double t0 =  0.0;
  constexpr double t1 =  5.0 * M_PI;
  const double dt = (t1-t0) / static_cast<double>(n_samples);
  double x, y;

  /*------------------------------------------------------------------
  | Create vertices
  ------------------------------------------------------------------*/
  timer.count("Vertex initialization");

  for (int i = 0; i < n_samples; ++i)
  {
    const double t = t0 + i*dt;
    fun.eval(t, x, y);
    vertices_sorted.push_back( {x, y} );
  }

  // Shuffle all vertices
  std::vector<size_t> perm_ids( n_samples );
  std::iota(perm_ids.begin(), perm_ids.end(), 0);

  auto rng = std::default_random_engine {};
  std::shuffle(perm_ids.begin(), perm_ids.end(), rng);

  for (int i = 0; i < n_samples; ++i)
    vertices.push_back( vertices_sorted[perm_ids[i]] );

  /*------------------------------------------------------------------
  | Insert tree data
  ------------------------------------------------------------------*/
  timer.count("OcTreeND object insertion");

  for (const auto& v : vertices)
    tree.insert( v, v.pos() );

  /*------------------------------------------------------------------
  | Query data points using the tree
  ------------------------------------------------------------------*/
  timer.count("OcTreeND object query");

  for (std::size_t i = 0; i < n_query; ++i)
  {
    const double t = t0 + i*dt;
    fun.eval( t, x, y );

    auto query = tree.query( {{x-dx,y-dy}, {x+dx,y+dy}} );
  }

  /*------------------------------------------------------------------
  | Query data points using brute force
  ------------------------------------------------------------------*/
  if (brute_force)
  {
    timer.count("Brute force search");

    std::size_t count = 0;

    for (std::size_t i = 0; i < n_query; ++i)
    {
      const double t = t0 + i*dt;
      fun.eval( t, x, y );

      for ( auto& v : vertices )
      {
        BBoxND<double,2> bbox = {{x-dx,y-dy}, {x+dx,y+dy}};
        
        if ( bbox.point_inside_touch( v.pos() ) )
          ++count;
      }
    }

  }

  /*------------------------------------------------------------------
  | Remove random tree data
  ------------------------------------------------------------------*/
  if ( object_removal )
  {
    timer.count("OcTreeND object removal");

    std::vector<size_t> removal_ids( n_samples / 2 );
    std::iota(removal_ids.begin(), removal_ids.end(), 0);
    std::shuffle(removal_ids.begin(), removal_ids.end(), rng);

    for (int i = 0; i < removal_ids.size(); ++i)
    {
      const auto& v = vertices[removal_ids[i]];
      tree.remove(v, v.pos() );
    }
  }

  /*------------------------------------------------------------------
  | Finalize benchmark - output times to user
  ------------------------------------------------------------------*/
  timer.count("Finalize");
  
  LOG(INFO) << "\n----------------------------------";

  LOG(INFO) << "n_samples: " << n_samples;

  for ( std::size_t i = 0; i < timer.size()-1; ++i )
    LOG(INFO) << timer.message(i)  << ": " 
              << timer.delta(i) << "s";

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  LOG(INFO) << "";
  LOG(INFO) << "";
  LOG(INFO) << "";
  LOG(INFO) << "";

  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };

  std::string file_name 
  { source_dir + "/auxiliary/benchmark_data/OcTree" };

  OcTreeNDWriter writer { tree };

  writer.write_to_vtu( file_name );


} // benchmark()


} // namespace OcTreeNDBenchmarks

/*********************************************************************
* Run benchmarks for OcTreeND.h
*********************************************************************/
void run_benchmarks_OcTreeND()
{
  std::size_t n_samples = 10000;
  std::size_t n_query   = 1000;
  bool object_removal   = true;
  bool brute_force      = true;

  OcTreeNDBenchmarks::benchmark<8>(n_samples, n_query, 
                                   object_removal,
                                   brute_force);
  
} // run_benchmarks_OcTreeND()
