/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein, 
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
#include "RTreeND.h"
#include "BBoxND.h"

#include "Timer.h"
#include "Helpers.h"
#include "Log.h"

namespace RTreeNDBenchmarks 
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
void benchmark(std::size_t n_samples, bool bulk_insertion=false)
{
  /*------------------------------------------------------------------
  | Initialize 
  ------------------------------------------------------------------*/
  std::vector<Vertex2d> vertices_sorted;
  std::vector<Vertex2d> vertices;
  Timer timer {};
  SpiralFunction fun( 2.5, -3.4, 6.5);
  RTreeND<Vertex2d,TREE_M,double,2> tree {};

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
  timer.count("RTreeND object insertion");

  if ( bulk_insertion )
  {
    std::vector<BBoxND<double,2>> v_bboxes {};

    for (const auto& v : vertices)
      v_bboxes.push_back( v.bbox() );

    tree.insert( vertices, v_bboxes );
  }
  else
  {
    for (const auto& v : vertices)
      tree.insert( v, v.bbox() );
  }

  /*------------------------------------------------------------------
  | Remove random tree data
  ------------------------------------------------------------------*/
  timer.count("RTreeND object removal");

  std::vector<size_t> removal_ids( n_samples / 2 );
  std::iota(removal_ids.begin(), removal_ids.end(), 0);
  std::shuffle(removal_ids.begin(), removal_ids.end(), rng);

  for (int i = 0; i < removal_ids.size(); ++i)
  {
    const auto& v = vertices[removal_ids[i]];
    tree.remove(v, v.bbox() );
  }

  /*------------------------------------------------------------------
  | Finalize benchmark - output times to user
  ------------------------------------------------------------------*/
  timer.count("Finalize");
  
  LOG(INFO) << "\n----------------------------------";

  LOG(INFO) << "n_samples: " << n_samples;

  LOG(INFO) << timer.message(1)  << ": " 
            << timer.delta(1) << "s";
  
  LOG(INFO) << timer.message(2)  << ": " 
            << timer.delta(2) << "s";

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  LOG(INFO) << "";
  LOG(INFO) << "";
  LOG(INFO) << "";
  LOG(INFO) << "";

  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };

  std::string file_name 
  { source_dir + "/auxiliary/benchmark_data/RTree" };

  //RTreeNDWriter writer { tree };

  //writer.write_to_txt( file_name );

  //writer.print(std::cout);


} // benchmark()


} // namespace RTreeNDBenchmarks

/*********************************************************************
* Run benchmarks for RTreeND.h
*********************************************************************/
void run_benchmarks_RTreeND()
{
  std::size_t n_samples = 10000;

  RTreeNDBenchmarks::benchmark<64>(n_samples, false);
  
} // run_benchmarks_RTreeND()
