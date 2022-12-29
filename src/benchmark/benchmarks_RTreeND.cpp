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
void benchmark(std::size_t n_samples, 
               std::size_t n_query, 
               bool object_removal=false,
               bool bulk_insertion=false,
               bool brute_force=false)
{
  /*------------------------------------------------------------------
  | Initialize 
  ------------------------------------------------------------------*/
  std::vector<Vertex2d> vertices_sorted;
  std::vector<Vertex2d> vertices;
  std::vector<Vertex2d> query_points;
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
  | Nearest-neighbor query points initialization
  ------------------------------------------------------------------*/
  timer.count("Query initialization");

  // Obtain extents of all vertices
  auto v_x_min = std::min_element(vertices.begin(), vertices.end(),
                   [](const Vertex2d& v1, const Vertex2d& v2)
                   { return v1.pos().x < v2.pos().x; });

  auto v_x_max = std::max_element(vertices.begin(), vertices.end(),
                   [](const Vertex2d& v1, const Vertex2d& v2)
                   { return v1.pos().x < v2.pos().x; });

  auto v_y_min = std::min_element(vertices.begin(), vertices.end(),
                   [](const Vertex2d& v1, const Vertex2d& v2)
                   { return v1.pos().y < v2.pos().y; });

  auto v_y_max = std::max_element(vertices.begin(), vertices.end(),
                   [](const Vertex2d& v1, const Vertex2d& v2)
                   { return v1.pos().y < v2.pos().y; });

  double x_min   = (*v_x_min).pos().x;
  double x_max   = (*v_x_max).pos().y;
  double y_min   = (*v_y_min).pos().x;
  double y_max   = (*v_y_max).pos().y;
  double dx_glob = x_max - x_min;
  double dy_glob = y_max - y_min;

  std::random_device               rand_dev;
  std::mt19937                     gen( rand_dev() );
  std::uniform_real_distribution<> distrib(0, 1);

  for (int i = 0; i < n_query; ++i)
  {
    double xq = x_min + dx_glob * distrib(gen);
    double yq = y_min + dy_glob * distrib(gen);

    query_points.push_back( {xq, yq} );
  }

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
  | Query data points
  ------------------------------------------------------------------*/
  typename RTreeND<Vertex2d,TREE_M,double,2>::ObjectDistFunction 
  sqr_dist_fun = [](const VecND<double,2> p, const Vertex2d& v)
  { return (v.pos()-p).norm_sqr(); };

  timer.count("RTreeND nearest query");

  std::vector<const Vertex2d*> winners_rtree {};

  for ( auto& q : query_points )
    winners_rtree.push_back( tree.nearest( q.pos(), sqr_dist_fun ) );


  if ( brute_force )
  {
    timer.count("Brute force nearest query");

    std::vector<const Vertex2d*> winners_bf {};

    for ( auto& q : query_points )
    {
      std::size_t i_winner = 0;
      double max_dist = std::numeric_limits<double>::max();

      for ( std::size_t i = 0; i < vertices.size(); ++i )
      {
        double dist = (vertices[i].pos() - q.pos()).norm_sqr();

        if ( dist < max_dist )
        {
          max_dist = dist;
          i_winner = i;
        }
      }

      winners_bf.push_back( &vertices[i_winner] );

      /*
      std::vector<size_t> ids( vertices.size() );
      std::iota(ids.begin(), ids.end(), 0);

      std::stable_sort(ids.begin(), ids.end(),
        [&vertices, &q](size_t i1, size_t i2)
      {
        double d1 = (vertices[i1].pos() - q.pos()).norm_sqr();
        double d2 = (vertices[i2].pos() - q.pos()).norm_sqr();
        return d1 < d2;
      });

      winners_bf.push_back( &vertices[ids[0]] );
      */
    }


    timer.count("Check for nearest query");
    for (std::size_t i = 0; i < query_points.size(); ++i)
    {
      ASSERT( winners_rtree[i] != nullptr, "ERROR_1");
      ASSERT( winners_bf[i] != nullptr, "ERROR_2");
      ASSERT( winners_rtree[i]->pos() == winners_bf[i]->pos(), "ERROR_3");
    }
  }

  /*------------------------------------------------------------------
  | Remove random tree data
  ------------------------------------------------------------------*/
  if ( object_removal )
  {
    timer.count("RTreeND object removal");

    std::vector<size_t> removal_ids( n_samples / 2 );
    std::iota(removal_ids.begin(), removal_ids.end(), 0);
    std::shuffle(removal_ids.begin(), removal_ids.end(), rng);

    for (int i = 0; i < removal_ids.size(); ++i)
    {
      const auto& v = vertices[removal_ids[i]];
      tree.remove(v, v.bbox() );
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
  std::size_t n_query   = 1000;
  bool object_removal   = false;
  bool bulk_insertion   = true;
  bool brute_force      = true;

  RTreeNDBenchmarks::benchmark<8>(n_samples, n_query, 
                                   object_removal,
                                   bulk_insertion,
                                   brute_force);
  
} // run_benchmarks_RTreeND()
