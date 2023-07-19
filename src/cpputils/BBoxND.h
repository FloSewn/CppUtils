/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <array>
#include <bitset>

#include "VecND.h"
#include "MathUtility.h"
#include "Helpers.h"

namespace CppUtils {




/*********************************************************************
* This class defines axis aligned boxes, which can be used for 
* algorithms requiring minimum bounded rectangles / bounding boxes
*********************************************************************/
template <typename T, std::size_t N>
class BBoxND
{
  using Vec = VecND<T,N>;

public:

  /*------------------------------------------------------------------ 
  | Constructors
  ------------------------------------------------------------------*/
  BBoxND() {}

  BBoxND(const Vec& ll, const Vec& ur)
  : lowleft_ { ll }
  , upright_ { ur }
  {
    scale_ = (ur-ll).product();

    ASSERT( scale_ >= 0, 
      "BBoxND: Invalid bounding box definition. "
      "Box size equals to " + std::to_string(scale_) );
  }

  /*------------------------------------------------------------------
  | Variables for precision 
  ------------------------------------------------------------------*/
  static inline std::size_t cmp_ulp = 2;

  /*------------------------------------------------------------------ 
  | Getters
  ------------------------------------------------------------------*/
  T scale() const { return scale_; }

  Vec& lowleft() { return lowleft_; }
  const Vec& lowleft() const { return lowleft_; }

  Vec& upright() { return upright_; }
  const Vec& upright() const { return upright_; }

  Vec center() const { T two {2}; return (lowleft_+upright_)/two; }

  /*------------------------------------------------------------------ 
  | Return all bounding box vertices in an array
  | (1<<N) returns a power of two, which in this case gives the 
  | number of bounding box vertices
  | 
  | References:
  | -----------
  | http://graphics.stanford.edu/~seander/bithacks.html\
  | #NextBitPermutation
  ------------------------------------------------------------------*/
  std::array<Vec,(1<<N)> vertices() const
  {
    std::array<Vec,(1<<N)> verts { {} };

    const Vec delta = upright_ - lowleft_;

    std::size_t nv = (1<<N);

    for (std::size_t i=0; i < nv; ++i)
    {
      // Create a bitset to access all bbox vertices, 
      // e.g. for N==2: { (0,0), (0,1), (1,1), (1,0) }
      std::bitset<N> bits = std::bitset<N>((i>>1) ^ i);

      // Fill bitset into a Vec
      Vec dir {};
      for (std::size_t j=0; j < N; ++j)
        dir[j] = bits[j];

      // Finally, compute the vertex coordinates
      verts[i] = lowleft_ + dir * delta;
    }

    return std::move(verts);
  }

  /*------------------------------------------------------------------ 
  | Compute the intersection size between two RTreeBBoxes
  ------------------------------------------------------------------*/
  T bbox_intersection(const BBoxND<T,N>& bbox) const
  {
    const Vec& a = lowleft_;
    const Vec& b = upright_;

    const Vec& p = bbox.lowleft();
    const Vec& q = bbox.upright();

    Vec ll {};
    Vec ur {};

    std::transform(a.cbegin(), a.cend(), p.cbegin(), ll.begin(),
                   [](T ai, T pi) { return std::max<T>(ai,pi); } );
     
    std::transform(b.cbegin(), b.cend(), q.cbegin(), ur.begin(),
                   [](T ai, T pi) { return std::min<T>(ai,pi); } );

    const Vec delta = ur - ll;

    bool intersects = std::all_of(delta.cbegin(), delta.cend(), 
                                  [](T i) { return i >= 0; });

    if ( intersects )
      return delta.product();

    return {};

  } // bbox_intersection()

  /*------------------------------------------------------------------ 
  | Compute the union size between two RTreeBBoxes
  ------------------------------------------------------------------*/
  T bbox_union(const BBoxND<T,N>& bbox) const
  {
    return scale_ + bbox.scale() - bbox_intersection(bbox);
  }

  /*------------------------------------------------------------------ 
  | Compose an BBoxND that includes both this and another BBox
  ------------------------------------------------------------------*/
  BBoxND<T,N> bbox_cover(const BBoxND<T,N>& bbox) const
  {
    Vec min {};
    Vec max {};

    std::transform(lowleft_.cbegin(), lowleft_.cend(), 
                   bbox.lowleft_.cbegin(), min.begin(),
                   [](T ai, T pi) { return std::min<T>(ai,pi); } );

    std::transform(upright_.cbegin(), upright_.cend(), 
                   bbox.upright_.cbegin(), max.begin(),
                   [](T ai, T pi) { return std::max<T>(ai,pi); } );

    return { min, max }; 
  }

  /*------------------------------------------------------------------ 
  | Check if two BBoxNDs intersect - 
  | No intersection is detected if both boxes share edges / faces
  ------------------------------------------------------------------*/
  bool bbox_intersect(const BBoxND<T,N>& bbox) const
  {
    const Vec delta_1 = lowleft_ - bbox.upright_;
    const Vec delta_2 = bbox.lowleft_ - upright_;

    bool is_1 = std::all_of(delta_1.cbegin(), delta_1.cend(), 
                            [](T i) { return i < 0; });

    bool is_2 = std::all_of(delta_2.cbegin(), delta_2.cend(), 
                            [](T i) { return i < 0; });

    return (is_1 && is_2);
  }

  /*------------------------------------------------------------------ 
  | Check if two BBoxNDs intersect - 
  | An intersection is detected too, if both boxes share edges / faces
  ------------------------------------------------------------------*/
  bool bbox_intersect_touch(const BBoxND<T,N>& bbox) const
  {
    const Vec delta_1 = lowleft_ - bbox.upright_;
    const Vec delta_2 = bbox.lowleft_ - upright_;

    bool is_1 = std::all_of(delta_1.cbegin(), delta_1.cend(), 
                            [](T i) { return i <= 0; });

    bool is_2 = std::all_of(delta_2.cbegin(), delta_2.cend(), 
                            [](T i) { return i <= 0; });

    return (is_1 && is_2);
  }

  /*------------------------------------------------------------------ 
  | Check if two BBoxNDs touch (share faces)
  ------------------------------------------------------------------*/
  bool bbox_touch(const BBoxND<T,N>& bbox,
                  std::size_t ulp = BBoxND<T,N>::cmp_ulp) const
  {
    const Vec delta_1 = lowleft_ - bbox.upright_;
    const Vec delta_2 = bbox.lowleft_ - upright_;

    bool is_1 = std::any_of(delta_1.cbegin(), delta_1.cend(), 
    [&](T i) { 
      const T v = std::fabs(i);
      return (v <= std::numeric_limits<T>::epsilon() * v * ulp);
    });

    bool is_2 = std::any_of(delta_2.cbegin(), delta_2.cend(), 
    [&](T i) { 
      const T v = std::fabs(i);
      return (v <= std::numeric_limits<T>::epsilon() * v * ulp);
    });

    return (is_1 || is_2);
  }

  /*------------------------------------------------------------------ 
  | Check if a given BBoxNDs is contained inside this BBoxNDs 
  ------------------------------------------------------------------*/
  bool bbox_inside(const BBoxND<T,N>& bbox) const
  {
    return (  point_inside(bbox.lowleft()) 
           && point_inside(bbox.upright()) );
  }

  /*------------------------------------------------------------------ 
  | Check if a given BBoxNDs is contained inside this BBoxNDs 
  ------------------------------------------------------------------*/
  bool bbox_inside_touch(const BBoxND<T,N>& bbox) const
  {
    return (  point_inside_touch(bbox.lowleft()) 
           && point_inside_touch(bbox.upright()) );
  }

  /*------------------------------------------------------------------ 
  | Check if a point coordinate is within the BBoxND
  ------------------------------------------------------------------*/
  bool point_inside(const Vec& p) const
  {
    const Vec delta_1 = p - lowleft_;
    const Vec delta_2 = upright_ - p;

    bool is_1 = std::all_of(delta_1.cbegin(), delta_1.cend(), 
                            [](T i) { return i > 0; });

    bool is_2 = std::all_of(delta_2.cbegin(), delta_2.cend(), 
                            [](T i) { return i > 0; });

    return (is_1 && is_2);
  }

  /*------------------------------------------------------------------ 
  | Check if a point coordinate is within the BBoxND or on its 
  | surface
  ------------------------------------------------------------------*/
  bool point_inside_touch(const Vec& p) const
  {
    const Vec delta_1 = p - lowleft_;
    const Vec delta_2 = upright_ - p;

    bool is_1 = std::all_of(delta_1.cbegin(), delta_1.cend(), 
                            [](T i) { return i >= 0; });

    bool is_2 = std::all_of(delta_2.cbegin(), delta_2.cend(), 
                            [](T i) { return i >= 0; });

    return (is_1 && is_2);
  }

  /*------------------------------------------------------------------ 
  | Check if a point coordinate touches to surface of the BBoxND
  ------------------------------------------------------------------*/
  bool point_touch(const Vec& p,
                  std::size_t ulp = BBoxND<T,N>::cmp_ulp) const
  {
    const Vec delta_1 = p - lowleft_;
    const Vec delta_2 = upright_ - p;

    bool is_1 = std::any_of(delta_1.cbegin(), delta_1.cend(), 
    [&](T i) { 
      const T v = std::fabs(i);
      return (v <= std::numeric_limits<T>::epsilon() * v * ulp);
    });

    bool is_2 = std::any_of(delta_2.cbegin(), delta_2.cend(), 
    [&](T i) { 
      const T v = std::fabs(i);
      return (v <= std::numeric_limits<T>::epsilon() * v * ulp);
    });

    return (is_1 || is_2);
  }

  /*------------------------------------------------------------------ 
  | Check if this BBoxND intersects with a sphere / circle
  | 
  | Reference:
  | ----------
  | https://stackoverflow.com/questions/401847/circle-rectangle-\
  | collision-detection-intersection
  ------------------------------------------------------------------*/
  bool sphere_intersect(Vec center, T radius)
  {
    Vec closest;

    for ( std::size_t i = 0; i < N; ++i )
      closest[i] = std::clamp( center[i], lowleft_[i], upright_[i] );

    Vec delta = center - closest; 

    return delta.norm_sqr() < radius * radius;
  }

  /*------------------------------------------------------------------ 
  | Compute the squared distance from a given point to the BBoxND
  | 
  | Reference:
  | ----------
  | https://stackoverflow.com/questions/5254838/calculating-\
  | distance-between-a-point-and-a-rectangular-box-nearest-point
  ------------------------------------------------------------------*/
  T point_dist_sqr(const Vec& p) const
  {
    Vec delta {};

    const Vec a = lowleft_ - p;
    const Vec b = p - upright_;

    std::transform(a.cbegin(), a.cend(), b.cbegin(), delta.begin(),
     [](T ai, T bi) { return std::max<T>( {ai, {}, bi} ); });

    return delta.norm_sqr();
  }

  /*------------------------------------------------------------------ 
  | Compute the distance from a given point to the BBoxND
  ------------------------------------------------------------------*/
  T point_dist(const Vec& p) const
  { return sqrt( point_dist_sqr(p) ); }


private:
  Vec  lowleft_ {};
  Vec  upright_ {};
  T    scale_   {};

}; // BBoxND

/*--------------------------------------------------------------------
| BBoxND Equality / ineuality
--------------------------------------------------------------------*/
template <typename T, std::size_t N>
inline bool operator==(const BBoxND<T,N>& a, const BBoxND<T,N>& b)
{
  const VecND<T,N> ll = a.lowleft()-b.lowleft();
  const VecND<T,N> ur = a.upright()-b.upright();
  return (ll.is_zero() && ur.is_zero());
}

/*--------------------------------------------------------------------
| ostream
--------------------------------------------------------------------*/
template <typename T, long unsigned int M>
std::ostream& operator<<(std::ostream& os, 
                         const BBoxND<T,M>& bbox)
{
  return os << bbox.lowleft() << " X " << bbox.upright();
}


} // namespace CppUtils
