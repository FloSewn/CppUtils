/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <functional>
#include <float.h>
#include <limits.h>
#include <vector>
#include <numeric>
#include <algorithm>

namespace CppUtils {

/*********************************************************************
* GLOBAL CONSTANTS
*********************************************************************/
constexpr double CPPUTILS_SMALL  = 1.0E-13; //DBL_EPSILON;
constexpr double CPPUTILS_MAX    = DBL_MAX;
constexpr double CPPUTILS_MIN    = DBL_MIN;

/*********************************************************************
* USEFUL FUNCTIONS
*********************************************************************/
template <typename T> 
static inline T ABS(T a)
{ return ( (a > 0) ? a : -a); }

template <typename T>
static inline bool EQ(T a, T b)
{ return ( ABS(a-b) < CPPUTILS_SMALL ); }

template <typename T>
static inline bool EQ0(T a)
{ return ( ABS(a) < CPPUTILS_SMALL ); }

template <typename T>
static inline T MIN(T a, T b) 
{ return a < b ? a : b; }

template <typename T>
static inline T MAX(T a, T b) 
{ return a > b ? a : b; }

template <typename T>
static inline T MOD(T n, T M)
{ return ((n % M) + M) % M;  }

template <typename T>
static inline T CLAMP(T x, T lower, T upper)
{ return MAX(lower, MIN(upper, x)); }

/*********************************************************************
* Argsort
*
* Reference: 
* ---------
*   https://stackoverflow.com/questions/1577475/c-sorting-\
*   and-keeping-track-of-indexes
*********************************************************************/
template <typename T>
std::vector<size_t> argsort(const std::vector<T> &v)
{
  // Initial index locations
  std::vector<size_t> index( v.size() );
  std::iota(index.begin(), index.end(), 0);

  // Use std::stable_sort() instead of sort to avoid unnecessary
  // index re-orderings when v contains elements of equal values
  std::stable_sort(index.begin(), index.end(), 
    [&v](size_t i1, size_t i2) 
  { 
    return v[i1] < v[i2]; 
  });

  return index;
}


} // namespace CppUtils
