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
#include <limits>
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
constexpr size_t CPPUTILS_ULP    = 2;

/*********************************************************************
* USEFUL FUNCTIONS
*********************************************************************/
template <typename T> 
static inline T ABS(T a)
{ return ( (a > T{}) ? a : -a); }

template <typename T>
static inline bool EQ0(T a, size_t ulp=CPPUTILS_ULP)
{ 
  T a2 = a*a; 
  return (a2 <= std::numeric_limits<T>::epsilon()*a2*ulp);
}

template <typename T>
static inline bool EQ(T a, T b, size_t ulp=CPPUTILS_ULP)
{ return EQ0(a-b, ulp); }


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
