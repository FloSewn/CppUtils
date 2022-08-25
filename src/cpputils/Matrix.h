/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>

namespace CppUtils {


/*********************************************************************
* Inspired by:
* -----------
*  https://stackoverflow.com/questions/20873768/more-efficient-way-\
*  with-multi-dimensional-arrays-matrices-using-c
*********************************************************************/
template<typename T, typename Allocator = std::allocator<T>>
class Matrix
{
public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Matrix(int w, int h) 
  : width_ { w }
  , height_ { h }
  , data_ (w*h, 0) 
  { }

  Matrix(T* data, int w, int h)
  : width_ { w }
  , height_ { h }
  , data_ (w*h, 0) 
  {
    std::copy(&data[0], &data[0] + w*h, const_cast<T*>(data_.data()));
  }

  Matrix(T** data, int w, int h)
  : width_ { w }
  , height_ { h }
  , data_ (w*h, 0) 
  {
    std::copy(data[0], data[0] + w*h, const_cast<T*>(data_.data()));
  }

  virtual ~Matrix() {}

  /*------------------------------------------------------------------
  | Copy / Move
  ------------------------------------------------------------------*/
  Matrix(const Matrix& m)
  : width_ { m.width_ }
  , height_ { m.height_ }
  , data_ { m.data_ }
  { }

  Matrix(Matrix&& m)
  : width_ { m.width_ }
  , height_ { m.height_ }
  , data_ { std::move(m.data_) }
  {}

  /*------------------------------------------------------------------
  | Operators
  ------------------------------------------------------------------*/
  inline Matrix& operator = (const Matrix& m)
  {
    width_  = m.width_;
    height_ = m.height_;
    data_.assign( m.data_.begin(), m.data_.end() );
    return *this;
  }

  inline T* operator[](const int i)
  { return data_.data() + height_*i; }

  inline const T* operator [](const int i) const
  { return data_.data() + height_*i; }

  /*------------------------------------------------------------------
  | Swap data 
  ------------------------------------------------------------------*/
  inline Matrix& swap(Matrix& m)
  {
    width_  = m.width_;
    height_ = m.height_;
    data_.swap( m.data_ );
    return *this;
  }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  inline int width() { return width_; }
  inline int width() const { return width_; }

  inline int height() { return height_; }
  inline int height() const { return height_; }

  inline std::size_t size() { return data_.size(); }
  inline std::size_t size() const { return data_.size(); }

private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  int width_  { 0 };
  int height_ { 0 };

  std::vector<T, Allocator> data_;

};

} // namespace CppUtils
