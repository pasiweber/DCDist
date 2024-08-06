#pragma once

#include <iostream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include <armadillo>

namespace pargeo {
  using namespace std;


  template <class _dataType, class _floatType> class _armaPoint {

    static constexpr _dataType empty = numeric_limits<_dataType>::max();

  public:
    typedef _floatType floatType;

    arma::Col<_dataType> x;
    int dim;

    _armaPoint() : x(), dim(0) {
    };
    //Constructors
    _armaPoint(int dim) : x(dim) { 

      x.fill(empty); 
      this->dim = dim;}

    _armaPoint(arma::Col<_dataType> *p): x(*p) { 
      dim = x.size();
    }

    _armaPoint(_armaPoint *p)  { 
      x(p->x); 
      dim = x.size();
    }

    template<class _tIn>
    _armaPoint(parlay::slice<_tIn*,_tIn*> p) {
      x.set_size(p.x.size());
      for(int i=0; i<p.x.size(); ++i) x[i] = (_dataType)p[i];
      dim = x.size();
      }

    void setEmpty() {x[0]=empty;}

    bool isEmpty() {return x[0]==empty;}

    _armaPoint operator+(_armaPoint op2) {
      return _armaPoint(x+op2.x);}

    _armaPoint operator-(_armaPoint op2) {
      return _armaPoint(x-op2.x);}

    _armaPoint operator*(_dataType dv) {
      return _armaPoint(x*dv);} // % is the element wise multiplication in armadillo

    _armaPoint operator/(_dataType dv) {
      return _armaPoint(x/dv);}

    _dataType& operator[](int i) {return x[i];}

    _dataType& at(int i) {return x[i];}

    friend bool operator==(_armaPoint a, _armaPoint b) {  
      return arma::all(a==b);
      }

    friend bool operator!=(_armaPoint a, _armaPoint b) {return !(a==b);}

    arma::Col<_dataType>* coords() {return &x;}


    inline _floatType dist(_armaPoint p) {
      return arma::norm(x-p.x, 2);
    }

    _dataType dot(_armaPoint p2) {
      return arma::dot(x, p2.x);}

    _armaPoint mult(_dataType c) {
      return _armaPoint(x*c);}

  };


  using ArmaPoint = _armaPoint<double, double>;

  using fArmaPoint = _armaPoint<float, float>;

  using lArmaPoint = _armaPoint<long, double>;

  template<class _A, class _B>
  _B point2Cast(_B p) {
    _B q;
    for (int i=0; i<p.dim; ++i) q[i] = p[i];
    return q;
  }
}

static std::ostream& operator<<(std::ostream& os, const pargeo::ArmaPoint v) {
  for (int i=0; i<v.x.size(); ++i)
    os << v.x[i] << " ";
  return os;
}

static std::ostream& operator<<(std::ostream& os, const pargeo::fArmaPoint v) {
  for (int i=0; i<v.x.size(); ++i)
    os << v.x[i] << " ";
  return os;
}

static std::ostream& operator<<(std::ostream& os, const pargeo::lArmaPoint v) {
  for (int i=0; i<v.x.size(); ++i)
    os << v.x[i] << " ";
  return os;
}
