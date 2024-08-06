#pragma once

#include <iostream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include <vector>

namespace pargeo {
  using namespace std;

  struct _empty3 {
    int arr[0]; // todo this produces a struct of size 0 but seems dangerous, need to check
  };

  template <class _dataType, class _floatType, class attributeType> class _vectorPoint {

    static constexpr _dataType empty = numeric_limits<_dataType>::max();

  public:

    int dim;
    typedef _floatType floatType;

    std::vector<_dataType> x;

    _vectorPoint() {}

    _vectorPoint(int dim) : dim(dim), x(dim) {
        for (int i=0; i< dim; ++i) x[i]=empty;
    }

    _vectorPoint(std::vector<_dataType>* p) : dim(p->size()), x(p->size()) { for (int i=0; i<dim; ++i) x[i]=p[i]; }

    _vectorPoint(_vectorPoint* p): dim(p->dim), x(p->dim) { for (int i=0; i<dim; ++i) x[i]=p->x[i]; }

    template<class _tIn>
    _vectorPoint(parlay::slice<_tIn*,_tIn*> p) : dim(p.size()), x(p.size()) {
      for(int i=0; i<dim; ++i) x[i] = (_dataType)p[i];}

    void setEmpty() {x[0]=empty;}

    bool isEmpty() {return x[0]==empty;}

    _vectorPoint operator+(_vectorPoint op2) {
      std::vector<_dataType> xx(dim);
      for (int i=0; i<dim; ++i) xx[i] = x[i]+op2.x[i];
      return _vectorPoint(xx);}

    _vectorPoint operator-(_vectorPoint op2) {
      std::vector<_dataType> xx(dim);
      for (int i=0; i<dim; ++i) xx[i] = x[i]-op2.x[i];
      return _vectorPoint(xx);}

    _vectorPoint operator*(_dataType dv) {
      std::vector<_dataType> xx(dim);
      for (int i=0; i<dim; ++i) xx[i] = x[i]*dv;
      return _vectorPoint(xx);}

    _vectorPoint operator/(_dataType dv) {
      std::vector<_dataType> xx(dim);
      for (int i=0; i<dim; ++i) xx[i] = x[i]/dv;
      return _vectorPoint(xx);}

    _dataType& operator[](int i) {return x[i];}

    _dataType& at(int i) {return x[i];}

    friend bool operator==(_vectorPoint a, _vectorPoint b) {
        if (a.dim != b.dim) return false;
        for (int ii = 0; ii < a.dim; ++ii) {
            if (a[ii] != b[ii]) return false;
        }
        return true;
    }

    friend bool operator!=(_vectorPoint a, _vectorPoint b) {return !(a==b);}

    std::vector<_dataType>* coords() {return x;}

    inline _dataType distSqr(_vectorPoint p) {
      _dataType xx=0;
      for (int i=0; i<dim; ++i) xx += (x[i]-p.x[i])*(x[i]-p.x[i]);
      return xx;}

    inline _floatType dist(_vectorPoint p) {
      return sqrt(distSqr(p));
    }

    _dataType dot(_vectorPoint p2) {
      _dataType r = 0;
      for(int i=0; i<dim; ++i) r += x[i]*p2[i];
      return r;}

    _vectorPoint mult(_dataType c) {
      _vectorPoint r;
      for(int i=0; i<dim; ++i) r[i] = x[i]*c;
      return r;}

  };

  using VectorPoint = _vectorPoint<double, double, _empty3>;

  using fVectorPoint = _vectorPoint<float, float, _empty3>;

  using lVectorPoint = _vectorPoint<long, double, _empty3>;

  template<class _A, class _B>
  _B point3Cast(_B p) {
    _B q;
    for (int i=0; i<p.dim; ++i) q[i] = p[i];
    return q;
  }
}

static std::ostream& operator<<(std::ostream& os, const pargeo::VectorPoint v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

static std::ostream& operator<<(std::ostream& os, const pargeo::fVectorPoint v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

static std::ostream& operator<<(std::ostream& os, const pargeo::lVectorPoint v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}
