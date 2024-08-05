#pragma once

#include <iostream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include "parlay/parallel.h"
#include "parlay/primitives.h"

namespace pargeo {
  using namespace std;

  struct _empty {
    int arr[0]; // todo this produces a struct of size 0 but seems dangerous, need to check
  };

  template <int _dim, class _dataType, class _floatType, class attributeType> class _point {

    static constexpr _dataType empty = numeric_limits<_dataType>::max();

  public:

    static constexpr int dim = _dim;
    typedef _floatType floatType;

    _dataType x[_dim];
    attributeType attribute;

    _point() { for (int i=0; i<_dim; ++i) x[i]=empty; }

    _point(_dataType* p) { for (int i=0; i<_dim; ++i) x[i]=p[i]; }

    _point(_point* p): attribute(p->attribute) { for (int i=0; i<_dim; ++i) x[i]=p->x[i]; }

    template<class _tIn>
    _point(parlay::slice<_tIn*,_tIn*> p) {
      for(int i=0; i<_dim; ++i) x[i] = (_dataType)p[i];}

    void setEmpty() {x[0]=empty;}

    bool isEmpty() {return x[0]==empty;}

    _point operator+(_point op2) {
      _dataType xx[_dim];
      for (int i=0; i<_dim; ++i) xx[i] = x[i]+op2.x[i];
      return _point(xx);}

    _point operator-(_point op2) {
      _dataType xx[_dim];
      for (int i=0; i<_dim; ++i) xx[i] = x[i]-op2.x[i];
      return _point(xx);}

    _point operator*(_dataType dv) {
      _dataType xx[_dim];
      for (int i=0; i<_dim; ++i) xx[i] = x[i]*dv;
      return _point(xx);}

    _point operator/(_dataType dv) {
      _dataType xx[_dim];
      for (int i=0; i<_dim; ++i) xx[i] = x[i]/dv;
      return _point(xx);}

    _dataType& operator[](int i) {return x[i];}

    _dataType& at(int i) {return x[i];}

    friend bool operator==(_point a, _point b) {
      for (int ii=0; ii<dim; ++ii) {
	      if (a[ii] != b[ii]) return false;}
      return true;}

    friend bool operator!=(_point a, _point b) {return !(a==b);}

    _dataType* coords() {return x;}

    inline _dataType distSqr(_point p) {
      _dataType xx=0;
      for (int i=0; i<_dim; ++i) xx += (x[i]-p.x[i])*(x[i]-p.x[i]);
      return xx;}

    inline _floatType dist(_point p) {
      return sqrt(distSqr(p));
    }

    _dataType dot(_point p2) {
      _dataType r = 0;
      for(int i=0; i<dim; ++i) r += x[i]*p2[i];
      return r;}

    _point mult(_dataType c) {
      _point r;
      for(int i=0; i<dim; ++i) r[i] = x[i]*c;
      return r;}

  };

  template<int dim>
  using point = _point<dim, double, double, _empty>;

  template<int dim>
  using fpoint = _point<dim, float, float, _empty>;

  template<int dim>
  using lpoint = _point<dim, long, double, _empty>;

  template<class _A, class _B>
  _B pointCast(_B p) {
    _B q;
    for (int i=0; i<p.dim; ++i) q[i] = p[i];
    return q;
  }
}

template<int dim>
static std::ostream& operator<<(std::ostream& os, const pargeo::point<dim> v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

template<int dim>
static std::ostream& operator<<(std::ostream& os, const pargeo::fpoint<dim> v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

template<int dim>
static std::ostream& operator<<(std::ostream& os, const pargeo::lpoint<dim> v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}
