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


  template <class _tData, class _tFloat> class _point2 {

    static constexpr _tData empty = numeric_limits<_tData>::max();

  public:
    typedef _tFloat floatT;

    arma::Col<_tData> x;
    int dim;

    _point2() : x(20), dim(20){
      x.fill(empty);
    };
    //Constructors
    _point2(int dim) : x(dim) { 

      x.fill(empty); 
      this->dim = dim;}

    _point2(arma::Col<_tData> *p): x(*p) { 
      dim = x.size();
    }

    _point2(_point2 *p)  { 
      x(p->x); 
      dim = x.size();
    }

    template<class _tIn>
    _point2(parlay::slice<_tIn*,_tIn*> p) {
      x.set_size(p.x.size());
      for(int i=0; i<p.x.size(); ++i) x[i] = (_tData)p[i];
      dim = x.size();
      }

    void setEmpty() {x[0]=empty;}

    bool isEmpty() {return x[0]==empty;}

    _point2 operator+(_point2 op2) {
      return _point2(x+op2.x);}

    _point2 operator-(_point2 op2) {
      return _point2(x-op2.x);}

    _point2 operator*(_tData dv) {
      return _point2(x*dv);} // % is the element wise multiplication in armadillo

    _point2 operator/(_tData dv) {
      return _point2(x/dv);}

    _tData& operator[](int i) {return x[i];}

    _tData& at(int i) {return x[i];}

    friend bool operator==(_point2 a, _point2 b) {  
      return arma::all(a==b);
      }

    friend bool operator!=(_point2 a, _point2 b) {return !(a==b);}

    arma::Col<_tData>* coords() {return &x;}


    inline _tFloat dist(_point2 p) {
      return arma::norm(x-p.x, 2);
    }

    _tData dot(_point2 p2) {
      return arma::dot(x, p2.x);}

    _point2 mult(_tData c) {
      return _point2(x*c);}


    _tFloat length() {
      return arma::norm(x);}
  };


  using point2 = _point2<double, double>;

  using fpoint2 = _point2<float, float>;

  using lpoint2 = _point2<long, double>;

  template<class _A, class _B>
  _B point2Cast(_B p) {
    _B q;
    for (int i=0; i<p.dim; ++i) q[i] = p[i];
    return q;
  }
}

static std::ostream& operator<<(std::ostream& os, const pargeo::point2 v) {
  for (int i=0; i<v.x.size(); ++i)
    os << v.x[i] << " ";
  return os;
}

static std::ostream& operator<<(std::ostream& os, const pargeo::fpoint2 v) {
  for (int i=0; i<v.x.size(); ++i)
    os << v.x[i] << " ";
  return os;
}

static std::ostream& operator<<(std::ostream& os, const pargeo::lpoint2 v) {
  for (int i=0; i<v.x.size(); ++i)
    os << v.x[i] << " ";
  return os;
}
