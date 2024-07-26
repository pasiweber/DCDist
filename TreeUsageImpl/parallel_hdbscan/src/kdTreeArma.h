#pragma once

#include "../include/hdbscan/point.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include <mlpack/core.hpp>
namespace pargeo {

template <int _dim, class _objT> class kdNode2 {

  typedef int intT;
  typedef double floatT;
  typedef pargeo::point<_dim> pointT;
  typedef kdNode<_dim, _objT> nodeT;

  // Data fields
  intT id;

  int k;

  pointT pMin, pMax;
  arma::vec pMin2, pMax2;

  nodeT *left;

  nodeT *right;

  nodeT *sib;

  parlay::slice<_objT **, _objT **> items;
  parlay::slice<arma::vec **, arma::vec **> items2; //Currently only using ONE pointer here instead of the pointer pointer.

  // Methods

  inline void minCoords(pointT &_pMin, pointT &p) {
    for (int i = 0; i < _pMin.dim; ++i)
      _pMin[i] = min(_pMin[i], p[i]);
  }

  inline void maxCoords(pointT &_pMax, pointT &p) {
    for (int i = 0; i < _pMax.dim; ++i)
      _pMax[i] = max(_pMax[i], p[i]);
  }

    //New
    inline void minCoords2(arma::vec &_pMin, arma::vec &p) {
        for (int i = 0; i < _pMin.n_elem; ++i){
            _pMin[i] = min(_pMin[i], p[i]);
        }
    }
    //New
    inline void maxCoords2(arma::vec &_pMax, arma::vec &p) {
        for (int i = 0; i < _pMax.n_elem; ++i){
            _pMax[i] = max(_pMax[i], p[i]);
        }
    }




  inline void boundingBoxSerial() {
    pMin = pointT(items[0]->coords());
    pMax = pointT(items[0]->coords());
    for (intT i = 0; i < size(); ++i) {
      minCoords(pMin, items[i][0]);
      maxCoords(pMax, items[i][0]);
    }
  }

    //New
    inline void boundingBoxSerial2() {
        pMin2 = arma::vec(*(items2[0]));
        pMax2 = arma::vec(*(items2[0]));

        for (intT i = 0; i < size(); ++i) {
            minCoords2(pMin2, items2[i][0]);
            maxCoords2(pMax2, items2[i][0]);
        }
    }


  inline void boundingBoxParallel() {
    intT P = parlay::num_workers() * 8;
    intT blockSize = (size() + P - 1) / P;
    pointT localMin[P];
    pointT localMax[P];
    for (intT i = 0; i < P; ++i) {
      localMin[i] = pointT(items[0]->coords());
      localMax[i] = pointT(items[0]->coords());
    }
    parlay::parallel_for(0, P, [&](intT p) {
      intT s = p * blockSize;
      intT e = min((intT)(p + 1) * blockSize, size());
      for (intT j = s; j < e; ++j) {
        minCoords(localMin[p], items[j][0]);
        maxCoords(localMax[p], items[j][0]);
      }
    });
    pMin = pointT(items[0]->coords());
    pMax = pointT(items[0]->coords());
    for (intT p = 0; p < P; ++p) {
      minCoords(pMin, localMin[p]);
      maxCoords(pMax, localMax[p]);
    }
  }



    inline void boundingBoxParallel2() {
        intT P = parlay::num_workers() * 8;
        intT blockSize = (size() + P - 1) / P;
        pointT localMin[P];
        pointT localMax[P];
        for (intT i = 0; i < P; ++i) {
            localMin[i] = pointT(items[0]->coords());
            localMax[i] = pointT(items[0]->coords());
        }
        parlay::parallel_for(0, P, [&](intT p) {
            intT s = p * blockSize;
            intT e = min((intT)(p + 1) * blockSize, size());
            for (intT j = s; j < e; ++j) {
            minCoords(localMin[p], items[j][0]);
            maxCoords(localMax[p], items[j][0]);
            }
        });
        pMin = pointT(items[0]->coords());
        pMax = pointT(items[0]->coords());
        for (intT p = 0; p < P; ++p) {
            minCoords(pMin, localMin[p]);
            maxCoords(pMax, localMax[p]);
        }
    }






  inline intT splitItemSerial(floatT xM) {
    if (size() < 2) {
      throw std::runtime_error("Error, kdTree splitting singleton.");
    }
    intT lPt = 0;
    intT rPt = size() - 1;
    while (lPt < rPt) {
      if (items[lPt]->at(k) >= xM) {
        while (items[rPt]->at(k) >= xM && lPt < rPt) {
          rPt--;
        }
        if (lPt < rPt) {
          swap(items[lPt], items[rPt]);
          rPt--;
        } else {
          break;
        }
      }
      lPt++;
    }
    if (items[lPt]->at(k) < xM)
      lPt++;
    return lPt;
  }

  inline bool itemInBox(pointT pMin1, pointT pMax1, _objT *item) {
    for (int i = 0; i < _dim; ++i) {
      if (pMax1[i] < item->at(i) || pMin1[i] > item->at(i))
        return false;
    }
    return true;
  }

  intT findWidest() {
    floatT xM = -1;
    for (int kk = 0; kk < _dim; ++kk) {
      if (pMax[kk] - pMin[kk] > xM) {
        xM = pMax[kk] - pMin[kk];
        k = kk;
      }
    }
    return k;
  }

  void constructSerial(nodeT *space, intT leafSize) {
    boundingBoxSerial();
    sib = NULL;
    if (size() <= leafSize) {
      left = NULL;
      right = NULL;
    } else {
      intT k = findWidest();
      floatT xM = (pMax[k] + pMin[k]) / 2;

      // Split items by xM (serial)
      intT median = splitItemSerial(xM);

      if (median == 0 || median == size()) {
        median = ceil(size() / 2.0);
      }

      /* if (!space[0].isEmpty() || !space[2*median-1].isEmpty()) { */
      /*   throw std::runtime_error("Error, kdNode overwrite."); */
      /* } */

      // Recursive construction
      space[0] = nodeT(items.cut(0, median), median, space + 1, leafSize);
      space[2 * median - 1] = nodeT(items.cut(median, size()), size() - median,
                                    space + 2 * median, leafSize);
      left = space;
      right = space + 2 * median - 1;
      left->sib = right;
      right->sib = left;
    }
  }

  void constructParallel(nodeT *space, parlay::slice<bool *, bool *> flags,
                         intT leafSize) {
    boundingBoxParallel();

    sib = NULL;
    if (size() <= leafSize) {
      left = NULL;
      right = NULL;
    } else {
      intT k = findWidest();
      floatT xM = (pMax[k] + pMin[k]) / 2;

      // Split items by xM in dim k (parallel)
      parlay::parallel_for(0, size(), [&](intT i) {
        if (items[i]->at(k) < xM)
          flags[i] = 1;
        else
          flags[i] = 0;
      });
      auto mySplit = parlay::internal::split_two(items, flags);
      auto splited = mySplit.first;
      intT median = mySplit.second;
      parlay::parallel_for(0, size(),
                           [&](intT i) { items[i] = splited[i]; }); // Copy back

      if (median == 0 || median == size()) {
        median = (size() / 2.0);
      }

      /* if (!space[0].isEmpty() || !space[2*median-1].isEmpty()) { */
      /*   throw std::runtime_error("Error, kdNode overwrite."); */
      /* } */

      // Recursive construction
      parlay::par_do(
          [&]() {
            space[0] = nodeT(items.cut(0, median), median, space + 1,
                             flags.cut(0, median), leafSize);
          },
          [&]() {
            space[2 * median - 1] =
                nodeT(items.cut(median, size()), size() - median,
                      space + 2 * median, flags.cut(median, size()), leafSize);
          });
      left = space;
      right = space + 2 * median - 1;
      left->sib = right;
      right->sib = left;
    }
  }

public:
  using objT = _objT;

  static constexpr int dim = _dim;

  inline nodeT *L() { return left; }

  inline nodeT *R() { return right; }

  inline nodeT *siblin() { return sib; }

  inline intT size() { return items.size(); }

  inline _objT *operator[](intT i) { return items[i]; }

  inline _objT *at(intT i) { return items[i]; }

  inline bool isLeaf() { return !left; }

  inline _objT *getItem(intT i) { return items[i]; }

  inline pointT getMax() { return pMax; }

  inline pointT getMin() { return pMin; }

  inline floatT getMax(int i) { return pMax[i]; }

  inline floatT getMin(int i) { return pMin[i]; }

  bool hasId() { return id != -1; }

  void setId(intT _id) { id = _id; }

  void resetId() { id = -1; }

  intT getId() { return id; }

  static const int boxInclude = 0;

  static const int boxOverlap = 1;

  static const int boxExclude = 2;

  floatT diag() {
    floatT result = 0;
    for (int d = 0; d < _dim; ++d) {
      floatT tmp = pMax[d] - pMin[d];
      result += tmp * tmp;
    }
    return sqrt(result);
  }

  inline floatT lMax() {
    floatT myMax = 0;
    for (int d = 0; d < _dim; ++d) {
      floatT thisMax = pMax[d] - pMin[d];
      if (thisMax > myMax) {
        myMax = thisMax;
      }
    }
    return myMax;
  }

  inline int boxCompare(pointT pMin1, pointT pMax1, pointT pMin2,
                        pointT pMax2) {
    bool exclude = false;
    bool include = true; // 1 include 2
    for (int i = 0; i < _dim; ++i) {
      if (pMax1[i] < pMin2[i] || pMin1[i] > pMax2[i])
        exclude = true;
      if (pMax1[i] < pMax2[i] || pMin1[i] > pMin2[i])
        include = false;
    }
    if (exclude)
      return boxExclude;
    else if (include)
      return boxInclude;
    else
      return boxOverlap;
  }

  kdNode2(parlay::slice<_objT **, _objT **> itemss, intT nn, nodeT *space,
         parlay::slice<bool *, bool *> flags, intT leafSize = 16)
      : items(itemss) {
    resetId();
    if (size() > 2000)
      constructParallel(space, flags, leafSize);
    else
      constructSerial(space, leafSize);
  }

  kdNode2(parlay::slice<_objT **, _objT **> itemss, intT nn, nodeT *space,
         intT leafSize = 16)
      : items(itemss) {
    resetId();
    constructSerial(space, leafSize);
  }
};

template <typename nodeT> inline double nodeDistance2(nodeT *n1, nodeT *n2) {
  using floatT = typename nodeT::objT::floatT;

  for (int d = 0; d < n1->dim; ++d) {
    if (n1->getMin(d) > n2->getMax(d) || n2->getMin(d) > n1->getMax(d)) {
      // disjoint at dim d, and intersect on dim < d
      floatT rsqr = 0;
      for (int dd = d; dd < n1->dim; ++dd) {
        floatT tmp = max(n1->getMin(dd) - n2->getMax(dd),
                         n2->getMin(dd) - n1->getMax(dd));
        tmp = max(tmp, (floatT)0);
        rsqr += tmp * tmp;
      }
      return sqrt(rsqr);
    }
  }
  return 0; // could be intersecting
}

template <typename nodeT> inline double nodeFarDistance2(nodeT *n1, nodeT *n2) {
  using floatT = typename nodeT::objT::floatT;
  floatT result = 0;
  for (int d = 0; d < n1->dim; ++d) {
    floatT tmp =
        max(n1->getMax(d), n2->getMax(d)) - min(n1->getMin(d), n2->getMin(d));
    result += tmp * tmp;
  }
  return sqrt(result);
}

template <int dim, class objT>
kdNode<dim, objT> *buildKdt2(parlay::sequence<objT> &P, bool parallel = true,
                            bool noCoarsen = false,
                            parlay::sequence<objT *> *items = nullptr) {
  typedef kdNode<dim, objT> nodeT;

  size_t n = P.size();

  if (!items) {
    items = new parlay::sequence<objT *>(n);
  }

  parlay::parallel_for(0, n, [&](size_t i) { items->at(i) = &P[i]; });

  parlay::slice<objT **, objT **> itemSlice = items->cut(0, items->size());

  auto root = (nodeT *)malloc(sizeof(nodeT) * (2 * n - 1));

  /* parlay::parallel_for(0, 2*n-1, [&](size_t i) { */
  /* 	root[i].setEmpty(); */
  /*   }); */

  if (parallel) {
    auto flags = parlay::sequence<bool>(n);
    auto flagSlice = parlay::slice(flags.begin(), flags.end());
    root[0] = nodeT(itemSlice, n, root + 1, flagSlice, noCoarsen ? 1 : 16);
  } else {
    root[0] = nodeT(itemSlice, n, root + 1, noCoarsen ? 1 : 16);
  }
  return root;
}

} // namespace pargeo
