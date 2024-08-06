// This code is part of the project "Fast Parallel Algorithms for Euclidean
// Minimum Spanning Tree and Hierarchical Spatial Clustering"
// Copyright (c) 2021 Yiqiu Wang, Shangdi Yu, Yan Gu, Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#pragma once

#include <limits> // numeric_limits
#include <algorithm> // nth_element
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "kdTreeVector.h"
#include "../include/hdbscan/vectorpoint.h"

namespace pargeo {

  namespace knnBuf {

    typedef int intT;
    typedef double floatT;

    template <typename T>
    struct elemVector {
      floatT cost;// Non-negative
      T entry;
      elemVector(floatT t_cost, T t_entry) : cost(t_cost), entry(t_entry) {}
      elemVector() : cost(std::numeric_limits<floatT>::max()) {}
      bool operator<(const elemVector& b) const {
        if (cost < b.cost) return true;
        return false;
      }
    };

    template <typename T>
    struct bufferVector {
      typedef parlay::slice<elemVector<T>*, elemVector<T>*> sliceT;
      intT k;
      intT ptr;
      sliceT buf;

      bufferVector(intT t_k, sliceT t_buf): k(t_k), ptr(0), buf(t_buf) {}

      inline void reset() {ptr = 0;}

      bool hasK() {return ptr >= k;}

      elemVector<T> keepK() {
        if (ptr < k) throw std::runtime_error("Error, kbufferVector not enough k.");
        ptr = k;
        std::nth_element(buf.begin(), buf.begin()+k-1, buf.end());
        return buf[k-1];
      }

      void sort() { // todo check
        if (ptr < k) throw std::runtime_error("Error, sorting kbufferVector without enough k.");
        parlay::sort(buf.cut(0, k));
      }

      void insert(elemVector<T> t_elem) {
        buf[ptr++] = t_elem;
        if (ptr >= buf.size()) keepK();
      }

      elemVector<T> operator[](intT i) {
        if (i < ptr) return buf[i];
        else return elemVector<T>();
      }
    };
  }

  using namespace knnBuf;
  using namespace std;

  template<typename nodeT, typename objT>
  void knnRangeHelperVector(nodeT* tree, objT& q, VectorPoint qMin, VectorPoint qMax, double radius, bufferVector<objT*>& out) {
    int relation = tree->boxCompare(qMin, qMax, tree->getMin(), tree->getMax());
    if(relation == tree->boxExclude) {
      return;
    } else if (relation == tree->boxInclude) {
      for (size_t i = 0; i < tree->size(); ++i) {
        objT* p = tree->getItem(i);
        out.insert(elemVector(q.dist(*p), p));
      }
    } else { // intersect
      if (tree->isLeaf()) {
	      for (size_t i = 0; i < tree->size(); ++ i) {
          objT* p = tree->getItem(i);
          double dist = q.dist(*p);
	        if (dist <= radius) {
            out.insert(elemVector(dist, p));
          }
	      }
      } else {
        knnRangeHelperVector<nodeT, objT>(tree->L(), q, qMin, qMax, radius, out);
        knnRangeHelperVector<nodeT, objT>(tree->R(), q, qMin, qMax, radius, out);
      }
    }
  }

  template<typename nodeT, typename objT>
  void knnRangeVector(nodeT* tree, objT& q, double radius, bufferVector<objT*>& out) {
    int dim = q.dim;
    VectorPoint qMin(dim), qMax(dim);
    for (size_t i=0; i<dim; i++) {
      auto tmp = q[i] - radius;
      qMin[i] = tmp;
      qMax[i] = tmp + radius * 2;
    }
    knnRangeHelperVector<nodeT, objT>(tree, q, qMin, qMax, radius, out);
  }

  // modify to take in the tree instead
  template<typename nodeT, typename objT>
  void knnHelperVector(nodeT* tree, objT& q, bufferVector<objT*>& out) {
    // find the leaf first
    int relation = tree->boxCompare(tree->getMin(), tree->getMax(), //This calls functions on the first element in the array
				    VectorPoint(q.coords()),
				    VectorPoint(q.coords()));

    if (relation == tree->boxExclude) {
      return;
    } else {
      if (tree->isLeaf()) {
        // basecase
        for (size_t i = 0; i < tree->size(); ++ i) {
          objT* p = tree->getItem(i);
          out.insert(elemVector(q.dist(*p), p));
        }
      } else {
        knnHelperVector<nodeT, objT>(tree->L(), q, out);
        knnHelperVector<nodeT, objT>(tree->R(), q, out);
      }
    }
    if (!out.hasK()) {
      if (tree->siblin() == NULL) {
	      throw std::runtime_error("Error, knnHelper reached root node without enough neighbors.");
      }
      for (size_t i=0; i<tree->siblin()->size(); ++i) {
        objT* p = tree->siblin()->getItem(i);
        out.insert(elemVector(q.dist(*p), p));
      }
    }else { // bufferVector filled to a least k
      if (tree->siblin() != NULL) {
        elemVector tmp = out.keepK();
        knnRangeVector<nodeT, objT>(tree->siblin(), q, tmp.cost, out);
      }
    }
  }

  template<class objT>
  parlay::sequence<size_t> kdTreeKnnVector(parlay::sequence<objT> &queries, size_t k, kdNodeVector<objT>* tree = nullptr, bool sorted = false) {
    using nodeT = kdNodeVector<objT>;
    bool freeTree = false;

    if (!tree) {
      freeTree = true;
      tree = buildKdtVector<objT>(queries, true);
    }
    auto out = parlay::sequence<elemVector<objT*>>(2*k*queries.size());
    auto idx = parlay::sequence<size_t>(k*queries.size());

    parlay::parallel_for(0, queries.size(), [&](size_t i) {
					      bufferVector buf = bufferVector<objT*>(k, out.cut(i*2*k, (i+1)*2*k));
					      knnHelperVector<nodeT, objT>(tree, queries[i], buf);
					      buf.keepK();
					      if (sorted) buf.sort();
					      for(size_t j=0; j<k; ++j) {
						idx[i*k+j] = buf[j].entry - queries.data();
						//cout << buf[j].cost << endl;
					      }
					      //cout << endl;
					    });
    if (freeTree) free(tree);
    return idx;
  }

  template<typename objT>
  parlay::sequence<size_t> bruteforceKnnVector(parlay::sequence<objT> &queries, size_t k) {
    auto out = parlay::sequence<elemVector<objT*>>(2*k*queries.size());
    auto idx = parlay::sequence<size_t>(k*queries.size());
    parlay::parallel_for(0, queries.size(), [&](size_t i) {
					      objT q = queries[i];
					      bufferVector buf = bufferVector<objT*>(k, out.cut(i*2*k, (i+1)*2*k));
					      for(size_t j=0; j<queries.size(); ++j) {
						objT* p = &queries[j];
						buf.insert(elemVector(q.dist(p), p));
					      }
					      buf.keepK();
					      for(size_t j=0; j<k; ++j) {
						idx[i*k+j] = buf[j].entry - queries.data();
						//cout << buf[j].cost << endl;
					      }
					      //cout << endl;
					    });
    return idx;
  }

} // End namespace pargeo
