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
#include "kdTree.h"
#include "kdTreeArma.h"

#include "../include/hdbscan/point.h"
#include "../include/hdbscan/armapoint.h"

namespace pargeo {

  namespace knnBuf {

    typedef int intT;
    typedef double floatT;

    template <typename T>
    struct elem2 {
      floatT cost;// Non-negative
      T entry;
      elem2(floatT t_cost, T t_entry) : cost(t_cost), entry(t_entry) {}
      elem2() : cost(std::numeric_limits<floatT>::max()) {}
      bool operator<(const elem2& b) const {
	if (cost < b.cost) return true;
	return false;}
    };

    template <typename T>
    struct buffer2 {
      typedef parlay::slice<elem2<T>*, elem2<T>*> sliceT;
      intT k;
      intT ptr;
      sliceT buf;

      buffer2(intT t_k, sliceT t_buf): k(t_k), ptr(0), buf(t_buf) {}

      inline void reset() {ptr = 0;}

      bool hasK() {return ptr >= k;}

      elem2<T> keepK() {
	if (ptr < k) throw std::runtime_error("Error, kbuffer2 not enough k.");
	ptr = k;
	std::nth_element(buf.begin(), buf.begin()+k-1, buf.end());
	return buf[k-1];
      }

      void sort() { // todo check
	if (ptr < k) throw std::runtime_error("Error, sorting kbuffer2 without enough k.");
	parlay::sort(buf.cut(0, k));
      }

      void insert(elem2<T> t_elem) {
	buf[ptr++] = t_elem;
	if (ptr >= buf.size()) keepK();
      }

      elem2<T> operator[](intT i) {
	if (i < ptr) return buf[i];
	else return elem2<T>();
      }
    };
  }

  using namespace knnBuf;
  using namespace std;

  template<typename nodeT, typename objT>
  void knnRangeHelper2(nodeT* tree, objT& q, point2 qMin, point2 qMax, double radius, buffer2<objT*>& out) {
    int relation = tree->boxCompare(qMin, qMax, tree->getMin(), tree->getMax());
    if(relation == tree->boxExclude) {
      return;
    } else if (relation == tree->boxInclude) {
      for (size_t i = 0; i < tree->size(); ++i) {
        objT* p = tree->getItem(i);
        out.insert(elem2(q.dist(*p), p));
      }
    } else { // intersect
      if (tree->isLeaf()) {
	      for (size_t i = 0; i < tree->size(); ++ i) {
          objT* p = tree->getItem(i);
          double dist = q.dist(*p);
	        if (dist <= radius) {
            out.insert(elem2(dist, p));
          }
	      }
      } else {
        knnRangeHelper2<nodeT, objT>(tree->L(), q, qMin, qMax, radius, out);
        knnRangeHelper2<nodeT, objT>(tree->R(), q, qMin, qMax, radius, out);
      }
    }
  }

  template<typename nodeT, typename objT>
  void knnRange2(nodeT* tree, objT& q, double radius, buffer2<objT*>& out) {
    int dim = q.dim;
    point2 qMin(dim), qMax(dim);
    for (size_t i=0; i<dim; i++) {
      auto tmp = q[i] - radius;
      qMin[i] = tmp;
      qMax[i] = tmp + radius * 2;
    }
    knnRangeHelper2<nodeT, objT>(tree, q, qMin, qMax, radius, out);
  }

  // modify to take in the tree instead
  template<typename nodeT, typename objT>
  void knnHelper2(nodeT* tree, objT& q, buffer2<objT*>& out) {
    // find the leaf first
    int relation = tree->boxCompare(tree->getMin(), tree->getMax(), //This calls functions on the first element in the array
				    point2(q.coords()),
				    point2(q.coords()));

    if (relation == tree->boxExclude) {
      return;
    } else {
      if (tree->isLeaf()) {
        // basecase
        for (size_t i = 0; i < tree->size(); ++ i) {
          objT* p = tree->getItem(i);
          out.insert(elem2(q.dist(*p), p));
        }
      } else {
        knnHelper2<nodeT, objT>(tree->L(), q, out);
        knnHelper2<nodeT, objT>(tree->R(), q, out);
      }
    }
    if (!out.hasK()) {
      if (tree->siblin() == NULL) {
	      throw std::runtime_error("Error, knnHelper reached root node without enough neighbors.");
      }
      for (size_t i=0; i<tree->siblin()->size(); ++i) {
        objT* p = tree->siblin()->getItem(i);
        out.insert(elem2(q.dist(*p), p));
      }
    }else { // buffer2 filled to a least k
      if (tree->siblin() != NULL) {
        elem2 tmp = out.keepK();
        knnRange2<nodeT, objT>(tree->siblin(), q, tmp.cost, out);
      }
    }
  }

  template<class objT>
  parlay::sequence<size_t> kdTreeKnn2(parlay::sequence<objT> &queries, size_t k, kdNode2<objT>* tree = nullptr, bool sorted = false) {
    using nodeT = kdNode2<objT>;
    bool freeTree = false;

    if (!tree) {
      freeTree = true;
      tree = buildKdt2<objT>(queries, true);
    }
    auto out = parlay::sequence<elem2<objT*>>(2*k*queries.size());
    auto idx = parlay::sequence<size_t>(k*queries.size());

    parlay::parallel_for(0, queries.size(), [&](size_t i) {
					      buffer2 buf = buffer2<objT*>(k, out.cut(i*2*k, (i+1)*2*k));
					      knnHelper2<nodeT, objT>(tree, queries[i], buf);
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
  parlay::sequence<size_t> bruteforceKnn2(parlay::sequence<objT> &queries, size_t k) {
    auto out = parlay::sequence<elem2<objT*>>(2*k*queries.size());
    auto idx = parlay::sequence<size_t>(k*queries.size());
    parlay::parallel_for(0, queries.size(), [&](size_t i) {
					      objT q = queries[i];
					      buffer2 buf = buffer2<objT*>(k, out.cut(i*2*k, (i+1)*2*k));
					      for(size_t j=0; j<queries.size(); ++j) {
						objT* p = &queries[j];
						buf.insert(elem2(q.dist(p), p));
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
