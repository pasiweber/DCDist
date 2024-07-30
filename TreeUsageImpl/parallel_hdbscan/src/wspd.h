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

#include <iostream>
#include "parlay/sequence.h"
#include "../include/hdbscan/point.h"
#include "kdTree.h"
#include "parBuf.h"




/*
ChatGPT analysis:

The provided code implements fast parallel algorithms for computing the Well-Separated Pair Decomposition (WSPD) of a point set using kd-trees. WSPD is useful in computational geometry for solving various problems, such as closest pair, nearest neighbor search, and hierarchical clustering.
--------------
Key Components
--------------

Namespace and Includes:

The code is part of the pargeo namespace.
It includes necessary headers for parallel processing, point representation, kd-tree structure, and buffer management.
Data Structures:

wsp: Represents a well-separated pair, containing pointers to two kd-tree nodes.
wspdNormalSerial: A serial implementation for computing WSPD.
wspdNormalParallel: A parallel implementation for computing WSPD.
Functions:

wspdSerial: Computes the WSPD serially.
wspdParallel: Computes the WSPD in parallel.
computeWspdSerial: Helper function to recursively compute WSPD serially.
computeWspdParallel: Helper function to recursively compute WSPD in parallel.
findPairSerial: Finds well-separated pairs serially.
findPairParallel: Finds well-separated pairs in parallel.
geomWellSeparated: Checks if two kd-tree nodes are well-separated geometrically.
------------------
Detailed Breakdown
------------------

Well-Separated Pair (WSP) Structure:
  Defines a simple struct wsp holding two pointers to kd-tree nodes (u and v).

Serial WSPD Computation:

wspdSerial function initializes the computation using computeWspdSerial.
wspdNormalSerial struct handles the serial operations for WSPD, storing the results in a sequence.
computeWspdSerial recursively processes the kd-tree nodes to find well-separated pairs.
findPairSerial checks pairs of nodes to determine if they are well-separated and processes them accordingly.
Parallel WSPD Computation:

wspdParallel function initializes the parallel computation using computeWspdParallel.
wspdNormalParallel struct handles parallel operations, utilizing a buffer for thread-safe results collection.
computeWspdParallel recursively processes the kd-tree nodes in parallel.
findPairParallel checks pairs of nodes to determine if they are well-separated and processes them in parallel.
Geometric Separation Check:

geomWellSeparated checks if two kd-tree nodes are well-separated based on their bounding boxes and a separation parameter s.

----------------------
Separation from pargeo
----------------------

Yes, this functionality could be separated from pargeo. The WSPD computation logic relies on kd-tree structures and parallel processing utilities, which can be implemented independently. The code can be modularized as follows:

Kd-tree Implementation: Move the kd-tree related code to a separate module.
Point and Geometry Utilities: Move point representation and geometric utility functions to a separate module.
Parallel Utilities: Use a general parallel computing library (e.g., OpenMP, TBB, C++17 parallel algorithms) for parallel operations.
WSPD Computation: Implement WSPD computation as an independent module that uses the kd-tree and geometry utilities.

*/




namespace pargeo {

  using namespace parlay;

  // ------------ Well-separated pair struct --------------

  template <typename nodeT>
  struct wsp {
    nodeT* u;
    nodeT* v;
    wsp(nodeT* uu, nodeT* vv): u(uu), v(vv) {}
  };

  // ------------ Functions to call --------------

  template<int dim>
  sequence<wsp<kdNode<dim, point<dim>>>> wspdSerial(kdNode<dim, point<dim>>* tree, double s = 2);

  template<int dim>
  sequence<wsp<kdNode<dim, point<dim>>>> wspdParallel(kdNode<dim, point<dim>>* tree, double s = 2);

  // ------------ Serial implementation --------------

  template<typename nodeT, typename floatT>
  bool geomWellSeparated(nodeT *u, nodeT *v, floatT s=2) {
    floatT circleDiam_u = 0;
    floatT circleDiam_v = 0;
    floatT circleDistance = 0;
    for (int d = 0; d < u->dim; ++ d) {
      floatT uTmpDiff = u->getMax(d) - u->getMin(d);
      floatT vTmpDiff = v->getMax(d) - v->getMin(d);
      floatT uTmpAvg = (u->getMax(d) + u->getMin(d))/2;
      floatT vTmpAvg = (v->getMax(d) + v->getMin(d))/2;
      circleDistance += (uTmpAvg - vTmpAvg) * (uTmpAvg - vTmpAvg);
      circleDiam_u += uTmpDiff * uTmpDiff;
      circleDiam_v += vTmpDiff * vTmpDiff;
    }
    circleDiam_u = sqrt(circleDiam_u);
    circleDiam_v = sqrt(circleDiam_v);
    floatT myRadius = max(circleDiam_u, circleDiam_v)/2;
    circleDistance = sqrt(circleDistance) - circleDiam_u/2 - circleDiam_v/2;
    return circleDistance >= (s * myRadius);
  }

  template <typename nodeT>
  struct wspdNormalSerial {
    using floatT = double;
    using pType = wsp<nodeT>;
    sequence<pType> *out;

    wspdNormalSerial(sequence<pType> *outt) : out(outt) {}

    inline void run(nodeT *u, nodeT *v) {out->emplace_back(u, v);}
    inline bool moveon(nodeT *u, nodeT *v) {return true;}
    inline bool start(nodeT *u) { return true; }
    inline bool wellSeparated(nodeT *u, nodeT *v, floatT s) {
      return geomWellSeparated(u, v, s);
    }
  };

  template<typename nodeT, typename opT>
  inline void findPairSerial(nodeT *u, nodeT *v, opT* f, double s) {
    if (!f->moveon(u, v)) return;

    if (f->wellSeparated(u, v, s)) {
      f->run(u, v);
    } else {
      if (u->isLeaf() && v->isLeaf()) {
	throw std::runtime_error("error, leaves not well separated");
      } else if (u->isLeaf()) {
	findPairSerial(v->L(), u, f, s);
	findPairSerial(v->R(), u, f, s);
      } else if (v->isLeaf()) {
	findPairSerial(u->L(), v, f, s);
	findPairSerial(u->R(), v, f, s);
      } else {
	if (u->lMax() > v->lMax()) {
	  findPairSerial(u->L(), v, f, s);
	  findPairSerial(u->R(), v, f, s);
	} else {
	  findPairSerial(v->L(), u, f, s);
	  findPairSerial(v->R(), u, f, s);}
      }}
  }

  template<typename nodeT, typename opT>
  inline void computeWspdSerial(nodeT *nd, opT *f, double s=2) {
    if (!nd->isLeaf() && f->start(nd)) {
      computeWspdSerial(nd->L(), f, s);
      computeWspdSerial(nd->R(), f, s);
      findPairSerial<nodeT, opT>(nd->L(), nd->R(), f, s);}
  }

  template<int dim>
  sequence<wsp<kdNode<dim, point<dim>>>> wspdSerial(kdNode<dim, point<dim>>* tree, double s) {
    using pointT = point<dim>;
    using nodeT = kdNode<dim, pointT>;
    using pairT = wsp<nodeT>;

    sequence<pairT> out;

    auto wg = wspdNormalSerial<nodeT>(&out);

    computeWspdSerial<nodeT, wspdNormalSerial<nodeT>>(tree, &wg, s);
    //cout << "#-wspd = " << out.size() << endl;
    return out;
  }

  // ------------ Parallel implementation --------------

  template <typename nodeT>
  struct wspdNormalParallel {
    using floatT = double;
    using pType = wsp<nodeT>;
    using bufT = parBuf<pType>;

    bufT **out;

    wspdNormalParallel(size_t n) {
      size_t procs = num_workers();
      out = (bufT**) malloc(sizeof(bufT*)*procs);
      parallel_for(0, procs, [&](size_t p) {
			       out[p] = new bufT(n/procs);
			     });
    }

    ~wspdNormalParallel() {
      size_t procs = num_workers();
      parallel_for(0, procs, [&](size_t p) {
			       delete out[p];});
      free(out);
    }

    inline void run(nodeT *u, nodeT *v) {
      auto tmp = out[worker_id()]->increment();
      tmp->u = u;
      tmp->v = v;
    }

    sequence<pType> collectPairs() {
      int procs = num_workers();
      return parBufCollect<pType>(out, procs);
    }

    inline bool moveon(nodeT *u, nodeT *v) {return true;}
    inline bool start(nodeT *u) { return true; }
    inline bool wellSeparated(nodeT *u, nodeT *v, floatT s) {
      return geomWellSeparated(u, v, s);
    }
  };

  template<typename nodeT, typename opT>
  inline void findPairParallel(nodeT *u, nodeT *v, opT* f, double s) {
    if (!f->moveon(u, v)) return;
    if (u->size() + v->size() < 2000) return findPairSerial(u,v,f,s);
    std::cout << "fja1" <<std::endl;
    if (f->wellSeparated(u, v, s)) {
      f->run(u, v);//need to be thread safe
    } else {
        if (u->isLeaf() && v->isLeaf()) {
          throw std::runtime_error("error, leaves not well separated");
        } else if (u->isLeaf()) {
            std::cout << "fja2" <<std::endl;
            par_do([&](){findPairParallel(v->L(), u, f, s);}, [&](){findPairParallel(v->R(), u, f, s);});
        } else if (v->isLeaf()) {
            std::cout << "fja3" <<std::endl;
            par_do([&](){findPairParallel(u->L(), v, f, s);}, [&](){findPairParallel(u->R(), v, f, s);});
        } else {
          std::cout << "fja4" <<std::endl;
            if (u->lMax() > v->lMax()) {
              par_do([&](){findPairParallel(u->L(), v, f, s);}, [&](){findPairParallel(u->R(), v, f, s);});
            } else {
              par_do([&](){findPairParallel(v->L(), u, f, s);}, [&](){findPairParallel(v->R(), u, f, s);});
            }
        }
      }
  }

  template<typename nodeT, typename opT>
  inline void computeWspdParallel(nodeT *nd, opT *f, double s=2) {
    std::cout << "oo1" << std::endl;

    if (nd->size() < 2000) {
      std::cout << "oo01" << std::endl;
      computeWspdSerial(nd, f, s);
    } else if (!nd->isLeaf() && f->start(nd)) {
      std::cout << "oo2" << std::endl;
      par_do([&](){computeWspdParallel(nd->L(), f, s);},
	     [&](){computeWspdParallel(nd->R(), f, s);});
      findPairParallel<nodeT, opT>(nd->L(), nd->R(), f, s);
    }
  }

  template<int dim>
  sequence<wsp<kdNode<dim, point<dim>>>> wspdParallel(kdNode<dim, point<dim>>* tree, double s) {
    using pointT = point<dim>;
    using nodeT = kdNode<dim, pointT>;
    using pairT = wsp<nodeT>;

    auto wg = wspdNormalParallel<nodeT>(tree->size());

    computeWspdParallel<nodeT, wspdNormalParallel<nodeT>>(tree, &wg, s);
    auto out = wg.collectPairs();
    //cout << "#-wspd-collected = " << out.size() << endl;

    return out;
  }

} // End namespace
