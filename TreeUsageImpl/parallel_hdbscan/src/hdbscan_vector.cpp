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

#include <tuple>
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "getTime.h"
#include "wspd.h"
#include "kdTree.h"
#include "kdTreeKnn.h"
#include "bccp.h"
#include "kruskal.h"
#include "wspdFilter.h"
#include "mark.h"

#include "../include/hdbscan/hdbscan.h"
#include "../include/hdbscan/vectorpoint.h"
#include "kdTreeVector.h"
#include "kdTreeKnnVector.h"


using namespace std;
using namespace parlay;
using namespace pargeo;
using namespace pargeo::hdbscanInternal;

parlay::sequence<pargeo::wghEdge> pargeo::hdbscan_vector(parlay::sequence<pargeo::VectorPoint> &S, size_t minPts) {
  using pointType = VectorPoint; //NEW
  using treeNode = kdNodeVector<pointType>; //NEW
  using floatType = typename pointType::floatType;

  if (S.size() < 2) {
    throw std::runtime_error("need more than 2 points");
  }

  timer t0;
  t0.start();

  //Build the kd-tree for k nearest neighbors
    treeNode* tree = buildKdtVector<pointType>(S, true, true); //This returns a c array of nodes //NEW

    cout << "build-tree-time = " << t0.get_next() << endl;

    //Get the k nearest neighbors over the tree using the kd-tree from above
    sequence<size_t> nns = kdTreeKnnVector<pointType>(S, minPts, tree, true); //NEW

  //It looks like they do it in the same way that we do - first get all k, then take the k'th as the core distance. 
  //They use some nested indexing - i.e. they get the k'th nearest neighbor and compute the distance to that point from the current point and insert it into their list of core distances.
  sequence<floatType> coreDist = sequence<floatType>(S.size());
  parallel_for (0, S.size(), [&](intT i) {
			       coreDist[i] = S[nns[i*minPts + minPts-1]].dist(S[i]);
            //  if(i % 20 == 0){
            //   std::cout << "i:" <<  i << ", " << coreDist[i] << std::endl;
            //  }
			     });

  cout << "core-dist-time = " << t0.get_next() << endl;

  sequence<floatType> cdMin = sequence<floatType>(tree->size() * 2);
  sequence<floatType> cdMax = sequence<floatType>(tree->size() * 2);
  parallel_for(0, tree->size()*2, [&](intT i) {
      cdMin[i] = std::numeric_limits<floatType>::max();
      cdMax[i] = std::numeric_limits<floatType>::lowest();
    });

  //WSPD which produces the edges over which Kruskal's algorithm will be run. 
  //This is based on their extended notion of things being well-separated that uses mutual reachability terminology as well.
  hdbscanInternal::nodeCD(tree, coreDist, cdMin, cdMax, tree, S.data());

  //  for(int i = 0; i < 20; i++) {
  //       std::cout << "i:" <<  i << ", " << cdMin[i] << ", " << cdMax[i] << std::endl;
  // }


  floatType rhoLo = -0.1;
  floatType beta = 2;
  size_t numEdges = 0;

  floatType wspdTime = 0;
  floatType kruskalTime = 0;
  floatType markTime = 0;
  edgeUnionFind<long> UF(S.size());

  t0.stop();

  //This is the (parallel) Memo-GFK (GeoFilterKruskal in this while loop)
  while (UF.numEdge() < S.size() - 1) {
    t0.start();

    floatType rhoHi;
    auto bccps = filterWspdParallel<treeNode>(beta, rhoLo, rhoHi, tree, &UF,
					   coreDist, cdMin, cdMax, S.data());
    wspdTime += t0.get_next();

    cout << "---" << endl;
    cout << " beta = " << beta << endl;
    cout << " rho = " << rhoLo << " -- " << rhoHi << endl;

    numEdges += bccps.size();

    if (bccps.size() <= 0) {
      beta *= 2;
      rhoLo = rhoHi;
      continue;}

    cout << " edges = " << bccps.size() << endl;

    struct wEdge {
      size_t u,v;
      floatType weight;
    };

    auto base = S.data();
    sequence<wEdge> edges = tabulate(bccps.size(), [&](size_t i) {
        auto bcp = bccps[i];
        wEdge e;
        e.u = get<0>(bcp) - base;
        e.v = get<1>(bcp) - base;
        e.weight = get<2>(bcp);
        return e;
      });
    batchKruskal(edges, S.size(), UF);
    cout << " mst-edges = " << UF.numEdge() << endl;
    kruskalTime += t0.get_next();
    mark<treeNode, pointType, edgeUnionFind<long>>(tree, &UF, S.data());
    markTime += t0.stop();

    beta *= 2;
    rhoLo = rhoHi;
  }
  // floatType sum = 0;
  // auto E = UF.getEdge();
  // for (auto e: E) {
  //   floatType w = S[e.u].dist(S[e.v]);
  //   w = max(w, coreDist[e.u]);
  //   w = max(w, coreDist[e.v]);
  //   sum += w;
  // }
  // cout << "edge-sum = " << sum << endl;

  cout << "wspd-time = " << wspdTime << endl;
  cout << "kruskal-time = " << kruskalTime << endl;
  cout << "mark-time = " << markTime << endl;
  return UF.getEdge();
}

