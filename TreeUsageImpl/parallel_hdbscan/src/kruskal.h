// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2019 Guy Blelloch and the PBBS team
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

/*
  Added a similar version of Kruskal's algorithm to process batches of edges
  for the pargeo project
 */


/*
  MY THOUGHTS:
  It seems that this file does not use any pargeo functionality, it simply defines everything within the pargeo namespace.
  Potentially, the only work that needs to be done is removing the namespaces.

*/




#pragma once

#include <iostream>
#include <limits.h>
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "speculativeFor.h"
#include "unionFind.h"

// **************************************************************
//    PARALLEL VERSION OF KRUSKAL'S ALGORITHM
// **************************************************************

namespace pargeo {

  using namespace std;

  namespace kruskalInternal {
    using vertexId = long;
    using edgeId = long;
    using edgeWeight = double;

    // need to tag each edge with an index so we can keep track of
    // which edges are added to the MST
    struct indexedEdge {
      vertexId u; vertexId v; edgeId id; edgeWeight w;
    indexedEdge(vertexId u, vertexId v, edgeId id, edgeWeight w)
      : u(u), v(v), id(id), w(w){}
      indexedEdge() {};
    };

    using reservation = reservation<edgeId>;

    struct UnionFindStep {
      parlay::sequence<indexedEdge> &E;
      parlay::sequence<reservation> &R;
      unionFind<vertexId> &UF;
      parlay::sequence<bool> &inST;

      UnionFindStep(parlay::sequence<indexedEdge> &E, unionFind<vertexId> &UF,
        parlay::sequence<reservation> &R, parlay::sequence<bool> &inST) : E(E), R(R), UF(UF), inST(inST) {}

      bool reserve(edgeId i) {
        vertexId u = E[i].u = UF.find(E[i].u);
        vertexId v = E[i].v = UF.find(E[i].v);
        if (u != v) {
          R[v].reserve(i);
          R[u].reserve(i);
          return true;
        } else return false;
      }

      bool commit(edgeId i) {
        vertexId u = E[i].u;
        vertexId v = E[i].v;
        if (R[v].check(i)) {
          R[u].checkReset(i);
          UF.link(v, u);
          inST[E[i].id] = true;
          return true;}
        else if (R[u].check(i)) {
          UF.link(u, v);
          inST[E[i].id] = true;
          return true; }
        else return false;
      }
    };

    // Union find step without MST mapping
    struct EdgeUnionFindStep {
      parlay::sequence<indexedEdge> &E;
      parlay::sequence<indexedEdge> EReal;
      parlay::sequence<reservation> &R;
      edgeUnionFind<vertexId> &UF;
      EdgeUnionFindStep(parlay::sequence<indexedEdge> &E,
			edgeUnionFind<vertexId> &UF,
			parlay::sequence<reservation> &R) : E(E), R(R), UF(UF) {
        EReal = parlay::sequence<indexedEdge>(E);
      }

      bool reserve(edgeId i) {
        vertexId u = E[i].u = UF.find(E[i].u);
        vertexId v = E[i].v = UF.find(E[i].v);
        if (u != v) {
          R[v].reserve(i);
          R[u].reserve(i);
          return true;
        } else return false;
      }

      bool commit(edgeId i) {
        vertexId u = E[i].u;
        vertexId v = E[i].v;
        vertexId uReal = EReal[i].u;
        vertexId vReal = EReal[i].v;
        if (R[v].check(i)) {
          R[u].checkReset(i);
          UF.link(v, u, vReal, uReal, EReal[i].w);
          return true;}
        else if (R[u].check(i)) {
          UF.link(u, v, uReal, vReal, EReal[i].w);
          return true; }
        else return false;
      }
    }; //End EdgeUnionFindStep Struct

  } // End namespace kruskalInternal

  template <typename edgeT>
  parlay::sequence<long> kruskal(parlay::sequence<edgeT> &E, size_t n) {
    using namespace kruskalInternal;

    //timer t("mst", true);
    /* size_t m = E.m; */
    size_t m = E.size();
    /* size_t n = E.n; */
    size_t k = min<size_t>(5 * n / 4, m);

    // equal edge weights will prioritize the earliest one
    auto edgeLess = [&] (indexedEdge a, indexedEdge b) {
      return (a.w < b.w) || ((a.w == b.w) && (a.id < b.id));};

    // tag each edge with an index
    auto IW = parlay::delayed_seq<indexedEdge>(m, [&] (size_t i) {
	    return indexedEdge(E[i].u, E[i].v, i, E[i].weight);
    });

    auto IW1 = parlay::sort(IW, edgeLess);
    //t.next("sort edges");

    parlay::sequence<bool> mstFlags(m, false);
    unionFind<vertexId> UF(n);
    parlay::sequence<kruskalInternal::reservation> R(n);
    UnionFindStep UFStep1(IW1, UF, R,  mstFlags);
    speculative_for<vertexId>(UFStep1, 0, IW1.size(), 20, false);
    //t.next("union find loop");

    parlay::sequence<edgeId> mst = parlay::pack_index<edgeId>(mstFlags);
    //t.next("pack out results");

    return mst;
  }

  template <typename edgeT>
  void batchKruskal(parlay::sequence<edgeT> &E,size_t n, edgeUnionFind<long>& UF) {
    using namespace kruskalInternal;

    size_t m = E.size();
    size_t k = min<size_t>(5 * n / 4, m);

    // equal edge weights will prioritize the earliest one
    auto edgeLess = [&] (indexedEdge a, indexedEdge b) {
      return (a.w < b.w) || ((a.w == b.w) && (a.id < b.id));
    };

    // tag each edge with an index
    auto IW = parlay::delayed_seq<indexedEdge>(m, [&] (size_t i) {
	    return indexedEdge(E[i].u, E[i].v, i, E[i].weight);
    });

    auto IW1 = parlay::sort(IW, edgeLess);

    parlay::sequence<kruskalInternal::reservation> R(n);
    EdgeUnionFindStep UFStep1(IW1, UF, R);
    speculative_for<vertexId>(UFStep1, 0, IW1.size(), 20, false);
  }
} // End namespace pargeo
