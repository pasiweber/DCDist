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

#include "../include/hdbscan/point.h"
#include "../include/hdbscan/hdbscan.h"

using namespace std;
using namespace parlay;
using namespace pargeo;
using namespace pargeo::hdbscanInternal;

template<int dim>
parlay::sequence<pargeo::wghEdge> pargeo::hdbscan(parlay::sequence<pargeo::point<dim>> &S, size_t minPts) {
    using pointType = point<dim>; //pargeo point
    using treeNode = kdNode<dim, point<dim>>;
    using floatType = typename pointType::floatType;

    if (S.size() < 2) {
        throw std::runtime_error("need more than 2 points");
    }

    timer t0;
    t0.start();

    //Build the kd-tree for k nearest neighbors and find k-nn
    treeNode* tree = buildKdt<dim, point<dim>>(S, true, true);
    cout << "build-tree-time = " << t0.get_next() << endl;

    sequence<size_t> nns = kdTreeKnn<dim, pointType>(S, minPts, tree, true); 
    sequence<floatType> coreDist = sequence<floatType>(S.size());
    parallel_for (0, S.size(), [&](intT i) {
        coreDist[i] = S[nns[i*minPts + minPts-1]].dist(S[i]);
        // if(i % 20 == 0){
        //     std::cout << "i:" <<  i << ", " << coreDist[i] << std::endl;
        // }
    });
    cout << "core-dist-time = " << t0.get_next() << endl;

    //Two arrays that serve as min and max for each internal node in the kd-tree - they get initialized here to max and min default values
    sequence<floatType> cdMin = sequence<floatType>(tree->size() * 2);
    sequence<floatType> cdMax = sequence<floatType>(tree->size() * 2);
    parallel_for(0, tree->size()*2, [&](intT i) {
        cdMin[i] = std::numeric_limits<floatType>::max();
        cdMax[i] = std::numeric_limits<floatType>::lowest();
    });
    hdbscanInternal::nodeCD(tree, coreDist, cdMin, cdMax, tree, S.data()); //Set cdmin and cdmax to actual values for each node

    // for(int i = 0; i < 20; i++) {
    //     std::cout << "i:" <<  i << ", " << cdMin[i] << ", " << cdMax[i] << std::endl;
    // }

    floatType rhoLo = -0.1;
    floatType beta = 2; //Starts as 2, doubles each round. We partition pairs into cardinality < beta and > beta.
    size_t numEdges = 0;

    floatType wspdTime = 0;
    floatType kruskalTime = 0;
    floatType markTime = 0;
    
    //UF is the union find structure
    edgeUnionFind<long> UF(S.size()); //Initialize S.size() elements all as separate roots
    t0.stop();

    //This while loop is the (parallel) Memo-GFK (Algorithm 3 from the paper)
    while (UF.numEdge() < S.size() - 1) { //The loop continues until we have as many edges as there are nodes
        t0.start();

        floatType rhoHi; // rhoHI is the upper bound on edges in S_l (those with cardinality lower than beta): minimum 𝑑(𝐴,𝐵) for all (𝐴,𝐵) ∈ 𝑆_u
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

        //Initialize the set of edges
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

    cout << "wspd-time = " << wspdTime << endl;
    cout << "kruskal-time = " << kruskalTime << endl;
    cout << "mark-time = " << markTime << endl;
    return UF.getEdge();
}

template sequence<wghEdge> pargeo::hdbscan<2>(sequence<point<2>> &, size_t);
template sequence<wghEdge> pargeo::hdbscan<200>(sequence<point<200>> &, size_t);
template sequence<wghEdge> pargeo::hdbscan<24>(sequence<point<24>> &, size_t);

