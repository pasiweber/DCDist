#ifndef DC_DIST_HPP
#define DC_DIST_HPP

#include <vector>
#include <string>
#include <key_structs.hpp>

#include <mlpack/core.hpp>

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include <../parallel_hdbscan/src/kdTree.h>
#include <../parallel_hdbscan/src/kdTreeKnn.h>
#include <../parallel_hdbscan/include/hdbscan/point.h>
#include <../parallel_hdbscan/include/hdbscan/hdbscan.h>
#include <../parallel_hdbscan/include/hdbscan/edge.h>

#include <../parallel_hdbscan/src/kdTreeArma.h>
#include <../parallel_hdbscan/src/kdTreeKnnArma.h>
#include <../parallel_hdbscan/include/hdbscan/armapoint.h>



//Nodes for the dc-tree
/*
    Implementation structure:
    1. Get k'th nearest neighbor for each point
    1b. Compute mutual reachability graph (explicit or not?)
    2. Compute the MST over the complete mutual reachability graph
    3. Compute the dc-tree from the MST


    TODO(S): 
    1. Figure out what specific kd-tree implementation to use. Potentially also a ball-tree implementation. 
    2. Figure out if this specific implementation has a knn implemented already. If not, this needs to be implemented over the tree. 
        a) Does not seem too bad - maintain bounded priority queue in traversal over the tree that "prunes" non-useful regions of the tree.
    3. Adapt the dualtreetraverser dualtree boruvka - it needs to be able to take as input the mutual reachability distance specifically.
    4. Figure out how to implement the dc-tree over the output of MST. (What does the MST output format look like?)
    5. Potentially incorporate a base implementation of the dc-tree as well. In general, how should the code be strucured in terms of options as to how to construct it?

    Game plan:
    1. Download mlpack for the project
    2. Include it for the cmake
    3. Make a mututal reachability class, this is both given as input parameter to the KDTree and the DualTreeBoruvka which computes the MST
        3b. In this we need to compute the KNN, and store the core distances within the mututal reachability class. This means that computing the distance will take constant time. 
        3c. The class then needs to satisfy the metric policy and can then in theory just be plugged in and we are good to go..?
    4. Figure out the format of the MST -> whatever it is, the dc-tree should be computed from it - classical sorting the edges and constructing the tree.
        4b. Should most likely be done using a union find algorithm.

    https://github.com/junjiedong/KDTree/blob/master/src/KDTree.h
    https://github.com/cdalitz/kdtree-cpp/blob/master/demo_kdtree.cpp
    
*/





void printSubtree(const std::string &prefix, const Node& tree);
void printTree(const Node& tree);


void printSubtree2(const std::string &prefix, const Node& tree);
void printTree2(const Node& tree);

Node* construct_dc_tree(const std::vector<std::vector<double>> &points);


void swap(double *const a, double *const b);

unsigned long long partition(std::vector<double> &arr, const unsigned long long low, const unsigned long long high);

double quickSelect(std::vector<double> &arr, const unsigned long long low, const unsigned long long high, const int k);


std::vector<double> compute_cdists(arma::mat &data, size_t k, std::string mode);

std::vector<double> extract_cdists(arma::mat distances, size_t k);

std::vector<double> naive_cdists_efficient(arma::mat &data, size_t k);
std::vector<double> parallel_cdists(arma::mat &data, size_t k);
std::vector<double> parallel_cdists2(arma::mat &data, size_t k);
std::vector<double> parallel_cdists20(arma::mat &data, size_t k);



template<const int dim>
parlay::sequence<pargeo::point<dim>> convertArmaMatToParlayPoints(const arma::mat& mat);


std::vector<double> parallel_cdists_arma(arma::mat &data, size_t k);
parlay::sequence<pargeo::point2> convertArmaMatToParlayPoints2(const arma::mat& mat);

#endif