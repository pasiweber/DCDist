#ifndef DC_DIST_HPP
#define DC_DIST_HPP

#include <vector>
#include <string>


//mlpack stuff
#include <mlpack.hpp>


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




typedef struct Node {
    struct Node* parent;
    double cost;  // For DCTree this is dc-distance, for... 
    int id; // The id of a potential leaf node (this can be used to assign points -> we just make the id in internal nodes the optimal center, we then cap by k in the iterations -> this will be separate from the tree hierarchy construction - O(n) for each k to output the solution.)
    std::vector<struct Node*> children;

    //Structure quick access
    int size;
    int low; //Fast array indexing low index
    int high; //Fast array indexing high index
    //Kcentroids annotations
    int k; // The k at which this node was created
    bool is_orig_cluster; //Used to show that this node is the one corresponding to the parent that gets points taken from it. Basically, it has the same center as the parent node.

    //HDBSCAN annotations
    bool is_cluster = false; //Used to extract clusters for optimization (HDBSCAN / HCF) algorithm over the tree
} Node;


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
std::vector<double> naive_cdists_efficient2(arma::mat &data, size_t k);
#endif