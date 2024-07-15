#ifndef DC_DIST_HPP
#define DC_DIST_HPP

#include <vector>
#include <string>



//Nodes for the dc-tree
/*
Thoughts: An inherent weakness is that we seemingly cannot predefine the size of the children vector
*/
typedef struct Node {
    struct Node* parent;
    double cost;  // For DCTree this is dc-distance, for... 
    int id; // The id of a potential leaf node (this can be used to assign points -> we just make the id in internal nodes the optimal center, we then cap by k in the iterations -> this will be separate from the tree hierarchy construction - O(n) for each k to output the solution.)
    std::vector<struct Node*> children;
    int size;
    bool is_cluster = false; //Used to extract clusters for optimization algorithm over the tree
} Node;


void printSubtree(const std::string &prefix, const Node& tree);
void printTree(const Node& tree);

Node* construct_dc_tree(const std::vector<std::vector<double>> &points);

#endif