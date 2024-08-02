#ifndef KEY_STRUCTS_HPP
#define KEY_STRUCTS_HPP

#include <vector>

typedef struct Node {
    struct Node* parent;
    double cost;  // For DCTree this is dc-distance, for... 
    int id; // The id of a potential leaf node (this can be used to assign points -> we just make the id in internal nodes the optimal center, we then cap by k in the iterations -> this will be separate from the tree hierarchy construction - O(n) for each k to output the solution.)
    int k; // The k at which this node was created
    std::vector<struct Node*> children;

    //Structure quick access
    int size;
    int low; //Fast array indexing low index
    int high; //Fast array indexing high index
    //Kcentroids annotations
    bool is_orig_cluster; //Used to show that this node is the one corresponding to the parent that gets points taken from it. Basically, it has the same center as the parent node.

    //HDBSCAN annotations
    bool is_cluster = false; //Used to extract clusters for optimization (HDBSCAN / HCF) algorithm over the tree
} Node;



//Cost-decrease annotation structs
typedef struct Annotation {
    double cost_decrease;
    int center;
    Annotation* parent;
    Node* tree_node = nullptr;
    bool has_leaf = true; //This is set to true until a node sets it to false. When we are "done" with the tree, we do another run over the annotations, adding leaves of those set to false, as a child of the node the annotation is currently pointing to
    int k; //The next k for this annotation when it becomes a new node / smaller cluster again.
    Node* orig_node; //This is used to extract a specific k solution efficiently
} Annotation;



#endif