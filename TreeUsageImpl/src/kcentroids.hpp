#ifndef K_CENTROIDS_HPP
#define K_CENTROIDS_HPP

#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>


//Nodes for the dc-tree
/*
Thoughts: An inherent weakness is that we seemingly cannot predefine the size of the children vector
*/
typedef struct Node {
    struct Node* parent;
    std::vector<struct Node*> children;
    double cost;  // For DCTree this is dc-distance, for... 
    int id; // The id of a potential leaf node
    int size;
} Node;

//Cost-decrease annotation structs
typedef struct Annotation {
    double cost_decrease;
    int center;
} Annotation;




template <typename CostFunction>
std::vector<Annotation*> annotate_tree(const Node &root, CostFunction f){
    std::vector<Annotation*> arr;
    int n = 2*root.size; 
    arr.reserve(n);
    annotate_tree_inner(root, f, arr);
    return arr;
}

bool is_a_root(const Node &tree);

/*
Inner tree traversal function that computes the actual annotations and appends them to the arr.
Currently I assume the tree stores sizes - otherwise it makes this code a lot uglier - we really need to store the sizes in advance, otherwise we need to do two traversals of the list of children.
TODO: Should I somehow already keep track of the corresponding clusters?
*/
template <typename CostFunction>
std::pair<double, int> annotate_tree_inner(const Node &tree, CostFunction f, std::vector<Annotation*> &arr){
    if((tree.children).size() == 0){
        Annotation* anno =  new Annotation{f((tree.parent)->cost), tree.id};
        arr.push_back(anno);
        std::pair<double, int> result = {0.0, tree.id};
        return result;
    } else{
        int size = tree.size;
        double parent_cost = 0;
        if(is_a_root(tree)){
            parent_cost = tree.cost; //If we are in the root we use the distance itself. This will still result in the annotation having the tied highest cost decrease as desired. TODO: Check that this is true.
        }
        else {
            parent_cost = (tree.parent)->cost;
        }
        double best_cost = std::numeric_limits<double>::infinity();
        int best_center = -1;
        double curr_sub_cost = 0.0;
        int curr_center = -1;
        for(const Node* child : tree.children){ //Run through the children of the current node and find the lowest
            std::pair<double, int> sub_result = annotate_tree_inner(*child, f, arr);
            curr_sub_cost = sub_result.first;
            curr_center = sub_result.second;
            int child_size = child->size;
            double curr_cost = curr_sub_cost + f(parent_cost) * (size-child_size);
            if(curr_cost <= best_cost){
                best_cost = curr_cost;
                best_center = curr_center;
            }
        }
        double cost_decrease = f(parent_cost) * size - best_cost;
        Annotation* anno = new Annotation{cost_decrease, best_center};
        arr.push_back(anno);
        std::pair<double, int> result = {best_cost, best_center};
        return result;
    }
}

bool compareByCost(const Annotation* anno1, const Annotation* anno2);
void print_annotations(std::vector<Annotation*> annotations);

template <typename CostFunction>
Node* create_hierarchy(Node& root, CostFunction f){
    std::vector<Annotation*> cost_decreases = annotate_tree(root, f);
    std::cout << "before:" << std::endl;
    print_annotations(cost_decreases);
    std::sort(cost_decreases.begin(), cost_decreases.end(), compareByCost);
    std::cout << "after:" << std::endl;
    print_annotations(cost_decreases);

    return &root;
}
/*
Issue: The tree has no balance guarantees - this makes it (at least as I'm thinking now) difficult to make guarantees on complexity better than O(n^2) or space of O(n^2)
TODO: Fix that

TODO2: Add ability to get clustering for specific k - i.e. implement the extraction of the clustering given the centers.
*/

#endif