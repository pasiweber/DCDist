#ifndef K_CENTROIDS_HPP
#define K_CENTROIDS_HPP

#include <vector>
#include <limits>
#include <iostream>



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
    std::cout << "here1" << std::endl;
    std::vector<Annotation*> arr;
    int n = 100; //TODO: We need to extract the size to preallocate
    arr.reserve(n);
    annotate_tree_inner(root, f, arr);
    std::cout << "here2" << std::endl;

    return arr;
}



/*
Inner tree traversal function that computes the actual annotations and appends them to the arr.
Currently I assume the tree stores sizes - otherwise it makes this code a lot uglier - we really need to store the sizes in advance, otherwise we need to do two traversals of the list of children.
TODO: Should I somehow already keep track of the corresponding clusters?
*/
template <typename CostFunction>
std::pair<double, int> annotate_tree_inner(const Node &tree, CostFunction f, std::vector<Annotation*> &arr){
    if((tree.children).size() == 0){
        std::cout << "here3" << std::endl;

        Annotation* anno =  new Annotation{(tree.parent)->cost, tree.id};
        arr.push_back(anno);
        std::pair<double, int> result = {0.0, tree.id};
        return result;
    } else{
        std::cout << "here4" << std::endl;

        int size = tree.size; //Handle these two - they segfault due to the tree structure.
        double parent_cost = (tree.parent)->cost;
        std::cout << "here4.1" << std::endl;
        double best_cost = std::numeric_limits<double>::infinity();
        int best_center = -1;
        double curr_sub_cost = 0.0;
        int curr_center = -1;
        std::cout << "here5" << std::endl;
        for(const Node* child : tree.children){ //Run through the children of the current node and find the lowest
            std::cout << "here6" << std::endl;
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
        std::cout << "here7" << std::endl;
        double cost_decrease = f(parent_cost) * size - best_cost;
        Annotation* anno = new Annotation{cost_decrease, best_center};
        arr.push_back(anno);
        std::pair<double, int> result = {best_cost, best_center};
        std::cout << "here8" << std::endl;
        return result;
    }
}
    
#endif