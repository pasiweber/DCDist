#ifndef K_CENTROIDS_HPP
#define K_CENTROIDS_HPP

#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>

#include <dc_dist.hpp>




//Cost-decrease annotation structs
typedef struct Annotation {
    double cost_decrease;
    int center;
    Annotation* parent;
    Node* tree_node = nullptr;
    bool has_leaf = true; //This is set to true until a node sets it to false. When we are "done" with the tree, we do another run over the annotations, adding leaves of those set to false, as a child of the node the annotation is currently pointing to
    Node* orig_node; //This is used to extract a specific k solution efficiently
} Annotation;


bool is_a_root(const Node &tree);


bool compareByCost(const Annotation* anno1, const Annotation* anno2);
void print_annotations(std::vector<Annotation*> annotations);





void assign_sizes(Node* root);
int assign_size_helper(Node* tree);
void delete_annotations(std::vector<Annotation*> &annotations);

template <typename CostFunction>
Node* create_hierarchy(Node& root, CostFunction f){
    std::vector<Annotation*> annotations = annotate_tree(root, f);
    std::sort(annotations.begin(), annotations.end(), compareByCost);
    //print_annotations(annotations);


    //Create the tree here using the pointers in the annotations.
    Node* root_pointer = new Node{nullptr, 0.0, annotations[0]->center};
    annotations[0]->tree_node = root_pointer; //Unroll first iteration of loop to avoid an if/else in loop.

    Annotation* curr_anno; 
    Annotation* parent_anno;
    Node* parent_node;
    std::vector<Annotation*> leaves_to_add;

    //First loop/pass adding the main structure of the tree
    for(int i = 1; i < annotations.size(); i++){
        curr_anno = annotations[i];
        Node* new_node = new Node{nullptr, 0.0, curr_anno->center};
        //Fix the linkage
        parent_anno = curr_anno->parent;
        parent_node = parent_anno->tree_node;
        double cost = parent_node->cost;
        if(parent_anno->has_leaf){ //This parent center will not have a corresponding leaf as we will now be adding stuff below it
            parent_anno->has_leaf = false; //We will now be adding things below this annotation/node -> this will not become a leaf in this pass
            leaves_to_add.push_back(parent_anno);
        }

        if(cost !=0 && cost != curr_anno->cost_decrease){ //If cost is not 0 and not the current annos cost decrease we need to make a new internal node for the parent. 
            //Add new node, replace in annotation
            Node* new_parent = new Node{parent_node, curr_anno->cost_decrease, parent_node->id};
            new_parent->children.push_back(new_node);
            new_node->parent = new_parent;
            parent_anno->tree_node = new_parent;
            parent_node->children.push_back(new_parent); //Link new parent to old parent
        } else{
            //Link to current node and update cost in parent(might just rewrite the same value but who cares (not me))
            parent_node->cost = curr_anno->cost_decrease;
            parent_node->children.push_back(new_node);
            new_node->parent = parent_node;
        }

        curr_anno->tree_node = new_node;
    }

    //Second loop/pass adding final leaves to the tree
    for(int i = 0; i < leaves_to_add.size(); i++){ 
        curr_anno = leaves_to_add[i];
        if(!(curr_anno->has_leaf)){
            Node* node = curr_anno->tree_node;
            Node* new_leaf = new Node{node, 0.0, curr_anno->center};
            node->children.push_back(new_leaf);
        }
    }
    assign_sizes(root_pointer); //Add number of leaves in subtree to each node
    delete_annotations(annotations); //Free the annotations from memory

    return root_pointer;
}




/*
Leaf annotations and inner node annotations
Should have pointers to parents and children
Each annotation should also have a pointer to the resulting node when it has been created.
The pointers to the parent should be updated to a pointer to the node as we keep going
Only annotations that have another center (with a different center) annotation as a parent should get child pointers and so on.
We keep the maximal annotations in a list and return a pointer + a center in the recursive call - in that way we can update in place on that location and only add children to one annotation.
*/
template <typename CostFunction>
std::vector<Annotation*> annotate_tree(Node &root, CostFunction f){
    std::vector<Annotation*> arr;
    arr.resize(root.size);
    annotate_tree_inner(root, f, arr);
    return arr;
}
/*
Inner tree traversal function that computes the actual annotations and appends them to the arr.
Currently I assume the tree stores sizes - otherwise it makes this code a lot uglier - we really need to store the sizes in advance, otherwise we need to do two traversals of the list of children.
TODO: Should I somehow already keep track of the corresponding clusters?
*/
template <typename CostFunction>
std::pair<double, int> annotate_tree_inner(Node &tree, CostFunction f, std::vector<Annotation*> &arr){
    if((tree.children).size() == 0){
        Annotation* anno =  new Annotation{f((tree.parent)->cost), tree.id}; //Removed children as a list in annotations - not needed as previously thought
        anno->orig_node = &tree;
        
        arr[tree.id-1] = anno; //The ids are 1-indexed
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
        std::vector<int> centers;

        for(Node* child : tree.children){ //Run through the children of the current node and find the lowest
            std::pair<double, int> sub_result = annotate_tree_inner(*child, f, arr);
            curr_sub_cost = sub_result.first;
            curr_center = sub_result.second;
            int child_size = child->size;
            double curr_cost = curr_sub_cost + f(tree.cost) * (size-child_size);
            if(curr_cost <= best_cost){
                best_cost = curr_cost;
                best_center = curr_center;
            }

            centers.push_back(curr_center);
        }
        double cost_decrease = f(parent_cost) * size - best_cost;
        //Update annotation list
        for(const int center : centers){ 
            if(center == best_center){
                continue;
            }
            arr[center-1]->parent = arr[best_center-1];
        }
        arr[best_center-1]->cost_decrease = cost_decrease;
        arr[best_center-1]->orig_node = &tree; //For extracting specific k solutions

        std::pair<double, int> result = {best_cost, best_center};
        return result;
    }
}

void annotate_tree_centers(std::vector<Annotation*> &annotations, int k);
std::vector<int> assign_points(Node &tree);
void assign_points_helper(Node &tree, std::vector<int> &arr, int curr_center);
void delete_tree(Node* tree); // https://stackoverflow.com/questions/4790564/finding-memory-leaks-in-a-c-application-with-visual-studio
void printLabels(std::vector<int> labels);
void cleanTree(Node &tree);

template <typename CostFunction>
std::vector<int> kcentroids(Node& root, CostFunction f, int k){
    std::vector<Annotation*> annotations = annotate_tree(root, f);
    std::sort(annotations.begin(), annotations.end(), compareByCost);

    annotate_tree_centers(annotations, k); //This puts the center annotation in id of internal nodes and cost of leaf nodes
    std::vector<int> res = assign_points(root);
    
    cleanTree(root); //Reset all internal node ids to -1 and all leaf node costs to 0
    printLabels(res);
    return res;
}
#endif