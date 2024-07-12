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
    double cost;  // For DCTree this is dc-distance, for... 
    int id; // The id of a potential leaf node
    std::vector<struct Node*> children;
    int size;
} Node;


//Cost-decrease annotation structs
typedef struct Annotation {
    double cost_decrease;
    int center;
    Annotation* parent;
    Node* tree_node = nullptr;
    bool has_leaf = true; //This is set to true until a node sets it to false. When we are "done" with the tree, we do another run over the annotations, adding leaves of those set to false, as a child of the node the annotation is currently pointing to
} Annotation;


bool is_a_root(const Node &tree);


bool compareByCost(const Annotation* anno1, const Annotation* anno2);
void print_annotations(std::vector<Annotation*> annotations);






template <typename CostFunction>
std::vector<Annotation*> annotate_tree_old(const Node &root, CostFunction f){
    std::vector<Annotation*> arr;
    int n = 2*root.size; 
    arr.reserve(n);
    annotate_tree_inner_old(root, f, arr);
    return arr;
}


/*
Inner tree traversal function that computes the actual annotations and appends them to the arr.
Currently I assume the tree stores sizes - otherwise it makes this code a lot uglier - we really need to store the sizes in advance, otherwise we need to do two traversals of the list of children.
TODO: Should I somehow already keep track of the corresponding clusters?
*/
template <typename CostFunction>
std::pair<double, int> annotate_tree_inner_old(const Node &tree, CostFunction f, std::vector<Annotation*> &arr){
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
            std::pair<double, int> sub_result = annotate_tree_inner_old(*child, f, arr);
            curr_sub_cost = sub_result.first;
            curr_center = sub_result.second;
            int child_size = child->size;
            double curr_cost = curr_sub_cost + f(tree.cost) * (size-child_size);
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




int split_detector(std::vector<Annotation*> annotations, int i);

void assign_sizes(Node* root);
int assign_size_helper(Node* tree);
void delete_annotations(std::vector<Annotation*> annotations);

template <typename CostFunction>
Node* create_hierarchy(Node& root, CostFunction f){
    std::vector<Annotation*> annotations = annotate_tree(root, f);
    std::sort(annotations.begin(), annotations.end(), compareByCost);
    print_annotations(annotations);


    //Create the tree here using the pointers in the annotations.
    /*
    TODO: I need a split detector - how many split at once for the same parent center.
    The code here also needs to then create the extra new nodes corresponding to same old center when split occurs. 
    I can potentially start by making a nodes cost the id, and then set it to a cost when we detect it is not a leaf.
    TODO Assign sizes to internal nodes when done

    TODO: Change code to use references instead of pointers for annotations.

    */
    
    Node* root_pointer = new Node{nullptr, 0.0, annotations[0]->center};
    annotations[0]->tree_node = root_pointer; //Unroll first iteration of loop to avoid an if/else in loop.

    Annotation* curr_anno; 
    Annotation* parent_anno;
    std::vector<Annotation*> leaves_to_add;

    //First loop/pass adding the main structure of the tree
    for(int i = 1; i < annotations.size(); i++){
        curr_anno = annotations[i];
        Node* new_node = new Node{nullptr, 0.0, curr_anno->center};
        //Fix the linkage
        parent_anno = curr_anno->parent;
        Node* p_node = parent_anno->tree_node;
        double cost = p_node->cost;
        if(parent_anno->has_leaf){ //This parent center will not have a corresponding leaf as we will now be adding stuff below it
            parent_anno->has_leaf = false; //We will now be adding things below this annotation/node -> this will not become a leaf in this pass
            leaves_to_add.push_back(parent_anno);
        }

        if(cost !=0 && cost != curr_anno->cost_decrease){ 
            //Add new node, replace in annotation
            Node* new_parent = new Node{p_node, curr_anno->cost_decrease, p_node->id};
            new_parent->children.push_back(new_node);
            new_node->parent = new_parent;
            parent_anno->tree_node = new_parent;
            p_node->children.push_back(new_parent); //Link new parent to old parent
        } else{
            //Link to current node and update cost in parent(might just rewrite the same value but who cares (not me))
            p_node->cost = curr_anno->cost_decrease;
            p_node->children.push_back(new_node);
            new_node->parent = p_node;
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
std::vector<Annotation*> annotate_tree(const Node &root, CostFunction f){
    std::cout << "annotating..." << std::endl;
    std::vector<Annotation*> arr;
    std::cout << "size:" << root.size << std::endl;
    arr.resize(root.size);
    std::cout << "calling helper" << std::endl;
    annotate_tree_inner(root, f, arr);
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
        Annotation* anno =  new Annotation{f((tree.parent)->cost), tree.id}; //Removed children as a list in annotations - not needed as previously thought
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

        for(const Node* child : tree.children){ //Run through the children of the current node and find the lowest
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
        std::pair<double, int> result = {best_cost, best_center};
        return result;
    }
}




void delete_tree(Node* tree);


#endif