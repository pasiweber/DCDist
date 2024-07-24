#ifndef K_CENTROIDS2_HPP
#define K_CENTROIDS2_HPP

#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>

#include <dc_dist.hpp>
#include <lp_objective.hpp>





/*
Wrapper class for the tree itself to make fast assignments possible. I want to template it with the objective function somehow. 

*/
template <typename CostFunction = KMedian> 
class KCentroidsTree
{

    //Fields
   private:
    std::vector<int> index_order; //Contains the list of point ids in leaves sorted by left to right position in the tree
    Node *curr_k_solution; //Internal tree nodes representing the current k solution
    int curr_k; //Used for quick increasing / decreasing from a current solution
   public:
    Node *tree; //Should be a pointer rather than a reference due to c++ nature.

    //Methods


    //Constructors //This should create the tree when the class is constructed
    KCentroidsTree(Node& root) : tree(create_hierarchy(root)), curr_k(1){
        int n = root.size;
        setup_quick_clusters(n, this->tree);
        this->curr_k_solution = this->tree;

    }

    Node* get_tree(){
        return this->tree;
    }
    /*
    TODO:
    I want to have a dummy node that points to all the nodes in current k solution. This is easily storable and makes changing the k solution super easy.
    So first get all the nodes corresponding to the k solution.
    For the current k solution I will also end up making the dummy node a dummy tree instead to handle the case where multiple nodes are essentially merged in the tree. Should not be any worse though.
    All these internal dummy nodes will have their k as -1, so the algorithm should still work just fine.
    Also this new dummy node; remember to delete it.
    */
    std::vector<int> get_k_solution(int k){ //This function will always do it top down for now
        Node *curr_solution = this->tree; //Remember to delete this node when done.
        
        int n = this->tree->size;
        std::cout << "index order:" << std::endl;
        printLabels(index_order);

        if(k>=n || k < 1){
            this->curr_k = n;
            this->curr_k_solution = this->tree;
            std::cout << "easy return... "<< std::endl;
            return index_order; //Just return each id as a unique point as default for any invalid parameter for now.
        } else {
            std::cout << "advanced return" << std::endl;
            Node *solution_holder = new Node{nullptr, 0.0};
            solution_holder->k = -1; //This has a default k value that should never be part of a solution
            printTree2(*(this->tree));
            get_k_solution_helper(k, curr_solution, solution_holder);

            std::vector<int> res(n);
            extract_labels(res, solution_holder);
            
            this->curr_k_solution = solution_holder;
            std::cout << "Final result:" << std::endl;
            printLabels(res);
            return res;
        }
    }

    /*
    The generated tree structure should never get more than two deep
    */
    void get_k_solution_helper(int k, Node *tree, Node *new_solution){
        if(tree->size == 1){
            std::cout << "pushing leaf " << tree->id << std::endl;
            new_solution->children.push_back(tree); // We might end up with leaves part of the solution.
        }
        int min_child_k = std::numeric_limits<int>::max();
        std::cout << "startval" << min_child_k << std::endl;
        int max_child_k = tree->k;
        for(Node *child : tree->children){
            if(max_child_k < child->k){
                max_child_k = child->k;
            }
            if(child->k < min_child_k){
                min_child_k = child->k;
            }
        } 
        std::cout << "min, max" << min_child_k <<"," << max_child_k << std::endl;

        if(max_child_k <= k){ //Recurse if no cut detected
            for(Node *child : tree->children){
                    get_k_solution_helper(k, child, new_solution);
            }
        } else{
            if(min_child_k > k){ //Simple cut
                new_solution->children.push_back(tree);
                std::cout << "simple pushing c, id"  << tree->cost << ", " << tree->id << std::endl;
                
            } else{ //If some children have larger k and some have smaller, we need to merge things.
                //First create the merge node:
                Node* merge_node = new Node{nullptr, 0.0};
                merge_node->k = -1;
                merge_node->size = -1; //Just a "flag" value
                for(Node *child : tree->children){
                    if(child->is_orig_cluster || child->k > k){
                        merge_node->children.push_back(child);
                        std::cout << "merging c,id, k" << child->cost << ", " << child->id<< ", "<< child->k <<std::endl;
                    } else{
                        new_solution->children.push_back(child);
                        std::cout << "pushing c, id"  << child->cost << ", " << child->id << std::endl;
                    }
                   new_solution->children.push_back(merge_node);
                }
            }
        }
    }



    void extract_labels(std::vector<int> &res, Node *solution){
        int label = 0;
        for(Node *child : solution->children){
            if(child->size != -1){
                for(int i = child->low; i <= child->high; i++){
                    int id = index_order[i];
                    res[id-1] = label;
                }
            } else {
                for(Node *child_2 : child->children){
                    for(int i = child->low; i <= child->high; i++){
                        int id = index_order[i];
                        res[id-1] = label;
                    }
                }
            }
            label++;
        }
    }





    std::vector<double> inc_k_solution();

    std::vector<double> dec_k_solution();



    private:
        /*
        Creates the hierarchy tree, and annotates each node with the k instance it was created at.
        */
        Node* create_hierarchy(Node& root){
            std::vector<Annotation*> annotations = annotate_tree(root);
            std::sort(annotations.begin(), annotations.end(), compareByCost);


            //Create the tree here using the pointers in the annotations.
            Node* root_pointer = new Node{nullptr, 0.0, annotations[0]->center};
            root_pointer->k = 1; //Annotation for getting k solutions quickly
            annotations[0]->tree_node = root_pointer; //Unroll first iteration of loop to avoid an if/else in loop.

            Annotation* curr_anno; 
            Annotation* parent_anno;
            Node* parent_node;
            std::vector<Annotation*> leaves_to_add;

            //First loop/pass adding the main structure of the tree
            for(int i = 1; i < annotations.size(); i++){
                curr_anno = annotations[i];
                Node* new_node = new Node{nullptr, 0.0, curr_anno->center};
                new_node->k = i+1;//Annotation for getting k solutions quickly
                //Fix the linkage
                parent_anno = curr_anno->parent;
                parent_node = parent_anno->tree_node;
                double cost = parent_node->cost;
                if(parent_anno->has_leaf){ //This parent center will not have a corresponding leaf as we will now be adding stuff below it
                    parent_anno->has_leaf = false; //We will now be adding things below this annotation/node -> this will not become a leaf in this pass
                    leaves_to_add.push_back(parent_anno); //TODO: We need to update the k that should be inserted for this leaf somehow
                }

                if(cost !=0 && cost != curr_anno->cost_decrease){ //If cost is not 0 and not the current annos cost decrease we need to make a new internal node for the parent. 
                    //Add new node, replace in annotation
                    Node* new_parent = new Node{parent_node, curr_anno->cost_decrease, parent_node->id};
                    int new_k = parent_anno->k;
                    new_parent->k = new_k; //Annotation for getting k solutions quickly
                    new_parent->is_orig_cluster = true; //Used to distinguish nodes when getting clusterings out for each k
                    parent_anno->k = i+1; //Annotation for getting k solutions quickly

                    new_parent->children.push_back(new_node);
                    new_node->parent = new_parent;
                    parent_anno->tree_node = new_parent;
                    parent_node->children.push_back(new_parent); //Link new parent to old parent
                } else{
                    //Link to current node and update cost in parent(might just rewrite the same value but who cares (not me))
                    parent_node->cost = curr_anno->cost_decrease;
                    parent_node->children.push_back(new_node);
                    new_node->parent = parent_node;
                    if(cost == 0){
                        parent_anno->k = i+1; //Annotation for getting k solutions quickly (when we eventually do the parent split into a new node, this is then the k used there)
                    }
                }

                curr_anno->tree_node = new_node;
            }

            //Second loop/pass adding final leaves to the tree
            for(int i = 0; i < leaves_to_add.size(); i++){ 
                curr_anno = leaves_to_add[i];
                if(!(curr_anno->has_leaf)){ //TODO: Remove this if statement - not needed I think
                    Node* node = curr_anno->tree_node;
                    Node* new_leaf = new Node{node, 0.0, curr_anno->center};
                    new_leaf->k = curr_anno->k;
                    new_leaf->is_orig_cluster = true; //Used to distinguish nodes when getting clusterings out for each k
                    node->children.push_back(new_leaf);
                }
            }
            assign_sizes(root_pointer); //Add number of leaves in subtree to each node
            delete_annotations(annotations); //Free the annotations from memory

            return root_pointer;
        }


        std::vector<Annotation*> annotate_tree(Node &root){
            std::vector<Annotation*> arr;
            arr.resize(root.size);
            annotate_tree_inner(root, arr);
            return arr;
        }

        /*
        Inner tree traversal function that computes the actual annotations and appends them to the arr.
        Currently I assume the tree stores sizes - otherwise it makes this code a lot uglier - we really need to store the sizes in advance, otherwise we need to do two traversals of the list of children.
        TODO: Should I somehow already keep track of the corresponding clusters?
        */
        std::pair<double, int> annotate_tree_inner(Node &tree, std::vector<Annotation*> &arr){
            if((tree.children).size() == 0){
                Annotation* anno =  new Annotation{CostFunction::Evaluate((tree.parent)->cost), tree.id}; //Removed children as a list in annotations - not needed as previously thought
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
                    std::pair<double, int> sub_result = annotate_tree_inner(*child, arr);
                    curr_sub_cost = sub_result.first;
                    curr_center = sub_result.second;
                    int child_size = child->size;
                    double curr_cost = curr_sub_cost + CostFunction::Evaluate(tree.cost) * (size-child_size);
                    if(curr_cost <= best_cost){
                        best_cost = curr_cost;
                        best_center = curr_center;
                    }

                    centers.push_back(curr_center);
                }
                double cost_decrease = CostFunction::Evaluate(parent_cost) * size - best_cost;
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



        void setup_quick_clusters(int n, Node* centroid_tree){
            std::vector<int> arr;
            arr.resize(n);
            quick_clusters_helper(centroid_tree, arr, 0);
            printLabels(arr);
            this->index_order = arr;
        }

        /*
        Insert the leaf ids in an array based on their ordering in the tree. This creates a continuous array segment for each internal node of its leaves.
        Annotate each inner node with this segment. Also annotates the inner nodes with the size.
        */
        int quick_clusters_helper(Node* tree, std::vector<int> &arr, int ctr){
            if(tree->children.size() == 0){ //leaf
                tree->size = 1;
                tree->low = ctr;
                tree->high = ctr;
                arr[ctr] = tree->id;
                return ctr+1; //Could also decide to update the ID here to be ctr instead.
            } else{ //
                tree->low = ctr; //This ctr index will be used by leftmost child
                for(Node *child : tree->children){
                    ctr = quick_clusters_helper(child, arr, ctr);
                }
                tree->high = ctr-1; //ctr is what was returned by rightmost child (but it used ctr-1 as its own index)
                tree->size = tree->high - tree->low + 1; //+1 since we include both ends of spectrum
                //std::cout << "cost " << tree->cost << ", low, high:" << tree->low<<","<<tree->high << ", size " << tree->size << std::endl;
                return ctr;
            }
        }


















        void printLabels(std::vector<int> labels){
            std::cout << "Labels:" << std::endl;
            std::cout << "[";
            for(int label : labels){
                std::cout << label << ", ";
            }

            std::cout << "]" << std::endl;

        }


};


#endif