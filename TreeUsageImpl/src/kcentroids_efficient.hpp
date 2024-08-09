#ifndef K_CENTROIDS_EFFICIENT_HPP
#define K_CENTROIDS_EFFICIENT_HPP

#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>

#include <dc_dist.hpp>
#include <lp_objective.hpp>

/*
    TODO: Remember to delete the old internal "fake" pointer nodes when we create a new k solution.
    Possible TODOs:
    1. Implement dec_k_solution
    2. Maintain the labels as a field in the class as well.
    3. Optimize the actual code even more for efficiency
    4. Accommodate the elbow method by outputting the list of annotations. This should be very simple since we already have the "pruned" list of annotations without duplicate centers.
*/
template <typename CostFunction = KMedian> 
class KCentroidsTree
{
    //Fields
   private:
    std::vector<int> index_order; //Indexes for smart access
    
    //Current solution metadata. A solution is represented by "fake" nodes that have the real nodes as children.
    Node *curr_k_solution; 
    int curr_k; 
   
   public:
    Node *tree; //The hierarchy itself

    //Methods


    /*
        Basic constructor, creates the hierarchy based on the template objective.
        Sets up smart pointers within the tree to instantly get children of a node.
    */
    KCentroidsTree(Node& root) : tree(create_hierarchy_new(root)), curr_k(1){
        int n = root.size;
        setup_quick_clusters(n, this->tree);
        this->curr_k_solution = this->tree;

    }
    /*
        Simple getter if people want to work directly with the generated tree-hierarchy instead of with this wrapper class.
    */
    Node* get_tree(){
        return this->tree;
    }

    /*
        For get_k_solution we always start from the top of the tree and do it top down instead of from a potential previous solution.
        It works by taking the nodes above the cut of any nodes with annotated k values > k.

    */
    std::vector<int> get_k_solution(int k){ //This function will always do it top down for now
        Node *curr_solution = this->tree; //Remember to delete this node when done.   
        return k_solution(k, curr_solution);
    }
    
    /*
        Same as get_k_solution except we now start from a previous solution.
        It works by taking the nodes above the cut of any nodes with annotated k values > k.
    */
    std::vector<int> inc_k_solution(){
        Node *curr_solution = this->curr_k_solution; //Remember to delete this node when done.
        int k = this->curr_k + 1;
        return k_solution(k, curr_solution);
    }

    
    private:

        /*
            Simple comparator function for sorting. This ordering with ">" ensures decreasing order, largest value first.
            "static" is required for std::sort
        */
        static bool compareByCost(const Annotation* anno1, const Annotation* anno2){ 
            return anno1->cost_decrease > anno2->cost_decrease;
        }       

        /*
            Simple cleanup function for annotations when not required anymore
        */
        void delete_annotations(std::vector<Annotation*> &annotations){
            for(Annotation* anno : annotations){
                delete anno;
            }
            return;
        }

        /*
            Main workhorse function that creates the tree itself from the annotations.
            It also annotates each node in the tree with the k value that created it. 
            This ensures it is easy to get out k solutions efficiently afterwards.
            Internal nodes end up having id of the center / cluster that its leaves will be part of.
        */
        Node* create_hierarchy_new(Node& root){
            std::vector<Annotation*> annotations = annotate_tree(root);
            std::sort(annotations.begin(), annotations.end(), compareByCost);

            Node* root_pointer = new Node{nullptr, -1.0, annotations[0]->center, 1}; //parent, cost, id, k
            annotations[0]->tree_node = root_pointer; //Unroll first iteration of loop to avoid an if/else in loop.

            Annotation* curr_anno; 
            Annotation* parent_anno;
            Node* parent_node;

            for(int i = 1; i < annotations.size(); i++){
                curr_anno = annotations[i];
                Node* new_node = new Node{nullptr, -1.0, curr_anno->center, i+1};//parent, cost, id, k
                curr_anno->tree_node = new_node;

                parent_anno = curr_anno->parent;
                parent_node = parent_anno->tree_node;
                double cost = parent_node->cost;

                if(cost != curr_anno->cost_decrease && cost >= 0){ //No more new nodes added to current parent, update the pointers
                    Node* new_parent = parent_node->children[0]; //parent, cost, id, k.
                    new_node->parent = new_parent;
                    parent_anno->tree_node = new_parent; //we are done with the previous parent node that annotation pointed to
                    new_parent->cost = curr_anno->cost_decrease;

                    Node* new_splitter = new Node{new_parent, -1.0, new_parent->id, i+1};
                    new_splitter->parent = new_parent; //TODO: Remove, right?
                    new_splitter->is_orig_cluster = true; //Used to distinguish nodes when getting clusterings out for each k //Can maybe remove this - we always know orig_cluster is the first if needed. 
                    
                    new_parent->children.push_back(new_splitter); 
                    new_parent->children.push_back(new_node); 

                } else if (cost < 0){
                    Node* new_splitter = new Node{parent_node, -1.0, parent_node->id, i+1}; //parent, cost, id, k
                    new_splitter->is_orig_cluster = true;  //Can maybe remove this - we always know orig_cluster is the first if needed. 
                    parent_node->cost = curr_anno->cost_decrease;
                    new_node->parent = parent_node;
                    
                    parent_node->children.push_back(new_splitter); 
                    parent_node->children.push_back(new_node);
                } else { //Cost_decrease == parent cost_decrease
                    parent_node->children.push_back(new_node);
                    new_node->parent = parent_node;
                }
            }
            delete_annotations(annotations); //Free the annotations from memory
            //printTree(*root_pointer);
            return root_pointer;
        }

        /*
            The function that creates the list of cost_decrease annotations. 
            It returns a list the length of number of leaves, as we only store maximal annotations for each center. 
        */
        std::vector<Annotation*> annotate_tree(Node &root){
            std::vector<Annotation*> arr;
            arr.resize(root.size);
            annotate_tree_inner(root, arr);
            return arr;
        }


        bool is_a_root(const Node &tree){
            if(tree.parent == nullptr){
                return true;
            }
            return false;
        }

        /*
            Inner tree traversal function that computes the actual annotations and appends them to the arr.
            Currently I assume the tree stores sizes - otherwise it makes this code a lot uglier - we really need to store the sizes in advance, otherwise we need to do two traversals of the list of children.
        */
        std::pair<double, int> annotate_tree_inner(Node &tree, std::vector<Annotation*> &arr){
            if((tree.children).size() == 0){
                Annotation* anno =  new Annotation{CostFunction::Evaluate((tree.parent)->cost), tree.id}; //Removed children as a list in annotations - not needed as previously thought
                
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

                std::pair<double, int> result = {best_cost, best_center};
                return result;
            }
        }

        /*
            Function that adds creates an array where each leaf id is inserted in postfix order (left to right visually).
            This means each internal node has a continuous area of that array for its children, and it gets low, high pointers to its segment of the array.
            Also assigns sizes to each internal node, i.e. how many leaves does it have.
        */
        void setup_quick_clusters(int n, Node* centroid_tree){
            std::vector<int> arr;
            arr.resize(n);
            quick_clusters_helper(centroid_tree, arr, 0);
            //printLabels(arr);
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

        /*
            Main work function for getting a specific k solution over the tree. 
            Updates the current internal solution, which is a at most 2 deep, flat, tree structure with "fake" nodes pointing to the real nodes of the k solution.
            Can take as input either the full tree or a lower k solution.
        */
        std::vector<int> k_solution(int k, Node *curr_solution){
            int n = this->tree->size;
            if(k>=n || k < 1){ //If k >= n all points should just be output
                this->curr_k = n;
                this->curr_k_solution = this->tree;
                return index_order; //Just return each id as a unique point as default for any invalid parameter for now.
            } else {
                Node *solution_holder = new Node{nullptr, 0.0};
                solution_holder->k = -1; //This has a default k value that should never be part of a solution
                get_k_solution_helper(k, curr_solution, solution_holder); //Constructs the solution

                std::vector<int> res(n);
                extract_labels(res, solution_holder); //Extract labels from solution into res
                
                this->curr_k_solution = solution_holder;
                this->curr_k = k;

                return res;
            }
        }


        /*
            The generated tree structure will never be more than two deep.
            It checks for edges crossing between lower and higher k than the search k and returns the nodes above the cut.
            Edge case is when you have multiple children, some lower and some higher than search k. Then we merge to original node those that are higher.
            Original node is the one corresponding to the center cluster that all other children took from. 
        */
        void get_k_solution_helper(int k, Node *tree, Node *new_solution){
            if(tree->size == 1){
                //std::cout << "pushing leaf " << tree->id << std::endl;
                //std::cout << tree->low << ", " << tree->high << std::endl;
                new_solution->children.push_back(tree); // We might end up with leaves part of the solution.
            }
            int min_child_k = std::numeric_limits<int>::max();
            int max_child_k = tree->k;
            for(Node *child : tree->children){
                if(max_child_k < child->k){
                    max_child_k = child->k;
                }
                if(child->k < min_child_k){
                    min_child_k = child->k;
                }
            } 

            if(max_child_k <= k){ //Recurse if no cut detected
                for(Node *child : tree->children){
                        get_k_solution_helper(k, child, new_solution);
                }
            } else{
                if(min_child_k > k){ //Simple cut
                    new_solution->children.push_back(tree);
                    //std::cout << "simple pushing c, id"  << tree->cost << ", " << tree->id << std::endl;
                    //std::cout << tree->low << ", " << tree->high << std::endl;
                } else{ //If some children have larger k and some have smaller, we need to merge things.
                    //First create the merge node:
                    Node* merge_node = new Node{nullptr, 0.0};
                    merge_node->k = -1;
                    merge_node->size = -1; //Just a "flag" value
                    //std::cout << "startfor" << std::endl;
                    for(Node *child : tree->children){
                        if(child->is_orig_cluster || child->k > k){
                            merge_node->children.push_back(child);
                            //std::cout << "merging c, id: " << child->cost << ", " << child->id <<std::endl;
                            //std::cout << tree->low << ", " << tree->high << std::endl;
                        } else{
                            new_solution->children.push_back(child);
                            //std::cout << "pushing c, id: "  << child->cost << ", " << child->id << std::endl;
                            //std::cout << tree->low << ", " << tree->high << std::endl;
                        }
                    }
                    //std::cout << "endfor" << std::endl;
                    new_solution->children.push_back(merge_node);
                }
            }
        }

        /*
            This function takes the solution container (which is a "pseudo" node) and loops through its children of real solution nodes
            to output the cluster labels. 
        */
        void extract_labels(std::vector<int> &res, Node *solution){
            int label = 0;
            for(Node *child : solution->children){
                if(child->size != -1){
                    //std::cout << "for label:" << label << ", we get " << child->low <<", " << child->high << std::endl;
                    for(int i = child->low; i <= child->high; i++){
                        int id = index_order[i];
                        res[id-1] = label;
                    }
                } else {
                    for(Node *child2 : child->children){
                        for(int i = child2->low; i <= child2->high; i++){
                            //std::cout << "for label:" << label << ", we get " << child2->low <<", " << child2->high << std::endl;
                            int id = index_order[i];
                            res[id-1] = label;
                        }
                    }
                }
                label++;
            }
        }

        /*
            Simple print function to print the output if desired
        */
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