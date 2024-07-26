#include <kcentroids.hpp>

#include <vector>
#include <limits>


bool is_a_root(const Node &tree){
    if(tree.parent == nullptr){
        return true;
    }
    return false;
}


bool compareByCost(const Annotation* anno1, const Annotation* anno2){
    return anno1->cost_decrease > anno2->cost_decrease;
}

void print_annotations(std::vector<Annotation*> annotations){
    std::cout << "Annotations:" << std::endl;
    std::cout << "[";
    for(Annotation* anno : annotations){
        std::cout << "(" << anno->cost_decrease << ", " << anno->center << "), "; 
    }
    std::cout << "]" << std::endl;

}




int split_detector(std::vector<Annotation*> annotations, int i){
    std::vector<Annotation*> arr;
    int parent = annotations[i]->center; 
    i++;
    int ctr = 0;
    int cd = 0;
    double first_cost_decrease = 0.0;
    for(i; i< annotations.size(); i++){
        cd = annotations[i]->cost_decrease;
        if(ctr == 0){
            first_cost_decrease = cd;
        } 
        if(cd == first_cost_decrease){
           ctr++;
        } else {
            break;
        }
    }
    return ctr;

}



void delete_tree(Node* tree){
    for (Node* child : tree->children) {
        delete_tree(child);
    }

    delete tree;

}

void assign_sizes(Node* tree){
    assign_size_helper(tree);
}

int assign_size_helper(Node* tree){
    if((tree->children).size() == 0){
        tree->size = 1;
        return 1;
    }{ 
        int size = 0;
        for(Node* child : tree->children){
            size += child->size;
        }
        tree->size = size;
        return size;
    }

}



void delete_annotations(std::vector<Annotation*> &annotations){
    for(Annotation* anno : annotations){
        delete anno;
    }
    return;
}


//This puts the center annotation in id of internal nodes and cost of leaf nodes
void annotate_tree_centers(std::vector<Annotation*> &annotations, int k){
    for(int i = 0; i < k; i++){
        Node* node = annotations[i]->orig_node;
        if(node->children.size() !=0){
            node->id = annotations[i]->center;
        } else{
            node->cost =annotations[i]->center;
        }
        //std::cout << "Annotating node with dist" << node->cost << " with center " << node->id << std::endl;
    }
    return;
}


std::vector<int> assign_points(Node &tree){
    int n = tree.size;
    std::vector<int> arr;
    arr.resize(n);
    assign_points_helper(tree, arr, -1);
    return arr;
}

void assign_points_helper(Node &tree, std::vector<int> &arr, int curr_center){
    if(tree.children.size()==0){
        //std::cout << "curr center:" << curr_center << std::endl;
        if(tree.cost!= 0){ //For leaf nodes the annotated center is stored in the cost
            arr[tree.id-1] = tree.cost;
        } else{
            arr[tree.id-1] = curr_center;
        }
    } else {
        int center = curr_center;
        if(tree.id != -1){
            center = tree.id;
        }
        for(Node* child : tree.children){
            assign_points_helper(*child, arr, center);
        }

    }
}


//Reset all internal node ids to -1 and all leaf node costs to 0
void cleanTree(Node &tree){
    if(tree.children.size()==0){
        tree.cost = 0.0;
    } else {
        tree.id = -1;
        for(Node* child : tree.children){
            cleanTree(*child);
        }
    }
}



