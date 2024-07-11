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
    std::cout << "[";
    for(Annotation* anno : annotations){
        std::cout << "(" << anno->cost_decrease << ", " << anno->center << "," << anno->parent << "), "; 
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

void assign_size(Node* tree){
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