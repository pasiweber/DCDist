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