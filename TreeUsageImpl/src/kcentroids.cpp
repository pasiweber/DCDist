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
        std::cout << "(" << anno->cost_decrease << ", " << anno->center << "), "; 
    }
    std::cout << "]" << std::endl;

}