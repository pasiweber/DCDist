
#include <vector>
#include <kcentroids.hpp>
#include <hcf.hpp>

HCF::HCF(int min_pts, int min_cluster_size, double (*obj_func)(Node*)) : 
tree(nullptr){  //Initialize tree here
    if(min_pts < 1 || min_cluster_size < 1){
        throw std::invalid_argument("min_pts and min_cluster_size should be 1 or higher");
    }
    this->min_pts = min_pts; 
    this->min_cluster_size = min_cluster_size;
    this->objective_function = obj_func;
}

void HCF::fit(Node *tree){
    std::cout << "Tree already provided..." << std::endl;
    printTree(*tree);

    compute_clustering(tree);
}


//The set of clusters will be a vector of tree references I think... must be faster this way
//I know how many internal nodes we have, I guess I can just maintain a bool array...
void HCF::compute_clustering(Node *tree){
    std::cout << "computing clustering on tree with size: " << tree->size << std::endl;
    bottom_up_cluster(tree);
    label_clusters(tree);

}
/*
    Skeleton structure of the general HCF structure. 
*/
double HCF::bottom_up_cluster(Node *tree){
    if((this->min_cluster_size) > (tree->size)){ //noise node
        tree->is_cluster = false;
        return 0.0;

    } else{ //internal node
        double total_cluster_objective = 0.0; //This will contain the sum of stabilities of best clusters from below
        int num_children = tree->children.size();

        for(Node *child : tree->children){
            total_cluster_objective += bottom_up_cluster(child);
        }
        double new_stability = this->objective_function(tree);

        if(new_stability >= total_cluster_objective){
            tree->is_cluster = true;
            return new_stability;
        }
        tree->is_cluster = false;
        return total_cluster_objective;
    }
}


