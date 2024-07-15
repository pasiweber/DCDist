#include <dc_hdbscan.hpp>
#include <kcentroids.hpp>

#include <vector>
#include <limits>
#include <stdexcept>



Dc_hdbscan::Dc_hdbscan(int min_pts, int min_cluster_size) : 
tree(nullptr){  //Initialize tree here
    if(min_pts < 1 || min_cluster_size < 1){
        throw std::invalid_argument("min_pts and min_cluster_size should be 1 or higher");
    }
    min_pts = min_pts; 
    min_cluster_size = min_cluster_size;

}



void Dc_hdbscan::fit(const std::vector<std::vector<double>> &points){
    std::cout << "Constructing dc-tree from points..." << std::endl;
    tree = construct_dc_tree(points);
    printTree(*tree);

    compute_clustering(tree);
}

void Dc_hdbscan::fit(Node *tree){
    std::cout << "Tree already provided..." << std::endl;
    printTree(*tree);

    compute_clustering(tree);
}

//The set of clusters will be a vector of tree references I think... must be faster this way
//I know how many internal nodes we have, I guess I can just maintain a bool array...
void Dc_hdbscan::compute_clustering(Node *tree){
    std::cout << "computing clustering" << std::endl;


}