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
    this->min_pts = min_pts; 
    this->min_cluster_size = min_cluster_size;

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
    std::cout << "computing clustering on tree with size: " << tree->size << std::endl;
    bottom_up_cluster(tree, true);

}
/*
    If I maintain metadata about split below I can avoid recursing through tree multiple times (although only constant factor)
    I just return the running sum of when points fall out within the current cluster region, then I only need to add the sum of the max value for which all the points are in the cluster.

    I should still do the final computation in a function -> makes it more clear what is happening

    Return (total stability, stability from below node(s))
*/
std::pair<double, double> Dc_hdbscan::bottom_up_cluster(Node *tree, bool merge_above){
    if((this->min_cluster_size) > (tree->size)){ //noise node
        tree->is_cluster = false;
        double stability_contribution = (tree->size)*(1.0/(tree->parent->cost)); //TODO check how this casts
        return {0,stability_contribution};

    } else{ //internal node
        double total_cluster_stability = 0.0; //This will contain the sum of staabilities of best clusters from below
        double total_region_contribution = 0.0; //This will contain the sum of the level points fall out of current cluster region
        int s_size = split_size(tree, this->min_cluster_size); 

        for(Node *child : tree->children){
            std::pair<double, double> res;
            if(s_size >= 2){ //If we have a split at this level, this means that we should compute the stability at below recursive level.
                res = bottom_up_cluster(child, true);
            } else{
                res = bottom_up_cluster(child, false);
            }
            total_cluster_stability += res.first;
            total_region_contribution += res.second;

        }

        if(tree->parent == nullptr){ //root node //Here we can implement allowing a singular cluster to be output if want
           return {0.0,0.0}; //Just return default values in parent - nothing left to do
        }

        if(merge_above){ //merge_above is technically not needed, however we avoid multiple calls to stability and are able to pass cluster region metadata in the way we currently do.
            double pcost = tree->parent->cost;
            double new_stability = stability(tree->size, pcost, total_region_contribution);
            std::cout << "stability at node " << tree->cost << " is: "<<new_stability << std::endl;

            total_region_contribution = tree->size / pcost; //Since there is a split above all these points fall out of above cluster at pdist //TODO: check conversion int->double
            if(new_stability >= total_cluster_stability){
                tree->is_cluster = true;
                return {new_stability, total_region_contribution};
            }
            return {total_cluster_stability, total_region_contribution};
        } else{

            tree->is_cluster = false;
            return {total_cluster_stability, total_region_contribution};
        }
    }
}


/*
    This method ensures that we do not count noise branches in the split size
*/
int Dc_hdbscan::split_size(Node *tree, int mcs){
    int size = 0;
    if(mcs != 1){
        for(Node *child : tree->children){
            if(child->size >= mcs){
                size++;
            }
        }

    } else{ //We might just remove this else branch if we decide to not include mcs=1
        for(Node *child : tree->children){
            if(child->size == 1 && this->cdists[tree->id] == tree->cost){
                continue;
            }
            size++;
        }
    }
    return size;
}


double Dc_hdbscan::stability(int size, double pdist, double fallout_sum){
    return fallout_sum - size/pdist; //TODO: check conversion
}
