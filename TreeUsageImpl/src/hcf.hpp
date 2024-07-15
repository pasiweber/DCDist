#ifndef HCF_HPP
#define HCF_HPP

#include <vector>
#include <kcentroids.hpp>

/*
    Generalized implementation of the Hierarchical Clustering Framework.
    Can be parameterized with:
        - Tree structure
        - Objective function
*/

//TODO: Template the function based on the type of tree construction it should do. Either the dc-tree class or the kcentroids tree class (I guess).
//I will start by just implementing the fit that takes the tree as input. Then we can always parameterize it afterwards however..
//The compute clustering does not care about how the tree is constructed either way.
class HCF
{
    //Fields
   private:
    int min_pts;
    int min_cluster_size;
    std::vector<double> cdists;
    Node *tree; //Should be a pointer rather than a reference due to c++ nature.
    int (*objective_function)(Node*, double, double); //The parameterized objective function. 
    
   public:
    std::vector<int> labels_;


    //Constructors
    HCF(int min_pts, int min_cluster_size, int (*obj_func)(Node*, double, double));


    //Methods

    /*
        Finds the HDBSCAN clustering of dataset in parameter points.
        1) Construct the tree if not already provided.
        2) Compute clusters bottom up
        3) Assign internal fields for the labels
    */
    void fit(const std::vector<std::vector<double>> &points);

    void fit(Node *tree);

   private:


    void compute_clustering(Node *tree);
    //TODO: Create a version without the merge above - this is only here to speed up HDBSCAN's stability specifically. 
    std::pair<double, double> bottom_up_cluster(Node *tree, bool merge_above); //This signature should be changed depending on how the actual assigning and maintenance of current clusters will be.

    void label_clusters(Node *tree);

    int label_clusters_helper(Node *tree, std::vector<int> &labels, bool within_cluster, int ctr);
};



#endif //HCF