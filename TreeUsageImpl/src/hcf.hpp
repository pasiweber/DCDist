#ifndef HCF_HPP
#define HCF_HPP

#include <vector>
#include <key_structs.hpp>

/*
    Generalized implementation of the Hierarchical Clustering Framework.
    Can be parameterized with:
        - Tree structure (TODO)
        - Objective function

    Questions in terms of implementation:
    1. Should we accumulate any metadata within the recursive calls giving the option of further optimizing the objective function calls?
    2. How should we parameterize w.r.t. tree construction?
    3. Generally how permissive should we be when it comes to the objective function?
    4. How should we handle mcs? In general, how much of the structure of the recursive structure should we enforce and how much should be templated?
    5. Currently the function passed cannot use the internal state of the HCF - should we allow this, if so, how?
*/


class HCF
{
    //Fields
   private:
    int min_pts;
    int min_cluster_size;
    std::vector<double> cdists;
    Node *tree; //Should be a pointer rather than a reference due to c++ nature.
    double (*objective_function)(Node*); //The parameterized objective function. 
    
   public:
    std::vector<int> labels_;


    //Constructors
    HCF(int min_pts, int min_cluster_size, double (*obj_func)(Node*));


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
    double bottom_up_cluster(Node *tree); //This signature should be changed depending on how the actual assigning and maintenance of current clusters will be.

    void label_clusters(Node *tree);

    int label_clusters_helper(Node *tree, std::vector<int> &labels, bool within_cluster, int ctr);
};



#endif //HCF