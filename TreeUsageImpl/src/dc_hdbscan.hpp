#ifndef DC_HDBSCAN_HPP
#define DC_HDBSCAN_HPP

#include <vector>
#include <kcentroids.hpp>

/*
    Efficient focused implementation of HDBSCAN over the dc-tree.
    Handles all parameter values (when tested) of min_cluster_size and min_pts. 

*/
class Dc_hdbscan
{

    //Fields
   private:
    int min_pts;
    int min_cluster_size;
    std::vector<double> cdists;
    Node *tree; //Should be a pointer rather than a reference due to c++ nature.

    
   public:
    std::vector<int> labels_;


    //Constructors
    Dc_hdbscan(int min_pts, int min_cluster_size);


    //Methods

    /*
        Finds the HDBSCAN clustering of dataset in parameter points.
        1) Construct the dc-tree if not already provided.
        2) Compute clusters bottom up
        3) Assign internal fields for the labels
    */
    void fit(double *data, unsigned long long n, int dim, int minPts);

    void fit(Node *tree);

   private:


    void compute_clustering(Node *tree);
    //TODO: Create a version without the merge above - this is only here to speed up HDBSCAN's stability specifically. 
    std::pair<double, double> bottom_up_cluster(Node *tree, bool merge_above); //This signature should be changed depending on how the actual assigning and maintenance of current clusters will be.

    int split_size(Node *tree, int mcs);

    double stability(int size, double pdist, double fallout_sum);

    void label_clusters(Node *tree);

    int label_clusters_helper(Node *tree, std::vector<int> &labels, bool within_cluster, int ctr);
};




#endif