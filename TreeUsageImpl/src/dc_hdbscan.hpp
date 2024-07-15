#ifndef DC_HDBSCAN_HPP
#define DC_HDBSCAN_HPP

#include <vector>
#include <kcentroids.hpp>

class Dc_hdbscan
{

    //Fields
   private:
    int min_pts;
    int min_cluster_size;
    std::vector<int> labels_;
    std::vector<double> cdists;
    Node *tree; //Should be a pointer rather than a reference due to c++ nature.

    
   public:
    
    //Constructors
    Dc_hdbscan(int min_pts, int min_cluster_size);


    //Methods


    /*
        Finds the HDBSCAN clustering of dataset in parameter points.
        1) Construct the dc-tree if not already provided.
        2) Compute clusters bottom up
        3) Assign internal fields for the labels
    */
    void fit(const std::vector<std::vector<double>> &points);

    void fit(Node *tree);

   private:


    void compute_clustering(Node *tree);

    void bottom_up_cluster(Node *tree); //This signature should be changed depending on how the actual assigning and maintenance of current clusters will be.
};




#endif