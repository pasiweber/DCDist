#include <dc_dist.hpp>
#include <iostream>
#include <string>

//mlpack stuff
#include <mlpack.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <mlpack/core.hpp>
#include <mlpack/core/tree/binary_space_tree.hpp>


//TODO: Implement
Node* construct_dc_tree(const std::vector<std::vector<double>> &points){
    return new Node{};
}




// Function to print the tree (for debugging purposes)
void printSubtree(const std::string &prefix, const Node& tree) {
    using std::cout;
    using std::endl;
    if ((tree.children).size() == 0) return;
    cout << prefix;
    size_t n_children = tree.children.size();
    cout << (n_children > 1 ? "├────" : "");

    for (size_t i = 0; i < n_children; ++i) {
        Node *c = tree.children[i];
        if (i < n_children - 1) {
            if (i > 0) { // added fix
                cout << prefix<< "├────"; // added fix
            } // added fix
            bool printStrand = n_children > 1 && !c->children.empty();
            std::string newPrefix = prefix + (printStrand ? "│\t" : "\t");
            if(c->children.empty()){
                std::cout << "(L" << c->id << ")\n";
            } else{
                std::cout << "(" << c->cost << ")\n";
            }
            printSubtree(newPrefix, *c);
        } else {
            cout << (n_children > 1 ? prefix : "") << "└────";
            if(c->children.empty()){
                std::cout << "(L" << c->id << ")\n";
            } else{
                std::cout << "(" << c->cost << ")\n";
            }
            printSubtree(prefix + "\t", *c);
        }
    }
}

void printTree(const Node& tree) {
    using std::cout;
    cout << tree.cost << "\n";
    printSubtree("", tree);
    cout << "\n";
}

/*
    Also createe pybindings so that this method can be tested directly in python alongside current methods there.
*/
std::vector<double> compute_cdists(arma::mat &data, size_t k, std::string mode){
    using namespace mlpack::metric;
    using namespace mlpack::tree;
    using namespace mlpack::neighbor;
    std::cout << "compute_cdists called" << std::endl;
    //KDTree<EuclideanDistance, DTBStat, arma::mat> tree(data);
    //KDTree<EuclideanDistance, EmptyStatistic, arma::mat> tree(data);
    
    NeighborSearch<NearestNeighborSort, EuclideanDistance, arma::mat, KDTree> *searcher;

    if(mode == "naive"){ //Brute force knn search
        searcher = new NeighborSearch<NearestNeighborSort, EuclideanDistance, arma::mat, KDTree>(data, NAIVE_MODE);
    } else { //KD tree knn search
        searcher = new NeighborSearch<NearestNeighborSort, EuclideanDistance, arma::mat, KDTree>(data); 
    }
    
    arma::Mat<size_t> neighbors;
    arma::mat distances;

    searcher->Search(k-1, neighbors, distances); //This does not include point itself in k, which is why we do k-1
    //std::cout << "The neighbors:" << std::endl; 
    //neighbors.print();

    //std::cout << "The distances:" << std::endl; 
    //distances.print();

    std::vector<double> cdists = extract_cdists(distances, k);

    return cdists; //This default return value needs to be here, otherwise armadillo becomes confused...
}



std::vector<double> extract_cdists(arma::mat distances, size_t k){
    std::vector<double> cdists;

    std::cout << "num_cols:" << distances.n_cols << std::endl;
    for(size_t i = 0; i < distances.n_cols; i++)
    {
        arma::vec col = distances.col(i);
        cdists.push_back(col(k-2));
    }
    
    return cdists;
}
