#include <dc_dist.hpp>
//#include <quickselect.hpp>
#include <iostream>

//mlpack stuff
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <mlpack/core/tree/binary_space_tree.hpp>

#include <graph_hdb_c.hpp>
#include <graph.hpp>


void assign_node_sizes(Node* tree){
    assign_node_size_helper(tree);
}

int assign_node_size_helper(Node* tree){
    if((tree->children).size() == 0){
        tree->size = 1;
        return 1;
    }{ 
        int size = 0;
        for(Node* child : tree->children){
            size += assign_node_size_helper(child);
        }
        tree->size = size;
        return size;
    }
}


//TODO: Maybe move stuff around? Right now my files are in weird places
Node* construct_dc_tree(double *data, unsigned long long n, int dim, int k){ //k is minPts for core dists
    double* mut_dists = calc_mutual_reachability_dist(data, n, dim, k);
    std::vector<MSTEdge> edges = calc_mst(n, mut_dists, k);
    Node* root = constructHierarchy(edges);
    assign_node_sizes(root);

    return root;
}




std::vector<double> compute_cdists(arma::mat &data, size_t k, std::string mode){
    mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> *searcher;

    searcher = new mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::EuclideanDistance, arma::mat, mlpack::KDTree>(data); 

    int n = data.n_cols;

    arma::Mat<size_t> neighbors;
    arma::mat distances;

    searcher->Search(k-1, neighbors, distances); //This does not include point itself in k, which is why we do k-1

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


// Function to print the tree (for debugging purposes)
void printSubtree2(const std::string &prefix, const Node& tree) {
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
                std::cout << "(L" << c->k << ")\n";
            } else{
                std::cout << "(" << c->k << ")\n";
            }
            printSubtree2(newPrefix, *c);
        } else {
            cout << (n_children > 1 ? prefix : "") << "└────";
            if(c->children.empty()){
                std::cout << "(L" << c->k << ")\n";
            } else{
                std::cout << "(" << c->k << ")\n";
            }
            printSubtree2(prefix + "\t", *c);
        }
    }
}

void printTree2(const Node& tree) {
    using std::cout;
    cout << tree.k << "\n";
    printSubtree2("", tree);
    cout << "\n";
}










