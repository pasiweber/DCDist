#include <dc_dist.hpp>
#include <quickselect.hpp>
#include <iostream>

//mlpack stuff
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
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


/*
    Also createe pybindings so that this method can be tested directly in python alongside current methods there.
*/
std::vector<double> compute_cdists(arma::mat &data, size_t k, std::string mode){

    std::cout << "compute_cdists called" << std::endl;
    //KDTree<EuclideanDistance, DTBStat, arma::mat> tree(data);
    //KDTree<EuclideanDistance, EmptyStatistic, arma::mat> tree(data);
    
    mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> *searcher;

    if(mode == "naive"){ //Brute force knn search QUICKSELECT WRITE IT YOURSELF TODO TODO
        searcher = new mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::EuclideanDistance, arma::mat, mlpack::KDTree>(data, mlpack::NAIVE_MODE);
    } else if(mode == "naive2"){
        return naive_cdists_efficient(data, k);
    } else if(mode == "parallel2"){
        return parallel_cdists2(data, k);
    } else if(mode == "parallel20"){
        return parallel_cdists20(data, k);
    } else if(mode == "parallel2000"){
        return parallel_cdists(data, k);
    } else if(mode == "parallel_arma"){
        return parallel_cdists_arma(data, k);

    }else { //KD tree knn search
        searcher = new mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::EuclideanDistance, arma::mat, mlpack::KDTree>(data); 
    }
    

    int n = data.n_cols;

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



/*
Two versions - quickselect and sorting efficiently.

1. Fill up n x n array with distances  

*/
std::vector<double> naive_cdists_efficient(arma::mat &data, size_t k){
    int n = data.n_cols;
    std::vector<std::vector<double>> dist_matrix(n, std::vector<double>(n)); //nxn dist matrix holder
    std::vector<double> cdists;
    cdists.resize(n);

    for(size_t i = 0; i < data.n_cols; i++)
    {
        const arma::vec &col_i = data.col(i);
        for(size_t j = i; j < data.n_cols; j++){
            double res = arma::norm(col_i-data.col(j), 2);
            dist_matrix[i][j] = res;
            dist_matrix[j][i] = res;
        }
        std::nth_element(dist_matrix[i].begin(), dist_matrix[i].begin() + k-1, dist_matrix[i].end());
        cdists[i] = dist_matrix[i][k-1];
    }

    return cdists;
}



std::vector<double> parallel_cdists(arma::mat &data, size_t k){
    int n = data.n_cols;
    const int dim = 2000; // Example dimension that is hardcoded for testing...
    const int minPts = k;
    
    using Point = pargeo::point<dim>; // Define the type of point
    using nodeT = pargeo::kdNode<dim, pargeo::point<dim>>;

    parlay::sequence<pargeo::point<dim>> points =  convertArmaMatToParlayPoints<dim>(data);
    
    nodeT* tree = pargeo::buildKdt<dim, pargeo::point<dim>>(points, true, true);
    parlay::sequence<size_t> nns = pargeo::kdTreeKnn<dim, Point>(points, minPts, tree, true); 

    std::vector<double> coreDist(n);

    for(int i = 0; i< points.size(); i++){
		coreDist[i] = points[nns[i*minPts + minPts-1]].dist(points[i]);
    }
    
    // parlay::parallel_for (0, points.size(), [&](int i) {
	// 		       coreDist[i] = points[nns[i*minPts + minPts-1]].dist(points[i]);
	// 		     });

    return coreDist;
}



std::vector<double> parallel_cdists2(arma::mat &data, size_t k){
    int n = data.n_cols;
    const int dim = 2; // Example dimension that is hardcoded for testing...
    const int minPts = k;
    
    using Point = pargeo::point<dim>; // Define the type of point
    using nodeT = pargeo::kdNode<dim, pargeo::point<dim>>;

    parlay::sequence<pargeo::point<dim>> points =  convertArmaMatToParlayPoints<dim>(data);
    
    nodeT* tree = pargeo::buildKdt<dim, pargeo::point<dim>>(points, true, true);
    parlay::sequence<size_t> nns = pargeo::kdTreeKnn<dim, Point>(points, minPts, tree, true); 

    std::vector<double> coreDist(n);

    for(int i = 0; i< points.size(); i++){
		coreDist[i] = points[nns[i*minPts + minPts-1]].dist(points[i]);
    }
    
    return coreDist;
}

std::vector<double> parallel_cdists20(arma::mat &data, size_t k){
    int n = data.n_cols;
    const int dim = 20; // Example dimension that is hardcoded for testing...
    const int minPts = k;
    
    using Point = pargeo::point<dim>; // Define the type of point
    using nodeT = pargeo::kdNode<dim, pargeo::point<dim>>;

    parlay::sequence<pargeo::point<dim>> points =  convertArmaMatToParlayPoints<dim>(data);
    
    nodeT* tree = pargeo::buildKdt<dim, pargeo::point<dim>>(points, true, true);
    parlay::sequence<size_t> nns = pargeo::kdTreeKnn<dim, Point>(points, minPts, tree, true); 

    std::vector<double> coreDist(n);

    for(int i = 0; i< points.size(); i++){
		coreDist[i] = points[nns[i*minPts + minPts-1]].dist(points[i]);
    }
    

    return coreDist;
}

std::vector<double> parallel_cdists_arma(arma::mat &data, size_t k){
    int n = data.n_cols;
    int dim = data.n_rows; // Example dimension that is hardcoded for testing...
    int minPts = k;
    
    using Point = pargeo::point2; // Define the type of point
    using nodeT = pargeo::kdNode2<pargeo::point2>;

    parlay::sequence<pargeo::point2> points =  convertArmaMatToParlayPoints2(data);
    
    nodeT* tree = pargeo::buildKdt2<pargeo::point2>(points, true, true);
    parlay::sequence<size_t> nns = pargeo::kdTreeKnn2<Point>(points, minPts, tree, true); 

    std::vector<double> coreDist(n);

    for(int i = 0; i< points.size(); i++){
		coreDist[i] = points[nns[i*minPts + minPts-1]].dist(points[i]);
    }
    
    // parlay::parallel_for (0, points.size(), [&](int i) {
	// 		       coreDist[i] = points[nns[i*minPts + minPts-1]].dist(points[i]);
	// 		     });

    return coreDist;
}








template<const int dim>
parlay::sequence<pargeo::point<dim>> convertArmaMatToParlayPoints(const arma::mat& mat) {

    // Create a parlay sequence of points
    parlay::sequence<pargeo::point<dim>> points(mat.n_cols);
    std::cout<< "num cols:" << mat.n_cols << std::endl;
    // Populate the sequence with points created from the matrix columns
    for (size_t j = 0; j < mat.n_cols; ++j) {
        double coords[dim];
        for (size_t i = 0; i < dim; ++i) {
            coords[i] = mat(i, j);
        }
        points[j] = pargeo::point<dim>(coords);
    }

    return points;
}


parlay::sequence<pargeo::point2> convertArmaMatToParlayPoints2(const arma::mat& mat) {
    // Create a parlay sequence of points
    parlay::sequence<pargeo::point2> points(mat.n_cols);
    std::cout<< "num cols:" << mat.n_cols << std::endl;
    // Populate the sequence with points created from the matrix columns
    for (size_t j = 0; j < mat.n_cols; ++j) {
        arma::Col<double> coords(mat.n_rows);
        for (size_t i = 0; i < mat.n_rows; ++i) {
            coords[i] = mat(i, j);
        }
        points[j] = pargeo::point2(&coords);
    }

    return points;
}