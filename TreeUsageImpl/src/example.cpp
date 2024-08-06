#include <iostream>
#include <vector>
#include <cmath>
#include <kcentroids.hpp>
#include <kcentroids2.hpp>
#include <dc_hdbscan.hpp>
#include <dc_dist.hpp>
#include <key_structs.hpp>

#include <armadillo>

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include <../parallel_hdbscan/src/kdTree.h>
#include <../parallel_hdbscan/src/kdTreeKnn.h>
#include <../parallel_hdbscan/include/hdbscan/point.h>
#include <../parallel_hdbscan/src/kdTreeArma.h>
#include <../parallel_hdbscan/src/kdTreeKnnArma.h>
#include <../parallel_hdbscan/include/hdbscan/armapoint.h>
#include <../parallel_hdbscan/include/hdbscan/hdbscan.h>


// Function to create a new node -
Node* addNode(Node* parent = nullptr, double cost=0.0, int id = -1, int size=1) {
    Node* newNode = new Node;
    newNode->cost = cost;
    newNode->parent = parent;
    newNode->id = id;
    newNode->size = size;

    if(parent != nullptr){
        parent->children.push_back(newNode);
    }

    return newNode;
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


parlay::sequence<pargeo::ArmaPoint> convertArmaMatToParlayPoints22(const arma::mat& mat) {
    // Create a parlay sequence of points
    parlay::sequence<pargeo::ArmaPoint> points(mat.n_cols);
    std::cout<< "num cols:" << mat.n_cols << std::endl;
    // Populate the sequence with points created from the matrix columns
    for (size_t j = 0; j < mat.n_cols; ++j) {
        arma::Col<double> coords(mat.n_rows);
        for (size_t i = 0; i < mat.n_rows; ++i) {
            coords[i] = mat(i, j);
        }
        points[j] = pargeo::ArmaPoint(&coords);
    }

    return points;
}

parlay::sequence<pargeo::VectorPoint> convertArmaMatToParlayPointsVector(const arma::mat& mat) {
    // Create a parlay sequence of points
    parlay::sequence<pargeo::VectorPoint> points(mat.n_cols);
    std::cout<< "num cols:" << mat.n_cols << std::endl;
    // Populate the sequence with points created from the matrix columns
    for (size_t j = 0; j < mat.n_cols; ++j) {
        std::vector<double> coords(mat.n_rows);
        for (size_t i = 0; i < mat.n_rows; ++i) {
            coords[i] = mat(i, j);
        }
        points[j] = pargeo::VectorPoint(&coords);
    }

    return points;
}

// Function to generate a tree with 11 leaf nodes
Node* generateTree11() {
    Node* root = addNode(nullptr, 100.0, -1, 11);  // Create the root node with an arbitrary cost

    // Create internal nodes and leaves
    Node* node1 = addNode(root, 90.0, -1, 4);
    Node* node2 = addNode(root, 80.0, -1, 4);
    Node* node3 = addNode(root, 28.0, -1, 3);

    Node* l1 = addNode(node3, 0.0, 1);
    Node* l2 = addNode(node3, 0.0, 2);
    Node* l3 = addNode(node3, 0.0, 3);


    Node* node4 = addNode(node1, 55.0, -1, 2);
    Node* node5 = addNode(node1, 35.0, -1, 2);

    Node* node6 = addNode(node2, 50.0, -1, 3);
    Node* node7 = addNode(node6, 40.0, -1, 2);

    Node* l4 = addNode(node4, 0.0, 4);
    Node* l5 = addNode(node4, 0.0, 5);
    Node* l6 = addNode(node5, 0.0, 6);
    Node* l7 = addNode(node5, 0.0, 7);

    Node* l8 = addNode(node2, 0.0, 8);
    Node* l9 = addNode(node6, 0.0, 9);
    Node* l10 = addNode(node7, 0.0, 10);
    Node* l11 = addNode(node7, 0.0, 11);

    return root;
}

Node* generateTree12() {
    Node* root = addNode(nullptr, 13.9, -1, 12);  // Create the root node with an arbitrary cost

    Node* l5 = addNode(root, 0.0, 5);
    Node* node94 = addNode(root, 9.4, -1, 11);

    Node* node22 = addNode(node94, 2.2, -1, 4);
    Node* node51 = addNode(node94, 5.1, -1, 7);

    Node* l4 = addNode(node22, 0.0, 4);
    Node* node20 = addNode(node22, 2.0, -1, 3);

    Node* l2 = addNode(node20, 0.0, 2);
    Node* node14 = addNode(node20, 1.4, -1, 2);

    Node* l1 = addNode(node14, 0.0, 1);
    Node* l3 = addNode(node14, 0.0, 3);

    Node* l8 = addNode(node51, 0.0, 8);
    Node* node45 = addNode(node51, 4.5, -1, 6);

    Node* node30 = addNode(node45, 3.0, -1, 3);
    Node* node32 = addNode(node45, 3.2, -1, 3);

    Node* l6 = addNode(node30, 0.0, 6);
    Node* l7 = addNode(node30, 0.0, 7);
    Node* l9 = addNode(node30, 0.0, 9);

    Node* l10 = addNode(node32, 0.0, 10);
    Node* l11 = addNode(node32, 0.0, 11);
    Node* l12 = addNode(node32, 0.0, 12);
    return root;
}


void printLabels(std::vector<int> labels){
    std::cout << "Labels:" << std::endl;
    std::cout << "[";
    for(int label : labels){
        std::cout << label << ", ";
    }

    std::cout << "]" << std::endl;

}


double kmedian(double x){
    return x;
}
double kmeans(double x){
    return std::pow(x,2);
}



void test_k_centroids(){
    Node* root = generateTree12();
    std::cout << "Generated dc-tree Tree:\n";
    std::vector<Annotation*> res = annotate_tree(*root, kmeans);
    //printTree(*root);



    Node* rootv2 = create_hierarchy(*root, kmeans);
    std::cout << "kmeans tree:" << std::endl;
    printTree(*rootv2);

    // Create an instance of KCentroidsTree
    KCentroidsTree<KMeans> tree(*root);

    // Use the tree
    Node* rootv3 = tree.get_tree();

    //std::cout << "kmeans tree new:" << std::endl;


    std::vector<int> res_new = tree.get_k_solution(6);
    
    //std::vector<int> labels = kcentroids(*root, kmeans, 6);
    std::cout << "next k" << std::endl;
    //std::vector<int> labels2 = kcentroids(*root, kmeans, 7);
    std::vector<int> res_new2 = tree.inc_k_solution();
    printTree(*rootv3);
    std::sort(res.begin(), res.end(), compareByCost);

    //print_annotations(res);

}

//How do I handle providing a tree in an elegant way?
void test_hdbscan(){
    Node* root = generateTree12();
    assign_sizes(root);
    printTree(*root);
    int mpts = 3; //TODO: Check whether the algorithm is robust to this parameter when we are able to construct the dc-tree
    int mcs = 2; //For now we only allow mcs>=2 since we do not have any cdists.

    Dc_hdbscan tree(mpts, mcs);
    tree.fit(root);
    std::vector<int> labels = tree.labels_;
    printLabels(labels);

}



void test_parallel_hdbscan(){
    const int dim = 3; // Example dimension
    const int minPts = 2;
    using Point = pargeo::point<dim>; // Define the type of point
    using nodeT = pargeo::kdNode<dim, pargeo::point<dim>>;

    // Initialize a parlay sequence of points
    parlay::sequence<Point> points(5); // Create a sequence of 5 points

    // Populate the sequence with some values
    for (int i = 0; i < 5; ++i) {
        double coords[dim] = {static_cast<double>(i), static_cast<double>(i+1), static_cast<double>(i+2)};
        points[i] = Point(coords);
    }

    nodeT* tree = pargeo::buildKdt<dim, pargeo::point<dim>>(points, true, true);
    parlay::sequence<size_t> nns = pargeo::kdTreeKnn<dim, Point>(points, minPts, tree, true); 

    parlay::sequence<float> coreDist = parlay::sequence<float>(points.size());
    parlay::parallel_for (0, points.size(), [&](int i) {
			       coreDist[i] = points[nns[i*minPts + minPts-1]].dist(points[i]);
			     });


    for (const float& value : coreDist) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}


void test_mlpack(){
    arma::mat data3 = { {0.0, 1.0, 20.0, 0.0, 1.0, 20.0, 0.0, 1.0, 20.0, 0.0, 1.0, 20.0, 0.0, 1.0, 20.0, 0.0, 1.0, 20.0, 0.0, 1.0, 20.0, 0.0, 1.0, 20.0},
                        {0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0},
                        {0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0},
                        {0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0},
                        {0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0},
                        {0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0},
                        {0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0},
                        {0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0},
                        {0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0},
                        {0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0},
                        {0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0},
                        {0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0},
                        
                        
                        {0.1, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0, 0.0, 2.0, 35.0},
                        {0.2, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0, 0.0, 3.0, 30.0},
                        {0.3, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0, 0.0, 1.5, 30.0},
                        {0.4, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0, 0.0, 5.5, 11.0},
                        {0.5, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0, 0.0, 10.0,28.0},
                        {0.6, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0, 0.0, 2.0, 30.0},
                        {0.7, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0, 0.0, 7.5, 11.0},
                        {0.8, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0, 0.0, 12.0,28.0},
                        {0.9, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0, 0.0, 4.5, 30.0},
                        {0.10, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0, 0.0, 5.5, 21.0},
                        {0.11, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0, 0.5, 11.0,28.0},
                        };

    arma::mat data2 = arma::trans(data3);
    const int minPts = 2;
    int n = data2.n_cols;
    parlay::sequence<pargeo::ArmaPoint> points2 =  convertArmaMatToParlayPoints22(data2);
    parlay::sequence<pargeo::point<24>> points =  convertArmaMatToParlayPoints<24>(data2);
    pargeo::ArmaPoint p1;
    pargeo::ArmaPoint p2(10); 

    parlay::sequence<pargeo::VectorPoint> pointsVector =  convertArmaMatToParlayPointsVector(data2);
    
    //std::cout << "size of points2: " << sizeof(points2[0]) << ", and size of points: " << sizeof(points[0]) << std::endl;
    //std::cout << "size of p1: " << sizeof(p1) << ", and size of p2: " << sizeof(p2) << std::endl;

    //VECTOR
    parlay::sequence<pargeo::wghEdge> E3;
    E3 = pargeo::hdbscan_vector(pointsVector, minPts);
    parlay::sequence<pargeo::dendroNode> dendro3 = pargeo::dendrogram(E3, n);
    parlay::sequence<double> A3(dendro3.size()*4);
    parlay::parallel_for(0, dendro3.size(), [&](size_t i){
                        A3[i*4+0] = std::get<0>(dendro3[i]);
                        A3[i*4+1] = std::get<1>(dendro3[i]);
                        A3[i*4+2] = std::get<2>(dendro3[i]);
                        A3[i*4+3] = std::get<3>(dendro3[i]);});

    std::cout << "1 ###########################################################################################" << std::endl;

    //PARGEO
    parlay::sequence<pargeo::wghEdge> E2;
    E2 = pargeo::hdbscan<24>(points, minPts);
    parlay::sequence<pargeo::dendroNode> dendro2 = pargeo::dendrogram(E2, n);
    parlay::sequence<double> A2(dendro2.size()*4);
    parlay::parallel_for(0, dendro2.size(), [&](size_t i){
                        A2[i*4+0] = std::get<0>(dendro2[i]);
                        A2[i*4+1] = std::get<1>(dendro2[i]);
                        A2[i*4+2] = std::get<2>(dendro2[i]);
                        A2[i*4+3] = std::get<3>(dendro2[i]);});

    std::cout << "2 ###########################################################################################" << std::endl;
    //ARMA
    parlay::sequence<pargeo::wghEdge> E;
    E = pargeo::hdbscan_arma(points2, minPts);
    parlay::sequence<pargeo::dendroNode> dendro = pargeo::dendrogram(E, n);
    parlay::sequence<double> A(dendro.size()*4);
    parlay::parallel_for(0, dendro.size(), [&](size_t i){
                        A[i*4+0] = std::get<0>(dendro[i]);
                        A[i*4+1] = std::get<1>(dendro[i]);
                        A[i*4+2] = std::get<2>(dendro[i]);
                        A[i*4+3] = std::get<3>(dendro[i]);});


    std::cout << "n: " << n<< ", dims: " << data2.n_rows <<std::endl;
    std::cout << "n: " << points2.size() << ", dims: " << points2[0].dim <<std::endl;

}


void simpleTest(){

    arma::Col<double> a = {4.0,5.0,6.0};

    arma::Col<double> c(a);

    arma::Col<double> d = {1.0,2.0,3.0};

    pargeo::ArmaPoint point(&d);

    arma::Col<double> *coords = point.coords();
    std::cout << point.dim << " dimensions" << std::endl;
    pargeo::ArmaPoint ArmaPoint(point.coords());
    std::cout << ArmaPoint.dim << " dimensions" << std::endl;


    
}



int main() {
    //test_k_centroids();
    //test_hdbscan();
    //test_parallel_hdbscan();
    //simpleTest();
    test_mlpack();
    return 0;
}
