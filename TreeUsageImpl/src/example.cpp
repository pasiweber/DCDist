#include <iostream>
#include <vector>
#include <cmath>
#include <kcentroids.hpp>
#include <kcentroids2.hpp>
#include <dc_hdbscan.hpp>
#include <dc_dist.hpp>
#include<mlpack.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <mlpack/core.hpp>
#include <mlpack/core/tree/binary_space_tree.hpp>



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
    print_annotations(res);
    printTree(*root);



    Node* rootv2 = create_hierarchy(*root, kmeans);
    std::cout << "kmeans tree:" << std::endl;
    printTree(*rootv2);

    // Create an instance of KCentroidsTree
    KCentroidsTree<KMeans> tree(*root);

    // Use the tree
    Node* rootv3 = tree.get_tree();

    //std::cout << "kmeans tree new:" << std::endl;
    //printTree(*rootv3);


    std::vector<int> res_new = tree.get_k_solution(4);
    
    std::vector<int> labels = kcentroids(*root, kmeans, 4);
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


void test_mlpack(){
    using namespace mlpack::metric; // ManhattanDistance
    using namespace mlpack::tree;
    using namespace mlpack::neighbor;

    // arma::mat data = {  {1.0, 1.0},
    //                     {2.0, 2.0},
    //                     {3.0, 3.0},
    //                     {1.5, 1.5},
    //                     {5.5, 5.5},
    //                     {10.0, 10.0}
    //                     };
     
    arma::mat data2 = {{0.0, 1.0},
                        {0.0, 2.0},
                        {0.0, 3.0},
                        {0.0, 1.5},
                        {0.0, 5.5},
                        {0.0, 10.0}
                        };


    //arma::vec query = data.row(0);
    arma::mat data3 = arma::trans(data2);
    compute_cdists(data3, 2, "kdtree");
    


}

int main() {
    test_k_centroids();
    //test_hdbscan();
    //test_mlpack();
    return 0;
}
