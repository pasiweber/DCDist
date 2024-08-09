#include <iostream>
#include <vector>
#include <cmath>
#include <kcentroids_efficient.hpp>
#include <dc_hdbscan.hpp>
#include <dc_dist.hpp>
#include <key_structs.hpp>

//#include <armadillo>

#include <graph_hdb_c.hpp>
#include <graph.hpp>

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



void test_mst(){

    double data[] = {0.5, 0.6,
                     2.0, 2.0,
                     3.0, 3.0,
                     4.0, 4.0,
                     50.0, 50.0,
                     21.0, 10.0,
                     36.0, 8.0,};

    unsigned long long n = 7;
    int dim = 2;
    int k = 2;

    auto mut_dists = calc_mutual_reachability_dist(data, n, dim, k);
    std::vector<MSTEdge> edges = calc_mst(n, mut_dists, k);

    std::cout << edges.size() << " MST edges (src, dst, weight): ";
    for(MSTEdge edge : edges){
        std::cout << "(" << edge.src << ", " << edge.dst << ", " << edge.weight << ") ";

    }
    std::cout << std::endl;

    Node* root = constructHierarchy(edges); //TODO: Get tree sizes annotated
    assign_node_sizes(root);

    printTree(*root);
    Dc_hdbscan tree(k, k);
    tree.fit(data, n, dim, k);

    std::vector<int> labels = tree.labels_;
    printLabels(labels);

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
    //printTree(*root);



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

    //print_annotations(res);

}





int main() {
    //test_k_centroids();
    //test_hdbscan();
    //test_parallel_hdbscan();
    //simpleTest();
    //test_mlpack();
    test_mst();
    return 0;
}
