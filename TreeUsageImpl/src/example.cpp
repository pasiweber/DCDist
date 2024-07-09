#include <iostream>
#include <vector>
#include <kcentroids.hpp>
#include <cmath>

// Function to create a new node -
Node* addNode(Node* parent = nullptr, double cost=0.0, int id = -1) {
    Node* newNode = new Node;
    newNode->cost = cost;
    newNode->parent = parent;
    newNode->id = id;

    if(parent != nullptr){
        parent->children.push_back(newNode);
    }

    return newNode;
}

// Function to generate a tree with 11 leaf nodes
Node* generateTree11() {
    Node* root = addNode();  // Create the root node with an arbitrary cost

    // Create internal nodes and leaves
    Node* node1 = addNode(root, 90.0);
    Node* node2 = addNode(root, 80.0);
    Node* node3 = addNode(root, 28.0);

    Node* l1 = addNode(node3, 0.0, 1);
    Node* l2 = addNode(node3, 0.0, 2);
    Node* l3 = addNode(node3, 0.0, 3);


    Node* node4 = addNode(node1, 55.0);
    Node* node5 = addNode(node1, 35.0);

    Node* node6 = addNode(node2, 50.0);
    Node* node7 = addNode(node6, 40.0);

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

// Function to print the tree (for debugging purposes)
void printTree(Node* node, int level = 0) {
    if (node == nullptr) return;
    for (int i = 0; i < level; ++i) std::cout << "  ";
    std::cout << "Node ID: " << node->id << ", Cost: " << node->cost << "\n";
    for (Node* child : node->children) {
        printTree(child, level + 1);
    }
}

double kmedian(double x){
    return x;
}
double kmeans(double x){
    return std::pow(x,2);
}


int main() {
    Node* root = generateTree11();
    std::cout << "Generated Tree:\n";
    std::vector<Annotation*> res = annotate_tree(*root, kmeans);
    printTree(root);

    // Clean up allocated memory (omitted for brevity but should be done in a real program)
    // ...

    return 0;
}
