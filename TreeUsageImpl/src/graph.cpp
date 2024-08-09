#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>

#include <graph_hdb_c.hpp>
#include <graph.hpp>
#include <iostream>
#include <algorithm>
#include <key_structs.hpp>
// gcc -shared -o graph.so graph.c -O3 -march=native -lm -fopenmp -fPIC

typedef struct Subset
{
    unsigned long long parent;
    unsigned long long rank;
    Node* root_node;
} Subset;

typedef struct Edge
{
    unsigned long long src;
    unsigned long long dst;
} Edge;



//This does path compression while finding the root, so that next time i is queried this will go a lot faster.
unsigned long long findParentIdx(Subset *subsets, unsigned long long i)
{
    while (subsets[i].parent != i)
    {
        subsets[i].parent = subsets[subsets[i].parent].parent;
        i = subsets[i].parent;
    }
    return subsets[i].parent;
}

void unifySets(Subset *subsets, unsigned long long x, unsigned long long y)
{
    unsigned long long xRoot = findParentIdx(subsets, x);
    unsigned long long yRoot = findParentIdx(subsets, y);

    if (xRoot != yRoot)
    {
        if (subsets[xRoot].rank > subsets[yRoot].rank)
        {
            subsets[yRoot].parent = xRoot;
        }
        else if (subsets[yRoot].rank > subsets[xRoot].rank)
        {
            subsets[xRoot].parent = yRoot;
        }
        else
        {
            subsets[yRoot].parent = xRoot;
            subsets[xRoot].rank++;
        }
    }
}

int compareEdges(const void *a, const void *b)
{
    return ((struct MSTEdge *)a)->weight > ((struct MSTEdge *)b)->weight;
}


bool compareByCost(const MSTEdge e1, const MSTEdge e2){
    return e1.weight <= e2.weight;
}

// int *constructHierarchy(struct MSTEdge *edges, unsigned long long numEdges, int minClusterPoints)
// {
//     qsort(edges, numEdges, sizeof(struct MSTEdge), compareEdges);
//     unsigned long long numPoints = numEdges + 1;

std::vector<MSTEdge> calc_mst(unsigned long long V, double *data, int minClusterPoints)
{
    //TODO: Make the data (mut reach dists) something it computes itself, essentially the function should take the points as input. 
    
    Subset *subsets = (Subset *) malloc(V * sizeof(Subset));
    Edge *nearest = (Edge *) calloc(V, sizeof(Edge));
    unsigned long long *rootParent = (unsigned long long *)calloc(V, sizeof(unsigned long long));

    MSTEdge *mstEdges = (MSTEdge *) calloc(V - 1, sizeof(MSTEdge));

    // initialize each element as a separate disjoint set
    // #pragma omp parallel for schedule(static, 32)
    for (unsigned long long i = 0; i < V; i++)
    {
        subsets[i].parent = i;
        subsets[i].rank = 0;
        rootParent[i] = -1;
        nearest[i].src = -1;
    }

    unsigned long long noOfTrees = V;

    unsigned long long v = 0;

    while (noOfTrees > 1)
    {

        for (unsigned long long i = 0; i < V; i++)
        {
            for (unsigned long long j = i + 1; j < V; j++)
            {
                if (data[i * V + j] != 0)
                {
                    // if 2 elements have the same parent at a point, then they will always have the same parent so
                    // no need to compute
                    if (rootParent[i] != -1 &&
                        (rootParent[i] == rootParent[j]))
                    {
                        continue;
                    }

                    // find parent
                    unsigned long long root1 = findParentIdx(subsets, i);
                    unsigned long long root2 = findParentIdx(subsets, j);

                    // ignore if the same parent (same disjoint set)
                    if (root1 == root2)
                    {
                        rootParent[i] = root1;
                        rootParent[j] = root2;
                        continue;
                    }

                    if (nearest[root1].src == -1 || data[nearest[root1].src * V + nearest[root1].dst] > data[i * V + j])
                    {
                        nearest[root1].src = i;
                        nearest[root1].dst = j;
                    }
                    if (nearest[root2].src == -1 || data[nearest[root2].src * V + nearest[root2].dst] > data[i * V + j])
                    {
                        nearest[root2].src = i;
                        nearest[root2].dst = j;
                    }
                }
            }
        }

        for (unsigned long long i = 0; i < V; i++)
        {
            if (nearest[i].src == -1)
            {
                continue;
            }

            unsigned long long root1 = findParentIdx(subsets, nearest[i].src);
            unsigned long long root2 = findParentIdx(subsets, nearest[i].dst);

            if (root1 == root2)
            {
                nearest[i].src = -1;
                continue;
            }

            mstEdges[v].src = nearest[i].src;
            mstEdges[v].dst = nearest[i].dst;
            mstEdges[v].weight = data[nearest[i].src * V + nearest[i].dst];
            nearest[i].src = -1;

            // unify trees/disjoint sets
            unifySets(subsets, root1, root2);
            noOfTrees--;
            v++;
        }
    }
    free(nearest);
    free(subsets);
    free(rootParent);

    //qsort(mstEdges, V-1, sizeof(MSTEdge), compareEdges); //qsort is slower than std::sort, so we should not use it.
    std::vector<MSTEdge> result(mstEdges, mstEdges + V-1); 
    std::sort(result.begin(), result.end(), compareByCost); //TODO: Reverse the order.
    
    free(mstEdges); // Remember to free the allocated memory
    return result;
}

//Construct tree from MST via union find with path compression
Node* constructHierarchy(std::vector<MSTEdge> edges)
{
    int numEdges = edges.size();
    unsigned long long numPoints = numEdges + 1;
    struct Subset *subsets = (struct Subset *)malloc(numPoints * sizeof(struct Subset)); //We need one subset struct for each point
    //These keep track of the parent and the rank, rank being how deep the point is relative to the parent in the current tree.
    std::vector<Node*>nodes(numPoints);

    for (int i = 0; i < numPoints; i++) //Note that all ID's are +1, 1-indexed, just to make it easier to look at the tree and prints from here at the same time.
    {
        subsets[i].parent = i;
        subsets[i].rank = 0;
        nodes[i] = new Node{nullptr, 0, i+1}; //i is the id of the point, 0 is the 0-cost that a leaf node has.
        subsets[i].root_node = nodes[i];
    }

    Node* root;
    for (unsigned long long i = 0; i < numEdges; i++)
    //Process all the edges and create internal nodes when needed
    { 
        double weight = edges[i].weight;
        unsigned long long src = edges[i].src;
        unsigned long long dst = edges[i].dst;
        
        // rootSrc and rootDst is just the representative of the connected component each point is currently part of. 
        unsigned long long rootSrc = findParentIdx(subsets, src); //goes to the leaf node id representing the current connected component 
        unsigned long long rootDst = findParentIdx(subsets, dst); //findParentIdx also does path compression

        // std::cout << "edge: " << weight << ", "<< src << "," << dst << std::endl;
        // std::cout << "rootSrc, rootDst: " << rootSrc << ", " << rootDst << std::endl; 

        Node* rootSrc_node = subsets[rootSrc].root_node;
        Node* rootDst_node = subsets[rootDst].root_node;
        double src_weight = rootSrc_node->cost;
        double dst_weight = rootDst_node->cost;
        // std::cout << "src weight: " << src_weight << std::endl;
        // std::cout << "dst weight: " << dst_weight << std::endl;

        if(weight == src_weight){ //dst cluster should be connected below src cluster.
            rootSrc_node->children.push_back(rootDst_node);
            subsets[rootDst].root_node = rootSrc_node;
            rootDst_node->parent = rootSrc_node;

            // std::cout << "old node" << std::endl;
        } else if (weight == dst_weight){ //src cluster should be connected below dst cluster.
            rootDst_node->children.push_back(rootSrc_node);
            subsets[rootSrc].root_node = rootDst_node;
            rootSrc_node->parent = rootDst_node;

            // std::cout << "old node" << std::endl;
        } else{ //We need to create a new node.
            // std::cout << "new node" << std::endl;
            Node* new_parent = new Node{nullptr, weight, -1};
            //Update two connected components to point to this
            rootSrc_node->parent = new_parent;
            rootDst_node->parent = new_parent;
            new_parent->children.push_back(rootSrc_node);
            new_parent->children.push_back(rootDst_node);
            subsets[rootSrc].root_node = new_parent;
            subsets[rootDst].root_node = new_parent;
            root = new_parent;
        }
        unifySets(subsets, edges[i].src, edges[i].dst);
    }

    // Clean up memory
    free(subsets);
    return root;
}

