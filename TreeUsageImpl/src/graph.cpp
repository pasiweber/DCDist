#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <graph_hdb_c.hpp>
#include <graph.hpp>
#include <iostream>
// gcc -shared -o graph.so graph.c -O3 -march=native -lm -fopenmp -fPIC

typedef struct Subset
{
    unsigned long long parent;
    unsigned long long rank;
} Subset;

typedef struct Edge
{
    unsigned long long src;
    unsigned long long dst;
} Edge;





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

// int *constructHierarchy(struct MSTEdge *edges, unsigned long long numEdges, int minClusterPoints)
// {
//     qsort(edges, numEdges, sizeof(struct MSTEdge), compareEdges);
//     unsigned long long numPoints = numEdges + 1;

std::vector<MSTEdge> calc_mst(unsigned long long V, double *data, int minClusterPoints)
{
    //TODO: Make the output MST a std::vector. 
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
    std::cout << "first edge: " << mstEdges[0].weight << std::endl;
    std::vector<MSTEdge> result(mstEdges, mstEdges + V-1); 

    return result;
}
