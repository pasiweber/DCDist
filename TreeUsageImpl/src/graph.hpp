#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <key_structs.hpp>
typedef struct MSTEdge
{
    unsigned long long src;
    unsigned long long dst;
    double weight;
} MSTEdge;

std::vector<MSTEdge> calc_mst(unsigned long long V, double *data, int minClusterPoints);

Node* constructHierarchy(std::vector<MSTEdge> edges);



#endif