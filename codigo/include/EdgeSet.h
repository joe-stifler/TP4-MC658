#ifndef TP4_EDGE_SET_H
#define TP4_EDGE_SET_H

#include <vector>
#include <Edge.h>

struct Chromosome { 
    std::vector<Edge> spanningTree;
    int fitness;
}; 

Chromosome* createChr(int V);

void addEdgeChr(Chromosome*, int, int);

#endif //TP4_EDGE_SET_H
