#ifndef TP4_EDGE_SET_H
#define TP4_EDGE_SET_H

#include <vector>
#include <Edge.h>

struct Chromosome {
    int fitness;
    std::vector<Edge> spanningTree;

    Chromosome() {
        fitness = 0;
    }

    Chromosome(int _n) {
        fitness = 0;
        spanningTree.resize(_n);
    }

    // Adds an edge to spanning tree
    void addEdge(Edge e) {
        fitness += e.w;
        spanningTree.push_back(e);
    }
};

void addEdgeChr(Chromosome*, int, int);

#endif //TP4_EDGE_SET_H
