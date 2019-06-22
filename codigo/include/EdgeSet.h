#ifndef TP4_EDGE_SET_H
#define TP4_EDGE_SET_H

#include <unordered_set>

struct Graph { 
    int V; 
    std::unordered_set<int>* adjList; 
}; 

Graph* createGraph(int V);

void addEdgeSet(Graph*, int, int);
  
bool searchEdge(Graph*, int, int);

#endif //TP4_EDGE_SET_H
