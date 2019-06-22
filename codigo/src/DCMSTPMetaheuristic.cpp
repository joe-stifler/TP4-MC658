/**
 * File: DCMSTPMetaheuristic.cpp
 *
 * Discipline: MC658
 * PED: Natanael Ramos
 * Professor: Cid C. de Souza
 * Data of creation: June 8, 2019
 * Author (RA 176665): Jose Ribeiro Neto <j176665@dac.unicamp.br>
 * Author (RA 171119): Felipe Lopes De Mello <f171119@dac.unicamp.br>
 *
 **/

#include <DCMSTPMetaheuristic.h>

#define POP_SIZE 500 // size of initial population

DCMSTPMetaheuristic::DCMSTPMetaheuristic(int _n, int _limitTime, clock_t _initialTime) : DCMSTP(_n, _limitTime, _initialTime) {
    population.resize(POP_SIZE);
    degreeTemp.resize(_n);
    disjointSets.initialize(_n);
}
  
// A utility function that creates a graph of  
// V vertices 
Graph* createGraph(int V) 
{ 
    Graph* graph = new Graph; 
    graph->V = V; 
  
    // Create an array of sets representing 
    // adjacency lists. Size of the array will be V 
    graph->adjList = new std:: unordered_set<int>[V]; 
  
    return graph; 
} 
  
// Adds an edge to an undirected graph 
void addEdgeSet(Graph* graph, int src, int dest) 
{ 
    // Add an edge from src to dest. A new 
    // element is inserted to the adjacent 
    // list of src. 
    graph->adjList[src].insert(dest); 
  
    // Since graph is undirected, add an edge 
    // from dest to src also 
    graph->adjList[dest].insert(src); 
} 
  
// A utility function to print the adjacency 
// list representation of graph 
/*void printGraph(Graph* graph) 
{ 
    for (int i = 0; i < graph->V; ++i) { 
        std::unordered_set<int> lst = graph->adjList[i]; 
        cout << endl << "Adjacency list of vertex "
             << i << endl; 
  
        for (auto itr = lst.begin(); itr != lst.end(); ++itr) 
            cout << *itr << " "; 
        cout << endl; 
    } 
} */
  
// Searches for a given edge in the graph 
bool searchEdge(Graph* graph, int src, int dest) 
{ 
    auto itr = graph->adjList[src].find(dest); 
    if (itr == graph->adjList[src].end()) 
        return false;
    return true;
}

void DCMSTPMetaheuristic::RandomKruskalX(Graph* individual) {

	std::vector<Edge> remainingEdges = edges; /* Edges not considered to be in the tree */

	Edge ek;   // Edge being considered
	int idx;   // Index of the edge being considered
	
	int edgesSpanTree = 0;
    auto currentEdge = edges.begin();

    disjointSets.clean();
    std::fill(degreeTemp.begin(), degreeTemp.end(), 0);
    
	/* initialize random seed: */
 	srand (time(NULL));
 	
    /* while a tree was not formed or not all the edges were used */
    while (edgesSpanTree < getNumVertices() - 1) {

 		/* random number between 0 and number of edges not considered */
		idx = rand() % (remainingEdges.size()-1);

		/* consider edge idx */
		ek = remainingEdges[idx];
		
        int u = ek.u;
        int v = ek.v;

		/* erase the considered edge */
  		remainingEdges.erase(remainingEdges.begin()+idx);

        if (degreeTemp[u] < degrees[u] && degreeTemp[v] < degrees[v]
                    && disjointSets.find(u) != disjointSets.find(v)) {

            ++degreeTemp[u];
            ++degreeTemp[v];
			
            bool treeUSaturated = true;
            bool treeVSaturated = true;

            /* Verify if the tree related with U and V are both nonsaturated */
            for (int w = 0; w < getNumVertices(); ++w) {
                if (degreeTemp[w] < degrees[w]) {
                    if (disjointSets.find(u) == disjointSets.find(w)) {
                        treeUSaturated = false;
                    }

                    if (disjointSets.find(v) == disjointSets.find(w)) {
                        treeVSaturated = false;
                    }
                }
            }

            /* components in (E_1 U {e_k}) are non saturated then */
            if (edgesSpanTree + 1 == getNumVertices() - 1 || treeUSaturated == false || treeVSaturated == false) {
                disjointSets.unionSets(u, v);

                addEdgeSet(individual, u, v);
                ++edgesSpanTree;
            } else {
                --degreeTemp[u];
                --degreeTemp[v];
            }
        }
    }
}

void DCMSTPMetaheuristic::initializePopulation() {
	int i;
	
	struct Graph* individual = createGraph(getNumVertices());
    /* Generate POP_SIZE random individuals */
	for(i = 0; i < POP_SIZE; i++) {
		RandomKruskalX(individual);
		population[i] = *individual;
	}
}

void DCMSTPMetaheuristic::solve() {

    int iters = 0;

    maxIters = 10000;
    
    /* Step 1: Initialize population with random DCMST */
	initializePopulation();

	//printf("teste: (%d)\n", population[0].V);
	
	/* Step 2: Recombine and mutate until stop criteria */
    while (iters < maxIters && GET_TIME(initialTime, clock()) < limitTime) {
        /* Step 2-A: Select individuals that will recombine */
        
        /* Step 2-B: Create new population with the parents chosen in 2-A*/
        
        /* Step 2-C: Mutate individuals resulted in 2-B */
        
        /* Step 2-D: Replace worst individuals with the best of resulted in 2-C */

        ++iters;
    }
}
