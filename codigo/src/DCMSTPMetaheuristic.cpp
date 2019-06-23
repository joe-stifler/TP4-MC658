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

#define POP_SIZE 10    // size of initial population
#define CROSS_PROB 0.8 // crossover probability

DCMSTPMetaheuristic::DCMSTPMetaheuristic(int _n, int _limitTime, clock_t _initialTime) : DCMSTP(_n, _limitTime, _initialTime) {
    population.resize(POP_SIZE);
    sons.resize(CROSS_PROB*POP_SIZE);
    degreeTemp.resize(_n);
    disjointSets.initialize(_n);
}

Chromosome* createChr(int V)
{
    Chromosome* chr = new Chromosome;
    //chr->spanningTree = new std::vector<Edge>;
    chr->spanningTree.resize(V);

  	chr->fitness = 0;

    return chr;
}

// Adds an edge to spanning tree
void addEdgeChr(Chromosome* chr, Edge e)
{
    chr->spanningTree.push_back(e);
    chr->fitness += e.w;
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
  
/* Searches for a given edge in the graph
bool searchEdge(Graph* graph, int src, int dest) 
{
    auto itr = graph->adjList[src].find(dest);
    if (itr == graph->adjList[src].end())
        return false;
    return true;
}*/

void DCMSTPMetaheuristic::RandomKruskalX(Chromosome* individual) {

	std::vector<Edge> remainingEdges = edges; /* Edges not considered to be in the tree */
	std::random_shuffle(remainingEdges.begin(), remainingEdges.end());
	
	Edge ek;   // Edge being considered
	int idx;   // Index of the edge being considered
	
	/* Inicializacao da arvore */
	int edgesSpanTree = 0;
    disjointSets.clean();
    std::fill(degreeTemp.begin(), degreeTemp.end(), 0);
    idx = 0;
    
    /* while a tree was not formed or not all the edges were used */
    while (edgesSpanTree < getNumVertices() - 1 && idx < remainingEdges.size()) {
		
		/* consider edge idx */
		ek = remainingEdges[idx];
        int u = ek.u;
        int v = ek.v;

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
                addEdgeChr(individual, ek);
                ++edgesSpanTree;
            } else {
                --degreeTemp[u];
                --degreeTemp[v];
            }
        }
		idx++;
    }
}

// Tempo medio de geracao de 1 individuo eh 0.008116 segundos no felipe pc na instancia tb2ct500_3
void DCMSTPMetaheuristic::initializePopulation() {
	int i;
	
	struct Chromosome* individual = createChr(getNumVertices());
    /* Generate POP_SIZE random individuals */
	for(i = 0; i < POP_SIZE; i++) {
		RandomKruskalX(individual);
		population[i] = *individual;
	}
}

/*
Pegar todas as arestas em comum entre os pais.
Pegar duas arestas aleatórias que pertencem a um e não a outro. Dessas duas escolher a de menor peso. Repetir esse passo até onde for possível.
Caso necessário, completar a árvore com arestas que não pertencem a nenhum pai. Novamente este passo pode ser bem difícil para nosso caso
(o grafo dele é completo e todos as restrições de grau são maiores que 2).
*/
Chromosome DCMSTPMetaheuristic::crossover(Chromosome x, Chromosome y) {
	int i;
	
	Chromosome xy;
    
	return xy;
}

/*
Muta apenas 1 indivíduo (o tamanho da população é 500).
A mutação é feita segundo a página 30 que são cálculos pesados talvez fora do escopo deste trabalho.
Uma opção é ordenar as arestas que não estão na árvore, pegar a de menor custo que é viável e incluí-la na árvore, criando um ciclo.
Remover deste ciclo uma aresta, criando uma árvore.
*/
void DCMSTPMetaheuristic::mutate() {
	// mutate sons
	return;
}

void DCMSTPMetaheuristic::solve() {

    int iters = 0, i, idx1, idx2;
    maxIters = 10000;
    std::vector<Chromosome> couple;
    Chromosome x, y;
    
    /* Step 1: Initialize population with random DCMST */
    // Tempo de inicilizar populacao aleatoria de tamanho 500 entre 4 e 5 segundos no felipe pc na instancia tb2ct500_3
	initializePopulation();
	
	for(i = 0; i < POP_SIZE; i++) {
		printf("Fitness do individuo %d: %d\n", i, population[i].fitness);
	}
    
	/* Step 2: Recombine and mutate until stop criteria */
	/* Better stop criteria propabably is when best individual isnt growing for too long */
    while (iters < maxIters && GET_TIME(initialTime, clock()) < limitTime) {
    	
    	/* Note: steps 2-X probably can be combined for optmization */
        /* Step 2-A: Select individuals that will recombine */
        /* New population should have size iguals to CROSS_PROB*POP_SIZE */
        std::random_shuffle(population.begin(), population.end());
        idx1 = 0;
        idx2 = POP_SIZE-1;
        for(i = 0; i < CROSS_PROB*POP_SIZE; i++) {
        	
        	/* Select 2 random individuals and choose the best*/
			if(population[idx1].fitness < population[idx2].fitness)
				x = population[idx1];
            else
            	x = population[idx2];
            
            /* Select 2 random individuals and choose the best*/
			if(population[idx1].fitness < population[idx2].fitness)
				y = population[idx1];
            else
            	y = population[idx2];
            	
			/* Step 2-B: Create new population with the parents chosen in 2-A*/	    	
			sons[i] = crossover(x, y);
			
			++idx1;
			--idx2;
		}
        /* Step 2-C: Mutate individuals resulted in 2-B */
        /* Mutate only 1 individual */
        mutate();
        
        /* Step 2-D: Replace worst individuals with the best of resulted in 2-C */
		/* Sort population by fitness */
		std::sort(population.begin(), population.end(), [&] (Chromosome &p1, Chromosome &p2) -> bool {
        	return p1.fitness < p2.fitness;
    	});
    	for(i = 0; i < CROSS_PROB*POP_SIZE; i++) {
    		population[i+(POP_SIZE-CROSS_PROB*POP_SIZE)] = sons[i];
    	}
    	
        ++iters;
    }
}
