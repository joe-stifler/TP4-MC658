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

#define POP_SIZE 500    // size of initial population
#define CROSS_PROB 0.8 // crossover probability

DCMSTPMetaheuristic::DCMSTPMetaheuristic(int _n, int _limitTime, clock_t _initialTime) : DCMSTP(_n, _limitTime, _initialTime) {
    degreeTemp.resize(_n);
    disjointSets.initialize(_n);
    population.resize(POP_SIZE);
}

void DCMSTPMetaheuristic::RandomKruskalX(Chromosome &individual, bool reduxEdges) {
	/* Inicializacao da arvore */
    int edgesSpanTree = 0;

    /* Initialize arrays */
    for (int i = 0; i < getNumVertices(); ++i) {
        disjointSets.rank[i] = 0;
        disjointSets.parent[i] = i;

        degreeTemp[i] = 0;
    }

    int kEdge = 0;
    auto currentEdge = edges.begin();

    /* while a tree was not formed or not all the edges were used */
    while (edgesSpanTree < getNumVertices() - 1 && currentEdge != edges.end()) {
        int u = currentEdge->u;
        int v = currentEdge->v;

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
                individual.addEdge(*currentEdge);

                ++edgesSpanTree;
            } else {
                --degreeTemp[u];
                --degreeTemp[v];
            }
        }

        ++kEdge;
        currentEdge++;
    }

    if (reduxEdges) {
        edges.resize(std::min(3 * kEdge, (int) edges.size()));
    }
}

// Tempo medio de geracao de 1 individuo eh 0.008116 segundos no felipe pc na instancia tb2ct500_3
void DCMSTPMetaheuristic::initializePopulation(bool heuristicInicialization) {
    
    std::sort(edges.begin(), edges.end(), [&] (Edge &e1, Edge &e2) -> bool {
        return e1.w < e2.w;
    });
	            
	if(heuristicInicialization == false) {

		/* Generate POP_SIZE random individuals */
		for(int i = 0; i < POP_SIZE; i++) {
		    if (i == 0) {
	            RandomKruskalX(population[i], true);
		    } else {
	            std::random_shuffle(edges.begin(), edges.end());
	            RandomKruskalX(population[i], false);
	        }
		}
		
		std::sort(edges.begin(), edges.end(), [&] (Edge &e1, Edge &e2) -> bool {
        	return e1.w < e2.w;
    	});
	}
	
	else {

		float alpha = 1.5/(500/POP_SIZE); //paper: alpha = 1.5 for pop_size = 500
		int k;
		int n = getNumVertices();
		std::default_random_engine generator;
		std::uniform_int_distribution<int> dist(0, getNumEdges() - 1);
				
	    RandomKruskalX(population[0], true);
	            
		for(int i = 1; i < POP_SIZE; i++) {
	
			k = (int) alpha*i*n/POP_SIZE; //
			
			/* Shuffle k cheapest edges */
			while(k > 0) {
				int randEdge = dist(generator);
				std::swap(edges[k],edges[randEdge]);
				k--;
			}
	        RandomKruskalX(population[i], false);
	        
	        std::sort(edges.begin(), edges.end(), [&] (Edge &e1, Edge &e2) -> bool {
	        	return e1.w < e2.w;
	   		});
	    }
	}
}

/*
Pegar todas as arestas em comum entre os pais.
Pegar duas arestas aleat�rias que pertencem a um e n�o a outro. Dessas duas escolher a de menor peso. Repetir esse passo at� onde for poss�vel.
Caso necess�rio, completar a �rvore com arestas que n�o pertencem a nenhum pai. Novamente este passo pode ser bem dif�cil para nosso caso
(o grafo dele � completo e todos as restri��es de grau s�o maiores que 2).
*/
Chromosome DCMSTPMetaheuristic::crossover(Chromosome &parent1, Chromosome &parent2) {
	Chromosome child;
    //std::default_random_engine generator;
    //std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    int edgesSpanTree = 0;

    std::vector<Edge> edgesPar1;
    std::vector<Edge> edgesPar2;

    std::set<std::pair<int, int>> usedEdges;
    std::set<std::pair<int, int>> usedEdgesPar1;
    std::set<std::pair<int, int>> usedEdgesPar2;

    /* Initialize arrays */
    for (int i = 0; i < getNumVertices(); ++i) {
        disjointSets.rank[i] = 0;
        disjointSets.parent[i] = i;

        degreeTemp[i] = 0;

        if (i < getNumVertices() - 1) {
            Edge e = parent1.spanningTree[i];
            usedEdgesPar1.insert(std::make_pair(e.u, e.v));

            e = parent2.spanningTree[i];
            usedEdgesPar2.insert(std::make_pair(e.u, e.v));
        }
    }

    /* Get all commom edges between spanningTree1 and spanningTree2 */
    for (int i = 0; i < getNumVertices() - 1; ++i) {
        Edge e = parent2.spanningTree[i];

        if (usedEdgesPar1.find(std::make_pair(e.u, e.v)) != usedEdgesPar1.end()) {
            /* insert edge belonging to spanningTree1 and spanningTree2*/
            child.addEdge(e);
            disjointSets.unionSets(e.u, e.v);
            usedEdges.insert(std::make_pair(e.u, e.v));

            ++degreeTemp[e.u];
            ++degreeTemp[e.v];

            ++edgesSpanTree;
        } else {
            edgesPar2.push_back(e); /* insert edge belonging only to spanningTree2 */
        }

        e = parent1.spanningTree[i];

        if (usedEdgesPar2.find(std::make_pair(e.u, e.v)) == usedEdgesPar2.end()) {
            edgesPar1.push_back(e); /* insert edge belonging only to spanningTree1 */
        }
    }

    std::random_shuffle(edgesPar1.begin(), edgesPar1.end());
    std::random_shuffle(edgesPar2.begin(), edgesPar2.end());

    auto currEdgePar1 = edgesPar1.begin();
    auto currEdgePar2 = edgesPar2.begin();

    /* Get randomly edges from spanningTree1 and spanningTree2 */
    while (edgesSpanTree < getNumVertices() - 1
                && (currEdgePar1 != edgesPar1.end()
                    || currEdgePar2 != edgesPar2.end())) {

        Edge e;

        if (currEdgePar1 != edgesPar1.end()
            && currEdgePar2 != edgesPar2.end()) {

            if (currEdgePar1->w < currEdgePar2->w) {
                e = *currEdgePar1;
            } else e = *currEdgePar2;

            ++currEdgePar1;
            ++currEdgePar2;
        } else if (currEdgePar1 != edgesPar1.end()) {
            e = *currEdgePar1;
            currEdgePar1++;
        } else {
            e = *currEdgePar2;
            currEdgePar2++;
        }

        /* Verify if edge was not already used and verify if degree constraints are satisfied */
        if (degreeTemp[e.u] < degrees[e.u] && degreeTemp[e.v] < degrees[e.v] &&
                                        disjointSets.find(e.u) != disjointSets.find(e.v)) {

            ++degreeTemp[e.u];
            ++degreeTemp[e.v];

            bool treeUSaturated = true;
            bool treeVSaturated = true;

            /* Verify if the tree related with U and V are both nonsaturated */
            for (int z = 0; z < getNumVertices(); ++z) {
                if (degreeTemp[z] < degrees[z]) {
                    if (disjointSets.find(e.u) == disjointSets.find(z)) {
                        treeUSaturated = false;
                    }

                    if (disjointSets.find(e.v) == disjointSets.find(z)) {
                        treeVSaturated = false;
                    }
                }
            }

            /* components in (E_1 U {e_k}) are non saturated then */
            if (edgesSpanTree + 1 == getNumVertices() - 1 || treeUSaturated == false || treeVSaturated == false) {
                child.addEdge(e);
                disjointSets.unionSets(e.u, e.v);
                usedEdges.insert(std::make_pair(e.u, e.v));

                ++edgesSpanTree;
            } else {
                --degreeTemp[e.u];
                --degreeTemp[e.v];
            }
        }
    }

    currEdgePar1 = edges.begin();
    while (edgesSpanTree < getNumVertices() - 1 && currEdgePar1 != edges.end()) {
        Edge &e = *currEdgePar1;

        if (usedEdges.find(std::make_pair(e.u, e.v)) == usedEdges.end()) {
            if (degreeTemp[e.u] < degrees[e.u] &&
                degreeTemp[e.v] < degrees[e.v] &&
                disjointSets.find(e.u) != disjointSets.find(e.v)) {

                ++degreeTemp[e.u];
                ++degreeTemp[e.v];

                bool treeUSaturated = true;
                bool treeVSaturated = true;

                /* Verify if the tree related with U and V are both nonsaturated */
                for (int z = 0; z < getNumVertices(); ++z) {
                    if (degreeTemp[z] < degrees[z]) {
                        if (disjointSets.find(e.u) == disjointSets.find(z)) {
                            treeUSaturated = false;
                        }

                        if (disjointSets.find(e.v) == disjointSets.find(z)) {
                            treeVSaturated = false;
                        }
                    }
                }

                /* components in (E_1 U {e_k}) are non saturated then */
                if (edgesSpanTree + 1 == getNumVertices() - 1 || treeUSaturated == false || treeVSaturated == false) {
                    child.addEdge(e);
                    disjointSets.unionSets(e.u, e.v);
                    usedEdges.insert(std::make_pair(e.u, e.v));

                    ++edgesSpanTree;
                } else {
                    --degreeTemp[e.u];
                    --degreeTemp[e.v];
                }
            }
        }

        currEdgePar1++;
    }

	return child;
}

/*
Muta apenas 1 indiv�duo (o tamanho da popula��o � 500).
A muta��o � feita segundo a p�gina 30 que s�o c�lculos pesados talvez fora do escopo deste trabalho.
Uma op��o � ordenar as arestas que n�o est�o na �rvore, pegar a de menor custo que � vi�vel e inclu�-la na �rvore, criando um ciclo.
Remover deste ciclo uma aresta, criando uma �rvore.
*/
void DCMSTPMetaheuristic::mutate(Chromosome &child) {
    std::default_random_engine generator;
    std::set<std::pair<int, int>> usedEdges;
    std::uniform_real_distribution<float> dist2(0, 1.0f);
    std::uniform_int_distribution<int> dist(0, getNumVertices() - 2);

    for(int i = 0; i < getNumVertices() - 1; ++i) {
        Edge &e = child.spanningTree[i];
        usedEdges.insert(std::make_pair(e.u, e.v));
    }

    int chosenEdge = dist(generator);

    Edge &e = child.spanningTree[chosenEdge];

    disjointSets.clean();
    std::fill(degreeTemp.begin(), degreeTemp.end(), 0);

    /* Find both trees and compute vertices degree for next step */
    for(int j = 0; j < getNumVertices() - 1; ++j) {
        if (chosenEdge != j) {
            Edge ej = child.spanningTree[j];
            degreeTemp[ej.u]++;
            degreeTemp[ej.v]++;

            disjointSets.unionSets(ej.u, ej.v);
        }
    }

    Edge newEdge = e;
    Edge oldEdge = e;

    newEdge.w = std::numeric_limits<int>::max();

    /* Select minimun cost edge that connects both
     * trees without violating degree constraints */
    auto currentEdge = edges.begin();
    for(int j = 0; j < getNumEdges(); ++j) {
        /* If selected edge is lighter than the removed one, switch */
        if (currentEdge->w < newEdge.w
            && degreeTemp[currentEdge->u] < degrees[currentEdge->u]
            && degreeTemp[currentEdge->v] < degrees[currentEdge->v]
            && disjointSets.find(currentEdge->u) != disjointSets.find(currentEdge->v)
            && usedEdges.find(std::make_pair(currentEdge->u, currentEdge->v)) == usedEdges.end()) {

            //if (dist2(generator) < CROSS_PROB) {
            newEdge = *currentEdge;
            //}
        }

        currentEdge++;
    }

    if (newEdge.u != oldEdge.u || newEdge.v != oldEdge.v) {
        z_ub += newEdge.w;
        z_ub -= oldEdge.w;

        child.spanningTree[chosenEdge] = newEdge;
        child.fitness += newEdge.w;
        child.fitness -= oldEdge.w;

        usedEdges.erase(std::make_pair(oldEdge.u, oldEdge.v));
        usedEdges.insert(std::make_pair(newEdge.u, newEdge.v));
    }
}

void DCMSTPMetaheuristic::solve() {
    int iters = 0;
    int maxIters = 100000;
    std::default_random_engine generator;
    int bestPrimal = std::numeric_limits<int>::max();
    std::uniform_int_distribution<int> dist(0, POP_SIZE - 1);

    /* Step 1: Initialize population with random DCMST */
    bool heuristicInicialization = true;
	initializePopulation(heuristicInicialization);

	for(int i = 0; i < POP_SIZE; i++) {
        bestPrimal = std::min(bestPrimal, population[i].fitness);
		//printf("Fitness do individuo %d: %d\n", i, population[i].fitness);
	}
	printf("Best primal da geracao %d: %d\n", iters, bestPrimal);	
	std::vector<Chromosome> child;
	child.resize(POP_SIZE*CROSS_PROB);
	
	/* Step 2: Recombine and mutate until stop criteria */
	/* Better stop criteria propabably is when best individual isnt growing for too long */
    while (iters < maxIters && GET_TIME(initialTime, clock()) < limitTime) {

		/* Sort generation by fitness */
		std::sort(population.begin(), population.end(), [&] (Chromosome &p1, Chromosome &p2) -> bool {
        	return p1.fitness < p2.fitness;
    	});
    	
        /* Step 2-A: Select individuals that will recombine */
        /* New population should have size iguals to CROSS_PROB*POP_SIZE */
        for(int i = 0; i < POP_SIZE*CROSS_PROB; i++) {
        	Chromosome parent1;
        	Chromosome parent2;
            int indexRand1 = dist(generator);
            int indexRand2;

            do indexRand2 = dist(generator);
            while(indexRand2 == indexRand1);
            
			if (population[indexRand1].fitness < population[indexRand2].fitness)
                parent1 = population[indexRand1];
            else parent1 = population[indexRand2];
			
			indexRand1 = dist(generator);
			do indexRand2 = dist(generator);
            while(indexRand2 == indexRand1);
            
			if (population[indexRand1].fitness < population[indexRand2].fitness)
                parent2 = population[indexRand1];
            else parent2 = population[indexRand2];
            
			/* Step 2-B: Create new population with the parents chosen in 2-A*/
            child[i] = crossover(population[indexRand1], population[indexRand2]);

            /*if (population[indexRand1].fitness < population[indexRand2].fitness)
                child[i] = crossover(population[indexRand1], population[indexRand2]);
            else child[i] = crossover(population[indexRand2], population[indexRand1]);*/

            /* Step 2-C: Mutate individuals resulted in 2-B */
            /* Mutate only 1 individual */
            if(i == 0)
            	mutate(child[i]);

            //bestPrimal = std::min(bestPrimal, child.fitness);

			/* Step 2-D: Replace worst population with the children */
			/* If it was done here children could reproduce with parents in next iteration */
			//population[i+(POP_SIZE-CROSS_PROB*POP_SIZE)] = child[i];
			
            /*printf("\tBest primal: %d --> %d\n", bestPrimal, child.fitness);
            if (child.fitness < population[indexRand1].fitness && child.fitness < population[indexRand2].fitness) {
                if (population[indexRand2].fitness > population[indexRand1].fitness) {
                    population[indexRand2] = child;
                } else population[indexRand1] = child;
            }*/
        }
        
        /* Step 2-D: Replace worst population with the children */
		for(int i = 0; i < POP_SIZE*CROSS_PROB; i++) {
			/* Should only replace if it does not duplicate an existing solution */
			population[i+(POP_SIZE-CROSS_PROB*POP_SIZE)] = child[i];
			bestPrimal = std::min(bestPrimal, child[i].fitness);
		}

		//if(iters % 100 == 0) {
		
			//for(int i = 0; i < POP_SIZE; i++) {
        		//bestPrimal = std::min(bestPrimal, population[i].fitness);
				//printf("Fitness do individuo %d: %d\n", i, population[i].fitness);
			//}
			printf("Best primal da geracao %d: %d\n", iters, bestPrimal);		
		//}
			
        ++iters;
    }

    printf("Best primal: %d\n", bestPrimal);
}
