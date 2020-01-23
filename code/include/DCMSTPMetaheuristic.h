/**
 * File: DCMSTPLagrangean.h
 *
 * Discipline: MC658
 * PED: Natanael Ramos
 * Professor: Cid C. de Souza
 * Data of creation: June 8, 2019
 * Author (RA 176665): Jose Ribeiro Neto <j176665@dac.unicamp.br>
 * Author (RA 171119): Felipe Lopes De Mello <f171119@dac.unicamp.br>
 *
 **/

#ifndef TP4_DCMSTPMETAHEURISTIC_H
#define TP4_DCMSTPMETAHEURISTIC_H

#include <set>
#include <cmath>
#include <memory>
#include <vector>
#include <random>
#include <algorithm>

#include <DCMSTP.h>
#include <EdgeSet.h>
#include <DisjointSets.h>

class DCMSTPMetaheuristic : public DCMSTP {
public:
    /*
     * Simple constructor
     *
     * Params:
     *     _n: total number of vertices in the graph
     * */
    DCMSTPMetaheuristic(int _n, int _limitTime, clock_t _initialTime);

    void solve() override;

    int getBestPrimal() override;

    void printSolution(std::string) override;

    void saveBestEdges(std::string) override;

private:
    int POP_SIZE;
    int bestPrimal;
    DisjointSets disjointSets;
    std::vector<int> degreeTemp;
    std::vector<Chromosome> sons;
    std::vector<Chromosome> population;
    std::vector<Edge> bestSpanningTree;

    void mutate(Chromosome &);
    void initializePopulation(bool);
    void RandomKruskalX(Chromosome &, bool);
    Chromosome crossover(Chromosome &, Chromosome &);
    bool testViability(std::vector<Edge> &spanningTreeAux, int &totalCost);
};

#endif // TP4_DCMSTPMETAHEURISTIC_H
