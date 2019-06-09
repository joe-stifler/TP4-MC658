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

#ifndef TP4_DCMSTPLAGRANGEAN_H
#define TP4_DCMSTPLAGRANGEAN_H

#include <cmath>
#include <memory>
#include <vector>
#include <algorithm>

#include <DCMSTP.h>
#include <DisjointSets.h>

class DCMSTPLagrangean : public DCMSTP {
public:
    /*
     * Simple constructor
     *
     * Params:
     *     _n: total number of vertices in the graph
     * */
    DCMSTPLagrangean(int _n, int _limitTime, clock_t _initialTime);

    void solveMstp();
    void solve() override;

private:
    int z_ub;
    float z_lb;
    int maxIters = 2000;
    float subgradientNorm;
    DisjointSets disjointSets;
    std::vector<int> subgradient;
    std::vector<Edge> spanningTree;
    std::vector<float> lagrangeanMultipliers;

    void calculateUb();
};

#endif // TP4_DCMSTPLAGRANGEAN_H
