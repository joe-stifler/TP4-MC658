/**
 * File: DCMSTPLagrangean.cpp
 *
 * Discipline: MC658
 * PED: Natanael Ramos
 * Professor: Cid C. de Souza
 * Data of creation: June 8, 2019
 * Author (RA 176665): Jose Ribeiro Neto <j176665@dac.unicamp.br>
 * Author (RA 171119): Felipe Lopes De Mello <f171119@dac.unicamp.br>
 *
 **/

#include <DCMSTPLagrangean.h>

DCMSTPLagrangean::DCMSTPLagrangean(int _n, int _limitTime, clock_t _initialTime) : DCMSTP(_n, _limitTime, _initialTime) {
    subgradient.resize(_n);
    spannigTree.resize(_n - 1);
    disjointSets.initialize(_n);
    lagrangeanMultipliers.resize(_n);
}

void DCMSTPLagrangean::solveMstp() {
    /* Sorts edge vector by edge weigth
     * (considering here the lagrangean
     * multipliers associated with the edge) */
    std::sort(edges.begin(), edges.end(), [&] (Edge &e1, Edge &e2) -> bool {
        float lamb1_e1 = lagrangeanMultipliers[e1.u];
        float lamb2_e1 = lagrangeanMultipliers[e1.v];

        float lamb1_e2 = lagrangeanMultipliers[e2.u];
        float lamb2_e2 = lagrangeanMultipliers[e2.v];

        return e1.w + lamb1_e1 + lamb2_e1 < e2.w + lamb1_e2 + lamb2_e2;
    });

    int edgesSpanTree = 0;
    auto currentEdge = edges.begin();

    disjointSets.clean();

    while (edgesSpanTree < getNumVertices() - 1 && currentEdge != edges.end()) {
        int u = currentEdge->u;
        int v = currentEdge->v;

        if (disjointSets.find(u) != disjointSets.find(v)) {
            disjointSets.unionSets(u, v);

            spannigTree[edgesSpanTree] = *currentEdge;
            ++edgesSpanTree;
        }

        currentEdge++;
    }
}

void DCMSTPLagrangean::solve() {

    /* Step 1: Initializes PI, Z_UB and lagrangean multipliers */
    int z_ub = 1;

    std::fill(lagrangeanMultipliers.begin(), lagrangeanMultipliers.end(), 0);

    while (GET_TIME(initialTime, clock()) < limitTime) {
        /* Step 2: Solves MSTP with lagrangean multipliers added to the edges */
        solveMstp();

        /* Step 3: Updates subgradient */
        for (int i = 0; i < getNumVertices(); ++i) {
            subgradient[i] = -1 * degrees[i];
        }

        for (int i = 0; i < getNumVertices() - 1; ++i) {
            subgradient[spannigTree[i].u] += 1;
            subgradient[spannigTree[i].v] += 1;
        }

        /* Step 4: Calculates step size */
        float stepSize = 1.0;

        /* Step 5: Update lagrangean multipliers */
        for (int i = 0; i < getNumVertices(); ++i) {
            lagrangeanMultipliers[i] = std::max(0.0f, lagrangeanMultipliers[i] + stepSize * subgradient[i]);
        }
    }
}