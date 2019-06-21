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
    degreeTemp.resize(_n);
    spanningTree.resize(_n - 1);
    disjointSets.initialize(_n);
    lagrangeanMultipliers.resize(_n);
}

void DCMSTPLagrangean::sortEdges() {
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
}

void DCMSTPLagrangean::kruskalx() {
    int edgesSpanTree = 0;
    auto currentEdge = edges.begin();

    disjointSets.clean();
    std::fill(degreeTemp.begin(), degreeTemp.end(), 0);

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
            for (int w = 0; w < getNumVertices() && (treeUSaturated || treeVSaturated); ++w) {
                if (degreeTemp[w] < degrees[w]) {
                    if (treeUSaturated && disjointSets.find(u) == disjointSets.find(w)) {
                        treeUSaturated = false;
                    } else if (treeVSaturated && disjointSets.find(v) == disjointSets.find(w)) {
                        treeVSaturated = false;
                    }
                }
            }

            /* components in (E_1 U {e_k}) are non saturated then */
            if (treeUSaturated == false && treeUSaturated == false) {
                disjointSets.unionSets(u, v);

                spanningTree[edgesSpanTree] = *currentEdge;
                ++edgesSpanTree;
            } else {
                --degreeTemp[u];
                --degreeTemp[v];
            }
        }

        currentEdge++;
    }
}

void DCMSTPLagrangean::kruskal() {
    int edgesSpanTree = 0;
    auto currentEdge = edges.begin();

    disjointSets.clean();

    /* while a tree was not formed or not all the edges were used */
    while (edgesSpanTree < getNumVertices() - 1 && currentEdge != edges.end()) {
        int u = currentEdge->u;
        int v = currentEdge->v;

        if (disjointSets.find(u) != disjointSets.find(v)) {
            disjointSets.unionSets(u, v);

            spanningTree[edgesSpanTree] = *currentEdge;
            ++edgesSpanTree;
        }

        currentEdge++;
    }
}

void DCMSTPLagrangean::updateSubgradLowerVals() {
    z_lb = 0.0f;
    std::fill(subgradient.begin(), subgradient.end(), 0);

    for (int i = 0; i < getNumVertices(); ++i) {
        /* Step 3: Updates subgradient */
        subgradient[i] -= 1 * degrees[i];

        if (i < getNumVertices() - 1) {
            /* Step 3: Updates subgradient */
            subgradient[spanningTree[i].u] += 1;
            subgradient[spanningTree[i].v] += 1;

            /* Step 4: Calculates current lower bound */
            Edge &e = spanningTree[i];

            z_lb += (e.w + lagrangeanMultipliers[e.u] + lagrangeanMultipliers[e.v]);
        }

        /* Step 4: Calculates current lower bound */
        z_lb -= lagrangeanMultipliers[i] * degrees[i];
    }

    /* calculates 2-norm of subgradient vector */
    subgradientNorm = std::inner_product(subgradient.begin(), subgradient.end(), subgradient.begin(), 0.0f);
}

void DCMSTPLagrangean::calculateUb() {
    sortEdges();
    kruskalx();

    z_ub = 0;

    for (int i = 0; i < getNumVertices(); ++i) {
        z_ub += spanningTree[i].w;
    }
}

void DCMSTPLagrangean::solve() {
    int iters = 0;
    float beta = 0.0f;
    float alpha = 2.0f;
    int maxIters = 1000;
    int min_z_ub = std::numeric_limits<int>::max();
    float max_z_lb = std::numeric_limits<float>::min();

    std::fill(lagrangeanMultipliers.begin(), lagrangeanMultipliers.end(), 0);

    /* Step 1: Finds a first valid solution to the DCMSTP */
    calculateUb();

    min_z_ub = std::min(min_z_ub, z_ub);

    while (iters < maxIters && GET_TIME(initialTime, clock()) < limitTime) {
        /* Step 2: Solves MSTP with lagrangean multipliers added to the edges */
        sortEdges();
        kruskal();

        /* Step 3 and 4 */
        updateSubgradLowerVals();

        max_z_lb = std::max(max_z_lb, z_lb);

        calculateUb();

        min_z_ub = std::min(min_z_ub, z_ub);

        if (iters % 100 == 0)
            printf("z_lb (%d) : %f  --  %d\n", iters, z_lb, z_ub);

        if (subgradientNorm < 1e-10) break;

        /* Step 5: Calculates step size */
        float stepSize = alpha * ((1.0f + beta) * min_z_ub - z_lb) / subgradientNorm;

        /* Step 6: Update lagrangean multipliers */
        for (int i = 0; i < getNumVertices(); ++i) {
            lagrangeanMultipliers[i] = std::max(0.0f, lagrangeanMultipliers[i] + stepSize * subgradient[i]);
        }

        ++iters;
    }

    printf("max lb, ub: %f -- %d\n", max_z_lb, z_ub);
}