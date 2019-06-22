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
    vertexInTree.resize(_n);
    spanningTree.resize(_n - 1);
    disjointSets.initialize(_n);
    disjointDegSets.initialize(_n);
    degSpanningTree.resize(_n - 1);
    spanningTreeTemp.resize(_n - 1);
    lagrangeanMultipliers.resize(_n);
}

void DCMSTPLagrangean::kruskalx(bool shouldRedux) {
    int kDeg = 0;
    int edgesSpanTree = 0;
    int edgesDegSpanTree = 0;
    auto currentEdge = edges.begin();

    z_ub = 0;
    z_lb = 0.0f;

    /* Initialize arrays */
    for (int i = 0; i < getNumVertices(); ++i) {
        disjointSets.rank[i] = 0;
        disjointSets.parent[i] = i;

        disjointDegSets.rank[i] = 0;
        disjointDegSets.parent[i] = i;

        degreeTemp[i] = 0;
        subgradient[i] = 0;

        vertexInTree[i] = false;
    }

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

    /* while a tree was not formed or not all the edges were used */
    while ((edgesSpanTree < getNumVertices() - 1 || edgesDegSpanTree < getNumVertices() - 1) && currentEdge != edges.end()) {
        int u = currentEdge->u;
        int v = currentEdge->v;

        if (edgesSpanTree < getNumVertices() - 1) {
            if (disjointSets.find(u) != disjointSets.find(v)) {
                disjointSets.unionSets(u, v);

                /* Step 3: Updates subgradient */
                subgradient[currentEdge->u] += 1;
                subgradient[currentEdge->v] += 1;

                if (vertexInTree[currentEdge->u] == false) {
                    vertexInTree[currentEdge->u] = true;

                    /* Step 3: Updates subgradient */
                    subgradient[currentEdge->u] -= 1 * degrees[currentEdge->u];

                    /* Step 4: Calculates current lower bound */
                    z_lb -= lagrangeanMultipliers[currentEdge->u] * degrees[currentEdge->u];
                }

                if (vertexInTree[currentEdge->v] == false) {
                    vertexInTree[currentEdge->v] = true;

                    /* Step 3: Updates subgradient */
                    subgradient[currentEdge->v] -= 1 * degrees[currentEdge->v];

                    /* Step 4: Calculates current lower bound */
                    z_lb -= lagrangeanMultipliers[currentEdge->v] * degrees[currentEdge->v];
                }

                /* Step 4: Calculates current lower bound */
                z_lb += (currentEdge->w + lagrangeanMultipliers[currentEdge->u] + lagrangeanMultipliers[currentEdge->v]);

                spanningTree[edgesSpanTree] = *currentEdge;
                ++edgesSpanTree;
            }
        }

        if (edgesDegSpanTree < getNumVertices() - 1) {
            ++kDeg;

            if (degreeTemp[u] < degrees[u] && degreeTemp[v] < degrees[v]
                        && disjointDegSets.find(u) != disjointDegSets.find(v)) {

                ++degreeTemp[u];
                ++degreeTemp[v];

                bool treeUSaturated = true;
                bool treeVSaturated = true;

                /* Verify if the tree related with U and V are both nonsaturated */
                for (int w = 0; w < getNumVertices(); ++w) {
                    if (degreeTemp[w] < degrees[w]) {
                        if (disjointDegSets.find(u) == disjointDegSets.find(w)) {
                            treeUSaturated = false;
                        }

                        if (disjointDegSets.find(v) == disjointDegSets.find(w)) {
                            treeVSaturated = false;
                        }
                    }
                }

                /* components in (E_1 U {e_k}) are non saturated then */
                if (edgesDegSpanTree + 1 == getNumVertices() - 1 || treeUSaturated == false || treeVSaturated == false) {
                    disjointDegSets.unionSets(u, v);

                    degSpanningTree[edgesDegSpanTree] = *currentEdge;
                    ++edgesDegSpanTree;

                    z_ub += currentEdge->w;
                } else {
                    --degreeTemp[u];
                    --degreeTemp[v];
                }
            }
        }

        currentEdge++;
    }

    if (shouldRedux) {
        edges.resize(std::min(2 * kDeg, (int) edges.size()));
    }
}

void DCMSTPLagrangean::improvementProcedure() {
	for(int i = 0; i < getNumVertices(); ++i) {
		Edge &e = degSpanningTree[i];

        disjointDegSets.clean();
        std::fill(degreeTemp.begin(), degreeTemp.end(), 0);

        /* Find both trees and compute vertices degree for next step */
        for(int j = 0; j < getNumVertices() - 1; ++j) {
		    if (i != j) {
                degreeTemp[degSpanningTree[j].u]++;
                degreeTemp[degSpanningTree[j].v]++;

                disjointDegSets.unionSets(degSpanningTree[j].u, degSpanningTree[j].v);
            }
		}

        /* Select minimun cost edge that connects both
         * trees without violating degree constraints */
        auto currentEdge = edges.begin();
        for(int j = 0; j < getNumEdges(); ++j) {
            /* If selected edge is lighter than the removed one, switch */
            if (currentEdge->w < degSpanningTree[i].w &&
			        disjointDegSets.find(currentEdge->u) != disjointDegSets.find(currentEdge->v)
				        && degreeTemp[currentEdge->u] < degrees[currentEdge->u] && degreeTemp[currentEdge->v] < degrees[currentEdge->v]) {

                z_ub += currentEdge->w;
                z_ub -= degSpanningTree[i].w;

                degSpanningTree[i] = *currentEdge;
            }

			currentEdge++;
		}
 	}
}

void DCMSTPLagrangean::solve() {
    int iters = 0;
    float beta = 0.0f;
    float alpha = 2.0f;
    int maxIters = 10000;
    int bestPrimal = std::numeric_limits<int>::max();
    float bestDual = std::numeric_limits<float>::min();

    std::fill(lagrangeanMultipliers.begin(), lagrangeanMultipliers.end(), 0);

    /* Step 1: Finds a first valid solution to the DCMSTP */
    kruskalx(true);
    improvementProcedure();

    bestDual = std::max(bestDual, z_lb);
    bestPrimal = std::min(bestPrimal, z_ub);

    while (iters < maxIters && GET_TIME(initialTime, clock()) < limitTime) {
        /* Step 2: Solves MSTP with lagrangean multipliers added to the edges */
        kruskalx(false);

        if (iters % 50 == 0) {
            improvementProcedure();
            printf("z_lb (%d) : %f  --  %d\n", iters, z_lb, z_ub);
        }

        bestDual = std::max(bestDual, z_lb);

        if (bestPrimal >= z_ub) {
            improvementProcedure();
            bestPrimal = std::min(bestPrimal, z_ub);
        }

        /* calculates 2-norm of subgradient vector */
        subgradientNorm = std::inner_product(subgradient.begin(), subgradient.end(), subgradient.begin(), 0.0f);

        if (subgradientNorm < 1e-10) break;

        /* Step 5: Calculates step size */
        float stepSize = alpha * ((1.0f + beta) * bestPrimal - z_lb) / subgradientNorm;

        /* Step 6: Update lagrangean multipliers */
        for (int i = 0; i < getNumVertices(); ++i) {
            lagrangeanMultipliers[i] = std::max(0.0f, lagrangeanMultipliers[i] + stepSize * subgradient[i]);
        }

        ++iters;
    }

    printf("max lb, ub: %f -- %d\n", bestDual, bestPrimal);
}
