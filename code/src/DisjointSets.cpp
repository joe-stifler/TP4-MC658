/**
 * File: DisjointSets.cpp
 *
 * Discipline: MC658
 * PED: Natanael Ramos
 * Professor: Cid C. de Souza
 * Data of creation: June 8, 2019
 * Author (RA 176665): Jose Ribeiro Neto <j176665@dac.unicamp.br>
 * Author (RA 171119): Felipe Lopes De Mello <f171119@dac.unicamp.br>
 *
 **/

#include <DisjointSets.h>

DisjointSets::DisjointSets() {}

DisjointSets::DisjointSets(int _n) {
    initialize(_n);
}

void DisjointSets::initialize(int _n) {
    rank.resize(_n);
    parent.resize(_n);
}

void DisjointSets::clean() {
    for (int i = 0; i < (int) parent.size(); ++i) {
        rank[i] = 0;
        parent[i] = i;
    }
}

int DisjointSets::find(int _u) {
    if (_u != parent[_u])
        parent[_u] = find(parent[_u]);

    return parent[_u];
}

void DisjointSets::unionSets(int _u, int _v) {
    _u = find(_u);
    _v = find(_v);

    if (rank[_u] > rank[_v]) parent[_v] = _u;
    else parent[_u] = _v;

    if (rank[_u] == rank[_v]) ++rank[_v];
}