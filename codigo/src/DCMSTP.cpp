/**
 * File: DCMSTP.cpp
 *
 * Discipline: MC658
 * PED: Natanael Ramos
 * Professor: Cid C. de Souza
 * Data of creation: June 8, 2019
 * Author (RA 176665): Jose Ribeiro Neto <j176665@dac.unicamp.br>
 * Author (RA 171119): Felipe Lopes De Mello <f171119@dac.unicamp.br>
 *
 **/

#include <DCMSTP.h>

DCMSTP::DCMSTP(int _n, int _limitTime, clock_t _initialTime) {
    degrees.resize(_n);

    initialTime = _initialTime;
    limitTime = (float) _limitTime;
}

void DCMSTP::addEdge(int _u, int _v, int _w) {
    edges.emplace_back(_u, _v, _w);
}

void DCMSTP::setVertexMaxDegree(int _u, int _d) {
    if (_u < degrees.size()) {
        degrees[_u] = _d;
    }
}

void DCMSTP::removeAllEdges() {
    edges.clear();
}

int DCMSTP::getNumEdges() {
    return edges.size();
}

int DCMSTP::getNumVertices() {
    return degrees.size();
}

void DCMSTP::solve() {
    printf("Solver Error: Virtual function call. Not implemented yet.");
}