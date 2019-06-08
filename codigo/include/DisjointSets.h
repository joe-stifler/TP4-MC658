/**
 * File: DisjointSets.h
 *
 * Discipline: MC658
 * PED: Natanael Ramos
 * Professor: Cid C. de Souza
 * Data of creation: June 8, 2019
 * Author (RA 176665): Jose Ribeiro Neto <j176665@dac.unicamp.br>
 * Author (RA 171119): Felipe Lopes De Mello <f171119@dac.unicamp.br>
 *
 **/

#ifndef TP4_DISJOINTSETS_H
#define TP4_DISJOINTSETS_H

#include <vector>

class DisjointSets {
public:
    DisjointSets();
    DisjointSets(int _n);

    int find();
    void clean();
    int find(int _u);
    void initialize(int _n);
    void unionSets(int _u, int _v);

private:
    std::vector<int> rank;
    std::vector<int> parent;
};


#endif //TP4_DISJOINTSETS_H
