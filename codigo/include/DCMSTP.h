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

#ifndef TP4_DCMSTP_H
#define TP4_DCMSTP_H

#include <ctime>
#include <memory>
#include <vector>
#include <cstdio>

#include <Edge.h>

#define GET_TIME(begin, end) ((end - begin) / (float) CLOCKS_PER_SEC) /* Gets the elapsed time from fim-ini */

class DCMSTP {
public:
    /*
     * Simple constructor
     *
     * Params:
     *     _n: total number of vertices in the graph
     * */
    explicit DCMSTP(int _n, int _limitTime, clock_t _initialTime);

    int getNumEdges();
    int getNumVertices();

    /* This function tries to
     * find a solution to the
     * problem */
    virtual void solve();

    /* Adds and edge to the current
     * graph*/
    void addEdge(int _u, int _v, int _w);

    /* Sets the maximum number of
     * incident edges in vertex _u*/
    void setVertexMaxDegree(int _u, int _d);

    /* This function removes all
     * current edges */
    void removeAllEdges();

protected:
    clock_t initialTime;
    float limitTime; /* Limit of execution in seconds */
    std::vector<Edge> edges; /* Current edges */
    std::vector<int> degrees; /* Current vertices degree */
};


#endif //TP4_DCMSTP_H
