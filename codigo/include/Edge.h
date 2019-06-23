//
// Created by joe on 08/06/19.
//

#ifndef TP4_EDGE_H
#define TP4_EDGE_H

struct Edge {
    int u; /* First vertex */
    int v; /* Second vertex */
    int w; /* Edge weigth */

    Edge() { u = v = w = -1; }

    /* Constructor */
    Edge(int _u, int _v, int _w) {
        u = _u;
        v = _v;
        w = _w;
    }

    bool operator <(const Edge &other) const {
        return w < other.w;
    }
};

#endif //TP4_EDGE_H
