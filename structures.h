/* Start Header ************************************************************************/
/*!
\file       structures.h
\author     Hu Zhiheng Henrey
\par        zhihenghenrey.hu@digipen.edu
\date       26/3/2026
\brief      Defines core data structures used in the APSC polygon
            simplification project, including the Vertex node for
            doubly-linked circular ring representation and the
            Candidate structure for area-preserving segment collapse.

Copyright (C) 2026 DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <tuple>

using namespace std;

const double EPS = 1e-12;

// ============================================================
// Vertex (doubly-linked circular list node)
// ============================================================
struct Vertex {
    double x, y;
    int ring_id;
    int uid;
    Vertex* prev;
    Vertex* next;
    bool valid;
    int version;  // Incremented when neighborhood changes; enables cheap PQ staleness check

    Vertex(double x_, double y_, int rid, int id)
        : x(x_), y(y_), ring_id(rid), uid(id),
        prev(nullptr), next(nullptr), valid(true), version(0) {}
};

// ============================================================
// Candidate collapse
// ============================================================
struct Candidate {
    Vertex* B;
    double ex, ey;
    double displacement;
    int version;  // Snapshot of B->version when this candidate was computed

    bool operator>(const Candidate& o) const {
        return displacement > o.displacement;
    }
};

// ============================================================
// Spatial Grid Index
// Contributor: CHU EN HUI VERA 2402441 chu.e@digipen.edu
// ============================================================
struct SpatialGrid {
    double cell_size;
    double min_x, min_y;
    unordered_map<long long, vector<pair<Vertex*, Vertex*>>> buckets;

    long long cellKey(double x, double y) {
        long long cx = (long long)floor((x - min_x) / cell_size);
        long long cy = (long long)floor((y - min_y) / cell_size);
        return (cx << 20) | (cy & 0xFFFFF);
    }

    void cellsForSegment(double x1, double y1, double x2, double y2, vector<long long>& out) {
        out.clear();
        out.push_back(cellKey(x1, y1));
        long long k2 = cellKey(x2, y2);
        if (k2 != out[0]) out.push_back(k2);

        double dx = x2 - x1, dy = y2 - y1;
        double len = sqrt(dx * dx + dy * dy);
        int steps = max(1, (int)(len / cell_size) + 1);
        for (int i = 1; i < steps; i++) {
            double t = (double)i / steps;
            long long k = cellKey(x1 + t * dx, y1 + t * dy);
            if (find(out.begin(), out.end(), k) == out.end())
                out.push_back(k);
        }
    }

    void addEdge(Vertex* u, Vertex* v) {
        vector<long long> cells;
        cellsForSegment(u->x, u->y, v->x, v->y, cells);
        for (long long k : cells)
            buckets[k].push_back({ u, v });
    }

    // Purge invalid edges from a single bucket (lazy cleanup).
    void cleanupBucket(long long key) {
        auto it = buckets.find(key);
        if (it == buckets.end()) return;
        auto& vec = it->second;
        size_t w = 0;
        for (size_t r = 0; r < vec.size(); ++r) {
            if (vec[r].first->valid && vec[r].second->valid)
                vec[w++] = vec[r];
        }
        vec.resize(w);
    }
};


#endif // STRUCTURES_H