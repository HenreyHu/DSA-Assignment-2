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

    Vertex(double x_, double y_, int rid, int id)
        : x(x_), y(y_), ring_id(rid), uid(id),
        prev(nullptr), next(nullptr), valid(true) {}
};

// ============================================================
// Candidate collapse
// ============================================================
struct Candidate {
    Vertex* B;
    int uidA, uidB, uidC, uidD; // Ensures exact sequence validation
    double ex, ey;
    double displacement;

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
    std::unordered_map<long long, std::vector<std::pair<Vertex*, Vertex*>>> buckets;

    // Safely pack 2 signed 32-bit integers into a 64-bit key
    long long cellKey(long long cx, long long cy) {
        long long ucx = cx & 0xFFFFFFFF;
        long long ucy = cy & 0xFFFFFFFF;
        return (ucx << 32) | ucy;
    }

    // Strict bounding-box collection guarantees no cells are skipped
    void cellsForBoundingBox(double minx, double miny, double maxx, double maxy, std::vector<long long>& out) {
        out.clear();
        long long min_cx = (long long)floor((minx - min_x) / cell_size);
        long long max_cx = (long long)floor((maxx - min_x) / cell_size);
        long long min_cy = (long long)floor((miny - min_y) / cell_size);
        long long max_cy = (long long)floor((maxy - min_y) / cell_size);

        for (long long cx = min_cx; cx <= max_cx; cx++) {
            for (long long cy = min_cy; cy <= max_cy; cy++) {
                out.push_back(cellKey(cx, cy));
            }
        }
    }

    void addEdge(Vertex* u, Vertex* v) {
        std::vector<long long> cells;
        double minx = std::min(u->x, v->x);
        double maxx = std::max(u->x, v->x);
        double miny = std::min(u->y, v->y);
        double maxy = std::max(u->y, v->y);
        cellsForBoundingBox(minx, miny, maxx, maxy, cells);

        for (long long k : cells) {
            buckets[k].push_back({ u, v });
        }
    }
};

#endif // STRUCTURES_H