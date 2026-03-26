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
    double ex, ey;
    double displacement;
    
    bool operator>(const Candidate& o) const {
        return displacement > o.displacement;
    }
};

#endif // STRUCTURES_H