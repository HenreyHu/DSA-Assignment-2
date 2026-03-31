#include "structures.h"

// ============================================================
// Global Variables
// Contributor: CHU EN HUI VERA 2402441 chu.e@digipen.edu
// ============================================================
vector<Vertex*> g_all_vertices;
vector<Vertex*> g_ring_heads;
vector<int> g_ring_sizes;
int g_total_verts = 0;
int g_next_uid = 0;
SpatialGrid g_grid;
priority_queue<Candidate, vector<Candidate>, greater<Candidate>> g_pq;
vector<vector<pair<double, double>>> g_orig_rings;

// ============================================================
// Helper Functions
// Contributor: CHU EN HUI VERA 2402441 chu.e@digipen.edu
// ============================================================
Vertex* makeVertex(double x, double y, int ring_id) {
    Vertex* v = new Vertex(x, y, ring_id, g_next_uid++);
    g_all_vertices.push_back(v);
    return v;
}

double ringSignedArea(Vertex* head) {
    if (!head) return 0;
    double area = 0;
    Vertex* v = head;
    do {
        area += v->x * v->next->y - v->next->x * v->y;
        v = v->next;
    } while (v != head);
    return 0.5 * area;
}

double polySignedArea(const vector<pair<double, double>>& pts) {
    double area = 0;
    int n = pts.size();
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += pts[i].first * pts[j].second - pts[j].first * pts[i].second;
    }
    return 0.5 * area;
}


double cross2d(double ax, double ay, double bx, double by) {
    return ax * by - ay * bx;
}

double signedTriArea(double px, double py, double qx, double qy, double rx, double ry) {
    return 0.5 * cross2d(qx - px, qy - py, rx - px, ry - py);
}

double side(double px, double py, double qx, double qy, double rx, double ry) {
    return cross2d(qx - px, qy - py, rx - px, ry - py);
}

bool segmentsCross(double ax, double ay, double bx, double by,
    double cx, double cy, double dx, double dy) {
    double d1 = side(cx, cy, dx, dy, ax, ay);
    double d2 = side(cx, cy, dx, dy, bx, by);
    double d3 = side(ax, ay, bx, by, cx, cy);
    double d4 = side(ax, ay, bx, by, dx, dy);

    return ((d1 > EPS && d2 < -EPS) || (d1 < -EPS && d2 > EPS)) &&
        ((d3 > EPS && d4 < -EPS) || (d3 < -EPS && d4 > EPS));
}

bool segmentIntersection(double ax, double ay, double bx, double by,
    double cx, double cy, double dx, double dy,
    double& ix, double& iy) {
    double d1 = side(ax, ay, bx, by, cx, cy);
    double d2 = side(ax, ay, bx, by, dx, dy);
    if (d1 * d2 > EPS) return false;

    double d3 = side(cx, cy, dx, dy, ax, ay);
    double d4 = side(cx, cy, dx, dy, bx, by);
    if (d3 * d4 > EPS) return false;

    double denom = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx);
    if (fabs(denom) < EPS) return false;

    double t = ((cx - ax) * (dy - cy) - (cy - ay) * (dx - cx)) / denom;
    ix = ax + t * (bx - ax);
    iy = ay + t * (by - ay);
    return true;
}

bool lineLineIntersect(double x1, double y1, double x2, double y2,
    double a, double b, double c,
    double& ix, double& iy) {
    double denom = a * (x2 - x1) + b * (y2 - y1);
    if (fabs(denom) < EPS) return false;
    double t = -(a * x1 + b * y1 + c) / denom;
    ix = x1 + t * (x2 - x1);
    iy = y1 + t * (y2 - y1);
    return true;
}

int main() {
    std::cout << "HELLO WORLD\n";
    return 0;
}

// ============================================================
// Part 4: APSC Algorithm Core
// Contributor: [CHUA JIA LIANG JOEL 2403273 c.jialiang@digipen.edu ]
// ============================================================

// APSC Placement (Kronenfeld et al. 2020)
bool computePlacement(Vertex* A, Vertex* B, Vertex* C, Vertex* D,
    double& ex, double& ey, double& displacement) {
    double xA = A->x, yA = A->y;
    double xB = B->x, yB = B->y;
    double xC = C->x, yC = C->y;
    double xD = D->x, yD = D->y;

    // Line E: ax + by + c = 0 (parallel to AD)
    double a = yD - yA;
    double b = xA - xD;
    double c = -yB * xA + (yA - yC) * xB + (yB - yD) * xC + yC * xD;

    double len_AD = sqrt(a * a + b * b);
    if (len_AD < EPS) return false;

    double side_B = a * xB + b * yB;
    double side_C = a * xC + b * yC;

    if (fabs(side_B) < EPS && fabs(side_C) < EPS) {
        ex = 0.5 * (xB + xC);
        ey = 0.5 * (yB + yC);
        displacement = 0;
        return true;
    }

    bool use_AB;
    if (side_B * side_C > EPS) {
        use_AB = (fabs(side_B) >= fabs(side_C) - EPS);
    }
    else {
        use_AB = (side_B * (-c) >= -EPS);
    }

    if (use_AB) {
        if (!lineLineIntersect(xA, yA, xB, yB, a, b, c, ex, ey)) {
            ex = 0.5 * (xB + xC);
            ey = 0.5 * (yB + yC);
        }
    }
    else {
        if (!lineLineIntersect(xC, yC, xD, yD, a, b, c, ex, ey)) {
            ex = 0.5 * (xB + xC);
            ey = 0.5 * (yB + yC);
        }
    }

    // Compute displacement
    double ix, iy;
    if (use_AB) {
        if (segmentIntersection(ex, ey, xD, yD, xB, yB, xC, yC, ix, iy)) {
            double a1 = fabs(signedTriArea(ex, ey, xB, yB, ix, iy));
            double a2 = fabs(signedTriArea(ix, iy, xC, yC, xD, yD));
            displacement = a1 + a2;
        }
        else {
            double area = signedTriArea(ex, ey, xB, yB, xC, yC) +
                signedTriArea(ex, ey, xC, yC, xD, yD);
            displacement = fabs(area);
        }
    }
    else {
        if (segmentIntersection(xA, yA, ex, ey, xB, yB, xC, yC, ix, iy)) {
            double a1 = fabs(signedTriArea(xA, yA, xB, yB, ix, iy));
            double a2 = fabs(signedTriArea(ix, iy, xC, yC, ex, ey));
            displacement = a1 + a2;
        }
        else {
            double area = signedTriArea(xA, yA, xB, yB, xC, yC) +
                signedTriArea(xA, yA, xC, yC, ex, ey);
            displacement = fabs(area);
        }
    }

    return true;
}

// Candidate management
void pushCandidate(Vertex* B) {
    if (!B || !B->valid) return;
    if (g_ring_sizes[B->ring_id] < 4) return;

    Vertex* A = B->prev;
    Vertex* C = B->next;
    Vertex* D = C->next;

    if (!A->valid || !C->valid || !D->valid) return;

    double ex, ey, disp;
    if (!computePlacement(A, B, C, D, ex, ey, disp)) return;

    Candidate cand;
    cand.B = B;
    cand.ex = ex;
    cand.ey = ey;
    cand.displacement = disp;
    g_pq.push(cand);
}

bool isValidCandidate(const Candidate& c) {
    if (!c.B->valid) return false;
    if (g_ring_sizes[c.B->ring_id] < 4) return false;

    Vertex* A = c.B->prev;
    Vertex* C = c.B->next;
    Vertex* D = C->next;

    if (!A->valid || !C->valid || !D->valid) return false;

    double ex, ey, disp;
    if (!computePlacement(A, c.B, C, D, ex, ey, disp)) return false;
    if (fabs(ex - c.ex) > EPS || fabs(ey - c.ey) > EPS) return false;

    return true;
}

// Topology check
bool topologyOK(Vertex* A, Vertex* B, Vertex* C, Vertex* D, double ex, double ey) {
    auto makeKey = [](int u1, int u2) -> long long {
        if (u1 > u2) swap(u1, u2);
        return ((long long)u1 << 20) | u2;
        };

    // Edges to exclude from checking
    long long excl1 = makeKey(A->uid, B->uid);
    long long excl2 = makeKey(B->uid, C->uid);
    long long excl3 = makeKey(C->uid, D->uid);

    auto checkSeg = [&](double x1, double y1, double x2, double y2) -> bool {
        vector<long long> cells;
        g_grid.cellsForSegment(x1, y1, x2, y2, cells);

        for (long long k : cells) {
            auto it = g_grid.buckets.find(k);
            if (it == g_grid.buckets.end()) continue;

            for (auto& edge : it->second) {
                Vertex* u = edge.first;
                Vertex* v = edge.second;
                if (!u->valid || !v->valid) continue;

                long long key = makeKey(u->uid, v->uid);
                if (key == excl1 || key == excl2 || key == excl3) continue;

                if (segmentsCross(x1, y1, x2, y2, u->x, u->y, v->x, v->y))
                    return false;
            }
        }
        return true;
        };

    if (!checkSeg(A->x, A->y, ex, ey)) return false;
    if (!checkSeg(ex, ey, D->x, D->y)) return false;
    return true;
}

// Perform collapse
double performCollapse(const Candidate& c) {
    Vertex* B = c.B;
    Vertex* A = B->prev;
    Vertex* C = B->next;
    Vertex* D = C->next;
    int ring_id = B->ring_id;

    Vertex* E = makeVertex(c.ex, c.ey, ring_id);

    A->next = E;
    E->prev = A;
    E->next = D;
    D->prev = E;

    B->valid = false;
    C->valid = false;

    if (g_ring_heads[ring_id] == B || g_ring_heads[ring_id] == C)
        g_ring_heads[ring_id] = E;

    g_ring_sizes[ring_id]--;
    g_total_verts--;

    g_grid.addEdge(A, E);
    g_grid.addEdge(E, D);

    if (g_ring_sizes[ring_id] >= 4) {
        pushCandidate(A->prev);
        pushCandidate(A);
        pushCandidate(E);
        pushCandidate(D);
    }

    return c.displacement;
}

// Main simplification loop
double runSimplification(int target) {
    double total_disp = 0;

    while (g_total_verts > target && !g_pq.empty()) {
        Candidate c = g_pq.top();
        g_pq.pop();

        if (!isValidCandidate(c)) continue;

        Vertex* A = c.B->prev;
        Vertex* C = c.B->next;
        Vertex* D = C->next;

        if (!topologyOK(A, c.B, C, D, c.ex, c.ey)) continue;

        total_disp += performCollapse(c);
    }

    return total_disp;
}

