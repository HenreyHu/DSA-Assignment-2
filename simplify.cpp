#include "structures.h"

using namespace std;

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

// ============================================================
// Part 4: APSC Algorithm Core
// Contributor: [CHUA JIA LIANG JOEL 2403273 c.jialiang@digipen.edu ]
// ============================================================

bool pointOnSegment(double px, double py, double ax, double ay, double bx, double by) {
    double dx = bx - ax, dy = by - ay;
    double lenSq = dx * dx + dy * dy;
    if (lenSq < EPS * EPS) return false;

    double cross_val = dx * (py - ay) - dy * (px - ax);
    if (cross_val * cross_val > 1e-8 * lenSq) return false;

    double t = ((px - ax) * dx + (py - ay) * dy) / lenSq;
    return t > 1e-6 && t < 1.0 - 1e-6;
}

bool computePlacement(Vertex* A, Vertex* B, Vertex* C, Vertex* D,
    double& ex, double& ey, double& displacement) {
    double xA = A->x, yA = A->y;
    double xB = B->x, yB = B->y;
    double xC = C->x, yC = C->y;
    double xD = D->x, yD = D->y;

    double a = yD - yA;
    double b = xA - xD;
    double c = -yB * xA + (yA - yC) * xB + (yB - yD) * xC + yC * xD;

    double len_AD_sq = a * a + b * b;
    if (len_AD_sq < EPS * EPS) return false;

    double k = a * xA + b * yA;
    double side_B = a * xB + b * yB - k;
    double side_C = a * xC + b * yC - k;
    double side_E = -c - k;

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
    else if (side_B * side_C < -EPS) {
        use_AB = (side_B * side_E >= -EPS);
    }
    else {
        if (fabs(side_B) < EPS) {
            use_AB = false;
        }
        else {
            use_AB = true;
        }
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
    cand.uidA = A->uid;
    cand.uidB = B->uid;
    cand.uidC = C->uid;
    cand.uidD = D->uid;
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

    // Reject if any node in the sequence has been modified since the candidate was generated
    if (!A->valid || !C->valid || !D->valid) return false;
    if (A->uid != c.uidA || c.B->uid != c.uidB || C->uid != c.uidC || D->uid != c.uidD) return false;

    return true;
}

bool topologyOK(Vertex* A, Vertex* B, Vertex* C, Vertex* D, double ex, double ey) {
    auto makeKey = [](long long u1, long long u2) -> long long {
        if (u1 > u2) swap(u1, u2);
        return (u1 << 32) | u2;
        };

    long long excl1 = makeKey(A->uid, B->uid);
    long long excl2 = makeKey(B->uid, C->uid);
    long long excl3 = makeKey(C->uid, D->uid);

    int exclV1 = A->uid, exclV2 = B->uid, exclV3 = C->uid, exclV4 = D->uid;

    auto checkSeg = [&](double x1, double y1, double x2, double y2) -> bool {
        vector<long long> cells;
        double minx = min(x1, x2), maxx = max(x1, x2);
        double miny = min(y1, y2), maxy = max(y1, y2);
        g_grid.cellsForBoundingBox(minx, miny, maxx, maxy, cells);

        for (long long ck : cells) {
            auto it = g_grid.buckets.find(ck);
            if (it == g_grid.buckets.end()) continue;

            for (auto& edge : it->second) {
                Vertex* u = edge.first;
                Vertex* v = edge.second;
                if (!u->valid || !v->valid) continue;

                long long key = makeKey(u->uid, v->uid);
                if (key == excl1 || key == excl2 || key == excl3) continue;

                if (segmentsCross(x1, y1, x2, y2, u->x, u->y, v->x, v->y))
                    return false;

                if (u->uid != exclV1 && u->uid != exclV2 &&
                    u->uid != exclV3 && u->uid != exclV4) {
                    if (pointOnSegment(u->x, u->y, x1, y1, x2, y2))
                        return false;
                }
                if (v->uid != exclV1 && v->uid != exclV2 &&
                    v->uid != exclV3 && v->uid != exclV4) {
                    if (pointOnSegment(v->x, v->y, x1, y1, x2, y2))
                        return false;
                }

                if (pointOnSegment(ex, ey, u->x, u->y, v->x, v->y))
                    return false;
            }
        }
        return true;
        };

    if (!checkSeg(A->x, A->y, ex, ey)) return false;
    if (!checkSeg(ex, ey, D->x, D->y)) return false;
    return true;
}

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

// ============================================================
// Part 5: Input/Output and Main
// Contributor: [ANANTHAN ELANGOVAN 2403536 ananthan.e@digipen.edu ]
// ============================================================

void readInput(const string& filename) {
    ifstream fin(filename);
    if (!fin) {
        cerr << "Error: Cannot open " << filename << endl;
        exit(1);
    }

    string line;
    getline(fin, line);  // Skip header

    unordered_map<int, vector<tuple<int, double, double>>> ring_data;
    int max_ring_id = -1;

    while (getline(fin, line)) {
        if (line.empty()) continue;
        if (line.back() == '\r') line.pop_back();

        stringstream ss(line);
        string token;

        getline(ss, token, ','); int ring_id = stoi(token);
        getline(ss, token, ','); int vertex_id = stoi(token);
        getline(ss, token, ','); double x = stod(token);
        getline(ss, token, ','); double y = stod(token);

        ring_data[ring_id].push_back({ vertex_id, x, y });
        max_ring_id = max(max_ring_id, ring_id);
    }

    for (auto& p : ring_data)
        sort(p.second.begin(), p.second.end());

    int num_rings = max_ring_id + 1;
    g_ring_heads.resize(num_rings, nullptr);
    g_ring_sizes.resize(num_rings, 0);
    g_orig_rings.resize(num_rings);

    double min_x = 1e18, min_y = 1e18, max_x = -1e18, max_y = -1e18;

    for (int r = 0; r < num_rings; r++) {
        auto& verts = ring_data[r];
        int n = verts.size();
        if (n == 0) continue;

        vector<Vertex*> ring_verts;
        for (auto& v : verts) {
            double x = get<1>(v), y = get<2>(v);
            ring_verts.push_back(makeVertex(x, y, r));
            g_orig_rings[r].push_back({ x, y });
            min_x = min(min_x, x); min_y = min(min_y, y);
            max_x = max(max_x, x); max_y = max(max_y, y);
        }

        for (int i = 0; i < n; i++) {
            ring_verts[i]->next = ring_verts[(i + 1) % n];
            ring_verts[i]->prev = ring_verts[(i - 1 + n) % n];
        }

        g_ring_heads[r] = ring_verts[0];
        g_ring_sizes[r] = n;
        g_total_verts += n;
    }

    double range = max(max_x - min_x, max_y - min_y);
    int cells = max(1, (int)sqrt((double)g_total_verts));
    g_grid.min_x = min_x;
    g_grid.min_y = min_y;
    g_grid.cell_size = range / cells + 1.0;

    for (int r = 0; r < num_rings; r++) {
        Vertex* head = g_ring_heads[r];
        if (!head) continue;
        Vertex* v = head;
        do {
            g_grid.addEdge(v, v->next);
            v = v->next;
        } while (v != head);
    }

    for (int r = 0; r < num_rings; r++) {
        Vertex* head = g_ring_heads[r];
        if (!head || g_ring_sizes[r] < 4) continue;
        Vertex* v = head;
        do {
            pushCandidate(v);
            v = v->next;
        } while (v != head);
    }
}

void writeOutput(double total_disp) {
    double input_area = 0;
    for (size_t r = 0; r < g_orig_rings.size(); r++)
        input_area += polySignedArea(g_orig_rings[r]);

    double output_area = 0;
    for (size_t r = 0; r < g_ring_heads.size(); r++)
        output_area += ringSignedArea(g_ring_heads[r]);

    cout << "ring_id,vertex_id,x,y" << endl;

    for (size_t r = 0; r < g_ring_heads.size(); r++) {
        Vertex* head = g_ring_heads[r];
        if (!head) continue;

        int vid = 0;
        Vertex* v = head;
        do {
            cout << r << "," << vid << ",";
            auto printCoord = [](double val) {
                if (fabs(val - round(val)) < 1e-9) {
                    cout << (long long)round(val);
                }
                else {
                    cout << defaultfloat << setprecision(10) << val;
                }
                };
            printCoord(v->x);
            cout << ",";
            printCoord(v->y);
            cout << endl;
            vid++;
            v = v->next;
        } while (v != head);
    }

    cout << scientific << setprecision(6);
    cout << "Total signed area in input: " << input_area << endl;
    cout << "Total signed area in output: " << output_area << endl;
    cout << "Total areal displacement: " << total_disp << endl;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input.csv> <target_vertices>" << endl;
        return 1;
    }

    readInput(argv[1]);
    double total_disp = (g_total_verts > stoi(argv[2])) ? runSimplification(stoi(argv[2])) : 0;
    writeOutput(total_disp);

    for (Vertex* v : g_all_vertices) delete v;
    return 0;
}