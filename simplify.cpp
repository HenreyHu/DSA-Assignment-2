#include "structures.h"

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