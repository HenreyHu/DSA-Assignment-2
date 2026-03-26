APSC Polygon Simplification
CSD2183 Data Structures - Project 2
Based on Kronenfeld et al. (2020)

================================================================================
BUILDING
================================================================================

    make

Requires: C++17 compatible compiler (g++ or clang++)

================================================================================
USAGE
================================================================================

    ./simplify <input.csv> <target_vertices>

Input format: CSV with columns ring_id,vertex_id,x,y
- Ring 0 is the exterior boundary
- Rings 1+ are interior holes
- Vertices listed in order around each ring

Output: CSV to stdout with simplified polygon vertices and area statistics

================================================================================
ALGORITHM OVERVIEW
================================================================================

The APSC algorithm iteratively simplifies polygons while:
1. Preserving area exactly per ring
2. Maintaining topology (no self-intersections or ring crossings)
3. Minimizing areal displacement (greedy selection)

Key Data Structures:
- Doubly-linked circular list: Each ring stored for O(1) neighbor access
- Priority queue (min-heap): Candidates ordered by displacement
- Spatial grid index: For efficient O(1) expected-time intersection queries

Collapse Operation:
For 4 consecutive vertices A -> B -> C -> D:
1. Compute placement line E parallel to AD (Equation 1b from paper)
2. Find optimal E at intersection with AB or CD
3. Replace B,C with E: A -> E -> D
4. This preserves signed area exactly

Complexity:
- Time: O(n log n) average case
- Space: O(n) for vertex storage and spatial index

================================================================================
FILES
================================================================================

simplify.cpp  - Main implementation (~560 lines)
Makefile      - Build configuration
README.txt    - This file

================================================================================
DEPENDENCIES
================================================================================

C++17 standard library only (no external dependencies)
Tested on Linux/Unix with g++ 11+

================================================================================
REFERENCES
================================================================================

Kronenfeld, B.J., Stanislawski, L.V., Buttenfield, B.P., & Brockmeyer, T. (2020).
Simplification of polylines by segment collapse: minimizing areal displacement
while preserving area. International Journal of Cartography, 6(1), 22-46.
