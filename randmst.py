import sys
import math
import random
from math import log
import itertools

# disjoint set union data structure
class DSU:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0]*n

    def find(self, x):
        if self.parent[x] != x:
            # Path compression - recursively find root and make x point directly to it for faster future lookups
            self.parent[x] = self.find(self.parent[x])
        # if x is its own parent, then x is the root so just return it
        return self.parent[x]

    def union(self, x, y):
        rx = self.find(x)
        ry = self.find(y)
        if rx != ry:
            # union by rank, following logic from class,
            # ensuring smaller trees are merged into larger ones to keep depth minimal
            if self.rank[rx] < self.rank[ry]:
                self.parent[rx] = ry
            elif self.rank[rx] > self.rank[ry]:
                self.parent[ry] = rx
            else:
                self.parent[ry] = rx
                self.rank[rx] += 1

# Kruskal's
def kruskal(n, edges):
    """
    Runs Kruskal's MST algo on the given list of edge and n
    Each edge is a tuple: (weight, u, v)
    Returns the total MST weight
    """
    # Sort edges by weight
    edges.sort(key=lambda e: e[0])
    dsu = DSU(n)
    mst_cost  = 0.0
    edges_used = 0

    # iterate over the sorted edges, adding lightest edges that don’t form a cycle (greedy approach)
    for w, u, v in edges:
        # If u and v are in different components, unite them
        if dsu.find(u) != dsu.find(v):
            dsu.union(u, v)
            mst_cost += w
            edges_used += 1
            # break when n-1 edges are used (an MST has n-1 edges)
            if edges_used == n - 1:
                break
    return mst_cost

###############################################################################
# Dimension 0: Complete graph on n vertices with edges uniform in [0,1]. Following
# the hint from the assignment, we try a function k(n) = log(n)/n (see writeup), and only
# generate edges with weight <= k(n). This handles large n, and we ignore edges >
# than this cutoff, since MST edges are extremely unlikely to be that large
def mst_dimension0(n, trials):
    if n < 2:
        return 0.0
    total_cost = 0.0

    # A small factor times k(n)
    cutoff = 2.0 * log(n) / n

    for _ in range(trials):
        edges = []
        # Generate edges only if weight <= cutoff
        for i in range(n):
            for j in range(i+1, n):
                w = random.random()
                if w <= cutoff:
                    edges.append((w, i, j))

        # Computes MST multiple times and return avg cost
        mst_cost = kruskal(n, edges)
        total_cost += mst_cost

    return total_cost / trials


###############################################################################
# Dimension 1: "Hypercube" adjacency. For each vertex v, we connect it to
# v + 2^i (if that is within n). then add uniform [0,1] weights for those edges
def mst_dimension1(n, trials):
    if n < 2:
        return 0.0
    total_cost = 0.0

    # max power of 2 that is < n
    max_i = int(math.floor(math.log2(n))) + 1

    # Generate edges with powers of 2 distance apart
    for _ in range(trials):
        edges = []
        for i in range(max_i):
            step = (1 << i) # computes 2^i using bit shifting
            for v in range(n):
                # Compute the vertex `w` at distance `step` from `v` and ensure within bounds
                w = v + step 
                if w < n:
                    weight = random.random()
                    edges.append((weight, v, w))

        mst_cost = kruskal(n, edges)
        total_cost += mst_cost

    return total_cost / trials


###############################################################################
# dim 2, 3, 4: Geometric MST
# We pick n points in d-dimensional unit hypercube, then connect "nearby" cells
def generate_points(n, d):
    """
    generate n points uniformly in the unit d-cube
    output a list of lists: points[i] = [x_0, x_1, ..., x_d-1]
    """
    return [[random.random() for _ in range(d)] for __ in range(n)]

def build_edges_geometric(points, d):
    """
    Bucket the points in a d-dimensional grid so that we only connect each point
    to points in the same or adjacent bucket
    this yields O(n) edges in practice, without losing edges that are likely to appear in the MST
    Long edges are very unlikely to be part of the MST because shorter edges are preferred in Kruskal’s algorithm,
    Also, in a geometric graph (where points are randomly placed in space), most edges in the MST are between nearby points
    """
    n = len(points)
    if n < 2:
        return []

    # we first want to determine the number of grid cells
    # we can set the number of grid cells per side ~ n^(1/d)
    # We round it so that side^d ~ n
    side = int(round(n ** (1.0 / d)))
    side = max(side, 1)  # avoid zero

    # now we assign points to buckets. Each point is mapped to a grid cell
    cell_size = 1.0 / side
    # Map cell coordinates -> list of (point_index).
    cell_map = {}

    def get_cell(coords):
        # convert d-dimensional point to a d-dimensional cell index (tuple)
        idx = []
        for c in coords:
            cindex = int(c / cell_size)
            # Edge case if c=1.0 due to float rounding
            if cindex == side:
                cindex = side - 1
            idx.append(cindex)
        return tuple(idx)

    # place points in cell_map
    for i, p in enumerate(points):
        cidx = get_cell(p)
        cell_map.setdefault(cidx, []).append(i)

    # List out edges
    edges = []

    # Precompute neighbor offsets in d dimensions/directions: all combinations in [-1,0,1]^d
    neighbor_offsets = list(itertools.product([-1, 0, 1], repeat=d))

    # computes edges only between points in the same or adjacent cells
    for cell_idx, point_indices in cell_map.items():
        # compare each cell with itself and neighbors
        for offset in neighbor_offsets:
            neighbor_idx = tuple(cell_idx[i] + offset[i] for i in range(d))
            # we'll only process neighbor cells where neighbor_idx >= cell_idx 
            # in lex order to avoid dupes
            if neighbor_idx < cell_idx:
                continue
            if neighbor_idx in cell_map:
                neighbor_points = cell_map[neighbor_idx]
                for i in range(len(point_indices)):
                    for j in range(len(neighbor_points)):
                        # if it's the same cell, only take j > i to avoid duplicates
                        if neighbor_idx == cell_idx and j <= i:
                            continue
                        pi = point_indices[i]
                        pj = neighbor_points[j]

                        # compute Euclidean distance
                        dist_sq = 0.0
                        for dd in range(d):
                            diff = points[pi][dd] - points[pj][dd]
                            dist_sq += diff * diff
                        dist = math.sqrt(dist_sq)

                        edges.append((dist, pi, pj))

    return edges

# run generate_points(n, d) to create random pts
# run build_edges_geometric(points, d) to construct the edge list
# run kruskal(n, edges) to find the MST
def mst_dimension_d(n, trials, d):
    total_cost = 0.0
    for _ in range(trials):
        pts   = generate_points(n, d)
        edges = build_edges_geometric(pts, d)
        mst_cost = kruskal(n, edges)
        total_cost += mst_cost
    return total_cost / trials


def main():
    """
    Expects command line:
      randmst 0 numpoints numtrials dimension
    """
    if len(sys.argv) != 5:
        print("Usage: randmst 0 numpoints numtrials dimension")
        return

    # flag is for formality / match problem statement
    flag       = int(sys.argv[1])
    numpoints  = int(sys.argv[2])
    numtrials  = int(sys.argv[3])
    dimension  = int(sys.argv[4])

    # use appropriate MST generator
    if dimension == 0:
        avg_mst = mst_dimension0(numpoints, numtrials)
    elif dimension == 1:
        avg_mst = mst_dimension1(numpoints, numtrials)
    else:
        # dimension in [2, 3, 4]
        avg_mst = mst_dimension_d(numpoints, numtrials, dimension)

    # Print in the requested format
    print(f"{avg_mst} {numpoints} {numtrials} {dimension}")

if __name__ == "__main__":
    main()
