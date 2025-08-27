# randmst-simulator


In this project I implement minimum spanning tree (MST) algorithms across multiple randomized graph models, analyzing how the expected MST weight scales as the number of vertices increases.

---

## Overview

The simulator generates several families of graphs and computes the MST weight over multiple trials:

- **Dimension 0:** Complete graphs with uniformly random edge weights \([0,1]\), with edge pruning for scalability.  
- **Dimension 1:** Hypercube-style graphs where vertices connect at distances of powers of two, with random edge weights.  
- **Dimensions 2–4:** Geometric graphs where vertices are uniformly random points in the unit square/cube/hypercube, edges weighted by Euclidean distance, and neighborhood bucketing to reduce edge counts.

MSTs are computed using Kruskal’s algorithm with an optimized union–find (disjoint set union) data structure supporting path compression and union by rank.  

---

## Usage


```bash
python randmst.py 0 numpoints numtrials dimension
```

- `0` : flag (for customization)
- `numpoints` : number of vertices in the graph
- `numtrials` : number of independent trials to average
- `dimension` :  
  - `0` → complete graph with random edge weights  
  - `1` → hypercube graph  
  - `2` → 2D geometric graph (points in unit square)  
  - `3` → 3D geometric graph (unit cube)  
  - `4` → 4D geometric graph (unit hypercube)  

Example: Run 5 trials on a 1024-vertex 3D geometric graph:

```bash
python randmst.py 0 1024 5 3
```

Output format:

```
<average_MST_weight> <numpoints> <numtrials> <dimension>
```

---

Developed as project for Harvard CS1240: Data Structures and Algorithms (Spring 2025).  
