# Sioux Falls User Equilibrium Solver (Julia)

Implemented a full User Equilibrium traffic assignment for the Sioux Falls network using the **Frank–Wolfe algorithm** with a **Dijkstra-based All-or-Nothing (AON)** shortest path assignment.

---

## Code Overview
- **`network.jl`** — loads the TNTP-format network and OD matrix, builds node and link structures.  
- **`dijk_aon.jl`** — computes shortest paths using Dijkstra’s algorithm with a priority queue and performs AON assignment to distribute OD flows.  
- **`frank_wolfe.jl`** — iteratively updates link flows and costs via the Frank–Wolfe method with precise line search and convergence tracking.  
- **`run_ue.jl`** — integrates everything: loads data, runs the solver, exports link flows (`ue_results.csv`), and plots convergence.

---

## Key Differences and Design Choices

- **Dict-based network representation** instead of adjacency matrices — offers dynamic access to outgoing links and easier handling for sparse networks.
- **`build_edge_index` mapping** — maintains a consistent edge ordering for vectorized updates.
- **Dijkstra with `PriorityQueue` (DataStructures.jl)** — more efficient than sequential relaxation.
- **Modular AON assignment (`aon_assign`)** — keeps Dijkstra reusable and separates concerns.
- **Functional Frank–Wolfe framework** — explicit flow extraction, cost updates, and convergence tracking functions.
- **Robust bisection line search** — numerically stable and precise step-size computation.
- **AON initialization** — reduces iterations by starting closer to equilibrium.
- **Normalized relative gap criterion** — uses frank–wolfe gap (dot(x - y, c) / dot(x, c)).
- **Modular code layout** — Network, Dijkstra/AON, Frank–Wolfe, Run separated for testing and maintainability.
- **Convergence visualization** — exports data and produces a clean convergence plot.

These refinements keep the mathematical core identical but make the implementation more modular and easier to extend.

---

