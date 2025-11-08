# run_ue.jl — runner that includes the other modules and executes UE
# put this file alongside network.jl, dijk_aon.jl, frank_wolfe.jl and data files

# -- change to your project folder (Windows example) --
cd(raw"C:\Users\jambh\.vscode\juliaproject\sioux")

using CSV, DataFrames, Printf
using Pkg
Pkg.add("Plots")

include("net.jl")
include("dijk.jl")
include("fw.jl")

println("Loading network...")
network, nodes = load_tntp_network("lpf_sioux.tntp")
println("Loaded nodes: ", length(keys(nodes)), " tails with outgoing links: ", sum(length(v) for v in values(network)))

println("Loading OD matrix...")
od_dict, od_table = load_od_matrix("SiouxFalls_trips.csv")
println("OD pairs: ", length(od_dict), " total trips: ", sum(values(od_dict)))

# quick sanity checks
if sum(values(od_dict)) == 0.0
    error("OD matrix sums to zero — check SiouxFalls_trips.csv")
end

tail_pos_map = build_edge_index(network)
println("Indexed edges: ", length(tail_pos_map))

# sample connectivity check
k = first(keys(od_dict))
origin, dest = k[1], k[2]
println("Checking dijkstra connectivity from $origin to $dest ...")
dist, prev = dijkstra(network, origin)
println("dist($dest) = ", get(dist, dest, Inf))
if isfinite(get(dist, dest, Inf))
    println("sample path: ", reconstruct_path(prev, dest))
end

println("\nRunning Frank-Wolfe User Equilibrium ...")
final_flows, history, tail_pos_map, conv = frank_wolfe_ue(network, od_dict; tol=1e-4, max_iter=10000, verbose=true)

println("\n--- Summary ---")
println("Edges: ", length(tail_pos_map))
println("Total flow (sum over links): ", round(sum(final_flows), digits=3))

println("\nSample first 20 edges:")
for i in 1:min(20, length(final_flows))
    (u,pos) = tail_pos_map[i]
    l = network[u][pos]
    @printf("Edge %3d : %3d -> %3d | flow=%8.4f | cost=%7.4f\n", i, u, l.head, l.flow, l.cost)
end

# write flows CSV
rows = Vector{NamedTuple}()
for (i,(u,pos)) in enumerate(tail_pos_map)
    l = network[u][pos]
    push!(rows, (edge_id=i, tail=u, head=l.head, flow=l.flow, cost=l.cost, capacity=l.capacity))
end
df_out = DataFrame(rows)
CSV.write("sioux_equilibrium_flows.csv", df_out)
println("\nWrote sioux_equilibrium_flows.csv")

# write convergence values
CSV.write("convergence.csv", DataFrame(rel_step = conv))
println("Wrote convergence.csv")
using Plots

using Plots

# --- Convergence Plot ---
iters = 1:length(conv)

plot(
    iters,
    conv,
    xlabel = "Iteration",
    ylabel = "Relative Step Size (log scale)",
    title = "Frank–Wolfe Convergence",
    lw = 2,
    legend = false,
    yscale = :log10,
    grid = true,
)

savefig("convergence_plot.png")
println("Saved convergence_plot.png")



