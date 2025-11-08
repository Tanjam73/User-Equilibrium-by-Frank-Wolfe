# run_ue.jl
using CSV, DataFrames, Printf, Dates
include("network.jl")
include("dijk.jl")
include("frank_wolfe.jl")

println("Loading network and OD matrix...")

# Adjust the file paths to your actual local files
net_path = raw"C:\Users\jambh\.vscode\juliaproject\sioux\lpf_sioux.tntp"
od_path  = raw"C:\Users\jambh\.vscode\juliaproject\sioux\SiouxFalls_trips.csv"

# Load network and OD matrix
network, nodes = load_tntp_network(net_path)
od_dict, od_df = load_od_matrix(od_path)

println("Network and OD loaded successfully.")
println("Network nodes: ", length(keys(network)))
println("OD pairs: ", length(od_dict))
println("------------------------------------")

# Run Frank–Wolfe UE solver
println("Running Frank–Wolfe User Equilibrium solver...")
x, conv = frank_wolfe_ue(network, od_dict; tol=1e-8, max_iter=3000, warm_start=true, verbose=true)

println("Completed UE assignment.")
println("Final relative gap: ", conv[end])
println("Total iterations: ", length(conv))

# Write UE results to CSV
edges = build_edge_index(network)
rows = [(tail=u, head=network[u][pos].head,
         flow=network[u][pos].flow,
         cost=network[u][pos].cost) for (u,pos) in edges]

CSV.write("ue_results.csv", DataFrame(rows))
println("Results written to ue_results.csv")

# Plot convergence (optional)
using Plots
plot(conv, yscale=:log10, xlabel="Iteration", ylabel="Relative Gap (log scale)",
     title="Convergence Plot", legend=false)
savefig("convergence_plot.png")
println("Convergence plot saved as convergence_plot.png")
