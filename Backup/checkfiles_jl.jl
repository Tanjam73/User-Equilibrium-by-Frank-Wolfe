include("net.jl")
include("dijk.jl")
include("fw.jl")

println("Checking type stability for all modules...")

net = load_network("lpf_sioux.tntp")
od = load_od("SiouxFalls_trips.csv")

edges, tail_pos_map = build_edge_index(net)
x = zeros(length(tail_pos_map))
y = zeros(length(tail_pos_map))

println("\nnet.jl")
@code_warntype load_network("lpf_sioux.tntp")
@code_warntype load_od("SiouxFalls_trips.csv")

println("\ndijk.jl")
@code_warntype dijkstra(net, 1)
@code_warntype all_or_nothing(net, tail_pos_map, od[1])

println("\nfw.jl")
@code_warntype line_search_alpha(net, tail_pos_map, x, y)
@code_warntype fw(net, od[1]; tol=1e-5, n=50)

println("\nType check complete.")



# (Quick test for one function only)


# include("dijk.jl")
# @code_warntype dijkstra(net, 1)