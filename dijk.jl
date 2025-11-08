using DataStructures, LinearAlgebra

# ------------------------------------------------------------
# Dijkstra’s Algorithm (PriorityQueue)
# ------------------------------------------------------------
"""
    dijkstra(net, source)

Computes shortest paths using link.cost on a Dict-based network.
Returns (dist, prev) as dictionaries.
"""
function dijkstra(net::Dict{Int,Vector{Link}}, source::Int)
    dist = Dict(n => Inf for n in keys(net))
    prev = Dict{Int,Int}()
    dist[source] = 0.0
    pq = PriorityQueue{Int,Float64}(); pq[source] = 0.0

    while !isempty(pq)
        u = dequeue_pair!(pq)[1]
        du = dist[u]
        for link in get(net, u, Link[])
            v, alt = link.head, du + link.cost
            if alt < get(dist, v, Inf)
                dist[v], prev[v], pq[v] = alt, u, alt
            end
        end
    end
    return dist, prev
end

# ------------------------------------------------------------
# Reconstruct path from predecessor map
# ------------------------------------------------------------
function reconstruct_path(prev::Dict{Int,Int}, target::Int)
    path = Int[]
    cur = target
    while haskey(prev, cur)
        push!(path, cur)
        cur = prev[cur]
    end
    push!(path, cur)
    reverse!(path)
    return path
end

# ------------------------------------------------------------
# Build flattened edge index (for FW and AoN)
# ------------------------------------------------------------
"""
    build_edge_index(net)

Returns a Vector{Tuple{Int,Int}} giving (tail_node, position_in_adjacency)
for every link in the network.
"""
function build_edge_index(net::Dict{Int,Vector{Link}})
    edges = Tuple{Int,Int}[]
    for u in sort(collect(keys(net)))
        for pos in eachindex(net[u])
            push!(edges, (u, pos))
        end
    end
    return edges
end

# ------------------------------------------------------------
# All-or-Nothing (AON) Assignment
# ------------------------------------------------------------
"""
    aon_assign(net, edges, od)

Computes AON assignment flows y::Vector{Float64}
for given network, edge map, and OD demand dict.
"""
function aon_assign(net::Dict{Int,Vector{Link}},
                    edges::Vector{Tuple{Int,Int}},
                    od::Dict{Tuple{Int,Int},Float64})

    m = length(edges)
    y = zeros(m)

    # Quick lookup (u,v) → edge index
    uv_to_eid = Dict{Tuple{Int,Int},Int}()
    for (eid, (u, pos)) in enumerate(edges)
        v = net[u][pos].head
        uv_to_eid[(u, v)] = eid
    end

    # Group destinations by origin
    grouped_od = Dict{Int,Vector{Tuple{Int,Float64}}}()
    for ((o, d), q) in od
        q > 0 || continue
        haskey(net, o) || continue
        push!(get!(grouped_od, o, []), (d, q))
    end

    # Run Dijkstra once per origin
    for (o, dests) in grouped_od
        dist, prev = dijkstra(net, o)
        for (d, q) in dests
            if !haskey(dist, d) || !isfinite(dist[d])
                continue
            end
            nodes = reconstruct_path(prev, d)
            if length(nodes) < 2; continue; end

            for i in 1:(length(nodes)-1)
                u, v = nodes[i], nodes[i+1]
                if haskey(uv_to_eid, (u, v))
                    y[uv_to_eid[(u, v)]] += q
                else
                    # fallback: find manually (very rare)
                    for (eid, (tu, tp)) in enumerate(edges)
                        link = net[tu][tp]
                        if tu == u && link.head == v
                            y[eid] += q
                            break
                        end
                    end
                end
            end
        end
    end
    println("DEBUG: Total assigned flow in AON = ", sum(y))

    return y
end
