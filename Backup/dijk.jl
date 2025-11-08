# dijk_aon.jl
using DataStructures   # PriorityQueue
using LinearAlgebra

# Dijkstra (priority queue) using link.cost on adjacency network Dict{Int,Vector{Link}}
"""
    dijkstra(network, source)

Returns (dist, prev) where
- dist::Dict{Int,Float64}: shortest distance from source
- prev::Dict{Int,Int}: predecessor node (simple tree)
"""
function dijkstra(network::Dict{Int, Vector{Link}}, source::Int)
    dist = Dict{Int, Float64}()
    prev = Dict{Int, Int}()
    for n in keys(network); dist[n] = Inf; end
    dist[source] = 0.0

    pq = PriorityQueue{Int, Float64}()     # node => priority (distance)
    pq[source] = 0.0

    while !isempty(pq)
        u = dequeue_pair!(pq)[1]
        du = dist[u]
        for link in network[u]
            v = link.head
            alt = du + link.cost
            if alt < get(dist, v, Inf)
                dist[v] = alt
                prev[v] = u
                pq[v] = alt
            end
        end
    end
    return dist, prev
end

# reconstruct path nodes (list of node ids) from prev map
function reconstruct_path(prev::Dict{Int,Int}, target::Int)
    if !(target in keys(prev)) && !haskey(prev, target)
        # maybe target==source (no predecessor) -> single node path
        return [target]
    end
    path = Int[]
    cur = target
    while cur in keys(prev)
        push!(path, cur)
        cur = prev[cur]
    end
    push!(path, cur)
    reverse!(path)
    return path
end

# All-or-nothing assignment using dijkstra per origin (group destinations)
"""
    aon_assign(network, tail_pos_map, od_dict)

Returns y::Vector{Float64} (flow on each indexed edge) assigned by AON shortest paths,
where tail_pos_map is Vector{Tuple{tail,pos_in_adj}} as produced by build_edge_index.
"""
function aon_assign(network::Dict{Int, Vector{Link}}, tail_pos_map::Vector{Tuple{Int,Int}}, od_dict::Dict{Tuple{Int,Int}, Float64})
    m = length(tail_pos_map)
    y = zeros(m)

    # build quick (u,v)->eid mapping (first occurrence)
    uv_to_eid = Dict{Tuple{Int,Int}, Int}()
    for (eid,(u,pos)) in enumerate(tail_pos_map)
        v = network[u][pos].head
        if !haskey(uv_to_eid, (u,v))
            uv_to_eid[(u,v)] = eid
        end
    end

    # group destinations by origin
    ods = Dict{Int, Vector{Tuple{Int,Float64}}}()
    for ((o,d), q) in od_dict
        if q <= 0.0; continue; end
        if !haskey(network, o)
            @warn "Origin $o not present in network — skipping"
            continue
        end
        push!(get!(ods, o, Vector{Tuple{Int,Float64}}()), (d, q))
    end

    for (o, dests) in ods
        dist, prev = dijkstra(network, o)
        for (d, q) in dests
            if !haskey(dist, d) || !isfinite(dist[d])
                @warn "Unreachable OD: $o -> $d"
                continue
            end
            nodes = reconstruct_path(prev, d)
            if length(nodes) < 2; continue; end
            for i in 1:length(nodes)-1
                u = nodes[i]; v = nodes[i+1]
                if haskey(uv_to_eid, (u,v))
                    y[uv_to_eid[(u,v)]] += q
                else
                    # fallback linear search in adjacency
                    found = false
                    for (pos, link) in enumerate(network[u])
                        if link.head == v
                            # find corresponding global index
                            for (eid, (tu,tp)) in enumerate(tail_pos_map)
                                if tu == u && tp == pos
                                    y[eid] += q
                                    found = true
                                    break
                                end
                            end
                            if found; break; end
                        end
                    end
                    if !found
                        @warn "Arc mapping fail for $u -> $v"
                    end
                end
            end
        end
    end

    return y
end

# helper: build edge index mapping and tail_pos_map
function build_edge_index(network::Dict{Int, Vector{Link}})
    tail_pos_map = Tuple{Int,Int}[]
    for u in sort(collect(keys(network)))
        for pos in eachindex(network[u])
            push!(tail_pos_map, (u, pos))
        end
    end
    return tail_pos_map
end
