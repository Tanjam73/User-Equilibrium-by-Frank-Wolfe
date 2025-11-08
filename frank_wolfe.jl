# frank_wolfe.jl
using LinearAlgebra, Printf, Dates


# Helpers: read/write/costs

function extract_flows(net::Dict{Int,Vector{Link}}, edges::Vector{Tuple{Int,Int}})
    x = zeros(length(edges))
    @inbounds for (i,(u,pos)) in enumerate(edges)
        x[i] = net[u][pos].flow
    end
    return x
end

function write_flows!(net::Dict{Int,Vector{Link}}, edges::Vector{Tuple{Int,Int}}, x::Vector{Float64})
    @inbounds for (i,(u,pos)) in enumerate(edges)
        l = net[u][pos]; f = x[i]
        cap = l.capacity > 0 ? l.capacity : 1e-9
        l.flow = f
        l.cost = l.free_flow_time * (1 + l.b * ((f / cap) ^ l.power))
        net[u][pos] = l
    end
end

function costs_from_flows(net::Dict{Int,Vector{Link}}, edges::Vector{Tuple{Int,Int}}, x::Vector{Float64})
    c = zeros(length(x))
    @inbounds for (i,(u,pos)) in enumerate(edges)
        l = net[u][pos]; cap = l.capacity > 0 ? l.capacity : 1e-9
        c[i] = l.free_flow_time * (1 + l.b * ((x[i] / cap) ^ l.power))
    end
    return c
end


# Line search (bisection on derivative)

function line_search_alpha(net::Dict{Int,Vector{Link}}, edges::Vector{Tuple{Int,Int}},
                           x::Vector{Float64}, y::Vector{Float64}; tol=1e-10, maxiter=200)
    m = length(x)
    g = α -> begin
        s = 0.0
        @inbounds for i in 1:m
            δ = y[i] - x[i]
            (u,pos) = edges[i]; l = net[u][pos]
            f = x[i] + α * δ
            cap = l.capacity > 0 ? l.capacity : 1e-9
            s += l.free_flow_time * (1 + l.b * ((f / cap) ^ l.power)) * δ
        end
        return s
    end

    a, b = 0.0, 1.0
    fa, fb = g(a), g(b)
    if abs(fa) < tol
        return 0.0
    end
    if fa * fb > 0
        cand = (0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0)
        vals = map(c -> abs(g(c)), cand)
        return cand[argmin(vals)]
    end

    for _ in 1:maxiter
        mid = 0.5*(a+b)
        fm = g(mid)
        if abs(fm) < tol
            return mid
        end
        if fa * fm <= 0
            b, fb = mid, fm
        else
            a, fa = mid, fm
        end
        if (b - a) <= tol
            return 0.5*(a+b)
        end
    end
    return 0.5*(a+b)
end

# ----------------------------
# Frank–Wolfe UE solver
# ----------------------------
"""
    

Solve User Equilibrium (Frank–Wolfe). Returns (x, conv_history).

- net: Dict{Int,Vector{Link}}
- od : Dict{(o,d)=>demand}
"""
function frank_wolfe_ue(net::Dict{Int,Vector{Link}}, od::Dict{Tuple{Int,Int},Float64};
                        tol=1e-8, max_iter=2000, warm_start=true, verbose=false)

    edges = build_edge_index(net)
    m = length(edges)

    # warm start with AON or zeros
    x = warm_start ? aon_assign(net, edges, od) : zeros(m)
    write_flows!(net, edges, x)

    conv = Float64[]
    rel_gap = Inf
    iter = 0
    t0 = now()

    while rel_gap > tol && iter < max_iter
        iter += 1

        # update link costs (for shortest paths)
        c_x = costs_from_flows(net, edges, x)
        @inbounds for i in 1:m
            (u,pos) = edges[i]
            net[u][pos].cost = c_x[i]
        end

        # AON assignment
        y = aon_assign(net, edges, od)
        sum(y) == 0 && error("AON assigned zero flow — check OD or network.")

        # line search (bisection)
        α = line_search_alpha(net, edges, x, y; tol=1e-10)
        if α == 0.0
            α = (iter == 1) ? 1.0 : (2.0 / (iter + 1))
        end

        # update flows
        x_new = x .+ α .* (y .- x)
        write_flows!(net, edges, x_new)

        # true FW relative gap
        cx = costs_from_flows(net, edges, x)   # cost at x
        gap = dot(x .- y, cx)
        rel_gap = abs(gap) / max(1.0, abs(dot(x, cx)))
        push!(conv, rel_gap)

        if verbose && (iter <= 5 || iter % 50 == 0)
            total_flow = sum(x_new)
            total_cost = dot(x_new, cx)
            t = (now() - t0).value / 1000
            @printf("iter %4d | α=%.6f | rel_gap=%.3e | flow=%.3f | cost=%.3f | t=%.3f s\n",
                    iter, α, rel_gap, total_flow, total_cost, t)
        end

        x = x_new
    end

    if rel_gap <= tol
        if verbose; println("Converged in $iter iters, rel_gap=$(round(rel_gap,digits=12))"); end
    else
        if verbose; println("Stopped at iter $iter, rel_gap=$(round(rel_gap,digits=8))"); end
    end

    return x, conv
end
