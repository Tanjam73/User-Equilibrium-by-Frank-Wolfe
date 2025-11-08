# frank_wolfe.jl
using LinearAlgebra, Printf, Dates

# ----------------------------
# Line search (bisection on derivative) — compact & robust
# ----------------------------
function line_search_alpha(network, tail_pos_map, x::Vector{Float64}, y::Vector{Float64}; tol=1e-6, maxiter=80)
    m = length(x)
    g = α -> begin
        s = 0.0
        for i in 1:m
            δ = y[i] - x[i]
            f = x[i] + α * δ
            (u,pos) = tail_pos_map[i]; l = network[u][pos]
            cap = l.capacity > 0 ? l.capacity : 1e-9
            c = l.free_flow_time * (1 + l.b * ((f / cap) ^ l.power))
            s += δ * c
        end
        return s
    end

    a, b = 0.0, 1.0
    fa, fb = g(a), g(b)
    if abs(fa) < tol
        return 0.0
    end
    if fa * fb > 0
        # pick best among candidates
        cand = [0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0]
        vals = map(x -> abs(g(x)), cand)
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
    end
    return 0.5*(a+b)
end

# ----------------------------
# Main Frank–Wolfe UE solver
# ----------------------------
"""
    frank_wolfe_ue(network, od_dict; tol=1e-8, max_iter=2000, warm_start=true, verbose=false)

Frank–Wolfe UE solver (edge-indexed). Returns (x, history, tail_pos_map, conv_hist).

Important:
- Uses true relative gap: gap = dot(x - y, c(x)); rel_gap = |gap| / max(1, dot(x, c(x)))
- Warm start with AON by default (accelerates convergence).
"""
function frank_wolfe_ue(network::Dict{Int, Vector{Link}}, od_dict::Dict{Tuple{Int,Int},Float64};
                        tol=1e-8, max_iter=2000, warm_start=true, verbose=false)

    tail_pos_map = build_edge_index(network)
    m = length(tail_pos_map)

    # initial flows (warm start recommended)
    x = warm_start ? aon_assign(network, tail_pos_map, od_dict) : zeros(m)
    write_flows!(network, tail_pos_map, x)

    history = Vector{Vector{Float64}}(); push!(history, copy(x))
    conv = Float64[]
    iter = 0
    rel_gap = Inf
    t0 = now()

    while rel_gap > tol && iter < max_iter
        iter += 1

        # update costs for current x
        c_x = costs_from_flows(network, tail_pos_map, x)
        for i in 1:m
            (u,pos) = tail_pos_map[i]
            network[u][pos].cost = c_x[i]
        end

        # All-or-nothing assignment
        y = aon_assign(network, tail_pos_map, od_dict)
        sum(y) == 0 && error("AON assigned zero flow — check OD data or connectivity")

        # optimal step size via compact bisection
        α = line_search_alpha(network, tail_pos_map, x, y; tol=1e-6)
        if α == 0.0
            α = (iter == 1) ? 1.0 : (2.0 / (iter + 1))   # robust fallback
        end

        # update flows
        x_new = x .+ α .* (y .- x)
        write_flows!(network, tail_pos_map, x_new)
        push!(history, copy(x_new))

        # true relative gap (benchmark-compatible)
        cx = costs_from_flows(network, tail_pos_map, x)   # cost at x
        gap = dot(x .- y, cx)
        rel_gap = abs(gap) / max(1.0, abs(dot(x, cx)))
        push!(conv, rel_gap)

        if verbose && (iter <= 5 || iter % 50 == 0)
            total_flow = sum(x_new)
            total_cost = dot(x_new, cx)
            t = (now() - t0).value / 1000
            @printf("iter %4d | α=%.6f | rel_gap=%.3e | flow=%.3f | cost=%.3f | time=%.3f s\n",
                    iter, α, rel_gap, total_flow, total_cost, t)
        end

        x = x_new
    end

    if rel_gap <= tol
        if verbose; println("Converged in $iter iterations (rel_gap=$(round(rel_gap,digits=12)))"); end
    else
        if verbose; println("Stopped at iter $iter (rel_gap=$(round(rel_gap,digits=6)))"); end
    end

    return x, history, tail_pos_map, conv
end
