# plot_convergence.jl
using Plots

"""
    plot_convergence(conv::Vector{Float64})

Plots convergence (relative gap vs. iteration) for Frank–Wolfe.
"""
function plot_convergence(conv::Vector{Float64})
    iters = 1:length(conv)
    p = plot(
        iters, conv,
        xlabel = "Iteration",
        ylabel = "Relative Gap",
        yscale = :log10,
        lw = 2,
        legend = false,
        grid = true,
        title = "Frank–Wolfe Convergence",
        framestyle = :box
    )
    display(p)
end

# Example usage:
# include("plot_convergence.jl")
# plot_convergence(conv)
