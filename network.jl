using CSV, DataFrames

struct Node
    id::Int
    tails::Vector{Int}
    heads::Vector{Int}
end

mutable struct Link
    tail::Int
    head::Int
    capacity::Float64
    length::Float64
    free_flow_time::Float64
    b::Float64
    power::Float64
    speed::Float64
    toll::Float64
    link_type::Int
    flow::Float64
    cost::Float64
end

function load_tntp_network(path::String)
    lines = readlines(path)
    clean_lines = filter(l -> begin
        s = strip(l)
        !isempty(s) && all(!startswith(s, c) for c in ["~", "<", "!"])
    end, lines)

    df = CSV.read(IOBuffer(join(clean_lines, "\n")), DataFrame;
                  delim='\t', ignorerepeated=true, header=true)

    for c in names(df)
        df[!, c] = strip.(replace.(string.(df[!, c]), ';' => ""))
    end

    tofloat(x, d=0.0) = tryparse(Float64, string(x)) |> y -> y === nothing ? d : y
    toint(x, d=0) = tryparse(Int, string(x)) |> y -> y === nothing ? d : y

    df.init_node = toint.(df.init_node)
    df.term_node = toint.(df.term_node)
    df.capacity = tofloat.(df.capacity, 1.0)
    df.length = tofloat.(df.length, 0.0)
    df.free_flow_time = tofloat.(df.free_flow_time, 1.0)
    df.b = tofloat.(df.b, 0.15)
    df.power = tofloat.(df.power, 4.0)
    df.speed = tofloat.(df.speed, 0.0)
    df.toll = tofloat.(df.toll, 0.0)
    df.link_type = toint.(df.link_type, 0)

    network = Dict{Int, Vector{Link}}()
    nodes = Dict{Int, Node}()

    for r in eachrow(df)
        u, v = r.init_node, r.term_node
        cap = r.capacity > 0 ? r.capacity : 1e-9
        cost = r.free_flow_time * (1 + r.b * (0 / cap)^r.power)
        link = Link(u, v, r.capacity, r.length, r.free_flow_time,
                    r.b, r.power, r.speed, r.toll, r.link_type, 0.0, cost)
        push!(get!(network, u, Link[]), link)
        get!(nodes, u, Node(u, Int[], Int[]))
        get!(nodes, v, Node(v, Int[], Int[]))
        push!(nodes[u].heads, v)
        push!(nodes[v].tails, u)
    end

    for n in keys(nodes)
        get!(network, n, Link[])
    end

    println("Network loaded: $(length(keys(nodes))) nodes, $(sum(length(v) for v in values(network))) links.")
    return network, nodes
end

using CSV, DataFrames

function load_od_matrix(path::String)
    df = CSV.read(path, DataFrame)

    # Normalize column names
    names_lower = Symbol.(lowercase.(String.(names(df))))
    rename!(df, names(df) .=> names_lower)

    # Detect correct column names
    col_origin = first(filter(n -> occursin("origin", String(n)), names_lower))
    col_dest   = first(filter(n -> occursin("dest", String(n)), names_lower))
    col_trips  = first(filter(n -> occursin("trip", String(n)), names_lower))

    # Rename uniformly
    rename!(df, Dict(col_origin => :origin, col_dest => :destination, col_trips => :trips))

    # Build OD dictionary
    od_dict = Dict{Tuple{Int, Int}, Float64}()
    for row in eachrow(df)
        od_dict[(Int(row.origin), Int(row.destination))] = Float64(row.trips)
    end

    println("Loaded OD pairs: ", length(od_dict))
    println("Total trips: ", round(sum(values(od_dict)), digits=2))

    return od_dict, df
end


