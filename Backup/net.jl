# network.jl
using CSV, DataFrames

# --- Node & Link types (compact, Julia-style) ---
struct Node
    id::Int
    tails::Vector{Int}    # incoming node ids (unused now, kept for clarity)
    heads::Vector{Int}    # outgoing node ids
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

"""
    load_tntp_network(path::String)

Reads a TNTP link file (rough-cleaning) and returns:
- network::Dict{Int, Vector{Link}} adjacency by tail node
- nodes::Dict{Int, Node} map of Node structs (heads filled)
"""
function load_tntp_network(path::String)
    # read file lines and drop blank/comment lines (convert to String to avoid SubString issues)
    raw_lines = readlines(path)
    data_lines = filter(line -> begin
        s = strip(String(line))
        !(isempty(s) || startswith(s, "~") || startswith(s, "<") || startswith(s, "!"))
    end, raw_lines)

    cleaned_data = join(data_lines, "\n")
    io = IOBuffer(cleaned_data)

    df = CSV.read(io, DataFrame; delim='\t', ignorerepeated=true, header=true)

    # helper to find first existing column name from candidates
    function findcol(df, candidates...)
        for c in candidates
            if Symbol(c) in names(df)
                return Symbol(c)
            end
            # also allow symbol candidates provided
            if c in names(df)
                return c
            end
        end
        return nothing
    end

    # identify columns (common TNTP names)
    col_tail = findcol(df, :init_node, :tail, :from, "init_node", "tail", "from")
    col_head = findcol(df, :term_node, :head, :to, "term_node", "head", "to")
    col_cap  = findcol(df, :capacity, :cap, "capacity", "cap")
    col_len  = findcol(df, :length, :len, "length", "len")
    col_fft  = findcol(df, :free_flow_time, :fft, "free_flow_time", "fft")
    col_b    = findcol(df, :b, "b")
    col_power= findcol(df, :power, "power")
    col_speed= findcol(df, :speed, "speed")
    col_toll = findcol(df, :toll, "toll")
    col_ltype= findcol(df, :link_type, :linktype, "link_type", "linktype")

    # create containers
    network = Dict{Int, Vector{Link}}()
    nodes = Dict{Int, Node}()

    # iterate rows and build Link + adjacency
    for r in eachrow(df)
        # read values with fallbacks
        tail = Int(get(r, col_tail, missing))
        head = Int(get(r, col_head, missing))
        cap  = try Float64(get(r, col_cap, 1.0)) catch; 1.0 end
        len  = try Float64(get(r, col_len, 0.0)) catch; 0.0 end
        fft  = try Float64(get(r, col_fft, 1.0)) catch; 1.0 end
        b    = try Float64(get(r, col_b, 0.15)) catch; 0.15 end
        power= try Float64(get(r, col_power, 4.0)) catch; 4.0 end
        speed= try Float64(get(r, col_speed, 0.0)) catch; 0.0 end
        toll = try Float64(get(r, col_toll, 0.0)) catch; 0.0 end
        ltype= try Int(get(r, col_ltype, 0)) catch; 0 end

        # guard capacity
        cap_nonzero = cap > 0.0 ? cap : 1e-9
        initial_flow = 0.0
        cost = fft * (1.0 + b * ((initial_flow / cap_nonzero) ^ power))

        link = Link(tail, head, cap, len, fft, b, power, speed, toll, ltype, initial_flow, cost)

        # insert in network adjacency
        if !haskey(network, tail)
            network[tail] = Vector{Link}()
        end
        push!(network[tail], link)

        # ensure nodes entries exist and update heads/tails
        get!(nodes, tail, Node(tail, Int[], Int[]))
        get!(nodes, head, Node(head, Int[], Int[]))
        push!(nodes[tail].heads, head)
        push!(nodes[head].tails, tail)
    end

    # ensure nodes with no outgoing links still exist (from heads)
    for (nid, _) in nodes
        if !haskey(network, nid)
            network[nid] = Vector{Link}()   # empty outgoing list
        end
    end

    return network, nodes
end

# --- simple OD loader (CSV with columns origin,destination,trips) ---
using CSV, DataFrames

function load_od_matrix(path::String)
    # Read CSV
    df = CSV.read(path, DataFrame)

    # Convert column names to lowercase symbols
    names_lower = Symbol.(lowercase.(String.(names(df))))
    rename!(df, names(df) .=> names_lower)

    # Try to detect correct column names
    possible_names = Dict(
        :origin => first(filter(n -> occursin("origin", String(n)), names_lower)),
        :destination => first(filter(n -> occursin("dest", String(n)), names_lower)),
        :trips => first(filter(n -> occursin("trip", String(n)), names_lower))
    )

    # Rename to uniform names
    rename!(df, Dict(
        possible_names[:origin] => :origin,
        possible_names[:destination] => :destination,
        possible_names[:trips] => :trips
    ))

    # Build OD dictionary
    od_dict = Dict{Tuple{Int, Int}, Float64}()
    for row in eachrow(df)
        od_dict[(Int(row.origin), Int(row.destination))] = Float64(row.trips)
    end

    println("✅ Loaded OD pairs: ", length(od_dict))
    println("   Total trips: ", sum(values(od_dict)))

    return od_dict, df
end

