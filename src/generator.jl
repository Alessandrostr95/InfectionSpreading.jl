using LightGraphs
using Distributions:Geometric
using Random:randperm

"""
    Returns a geometric torus graph.
    The parameters rows and columns must be at least 3.
"""
function torus_graph(rows::Integer, columns::Integer)::SimpleGraph
    (rows < 3 || columns < 3) && throw(DomainError("rows and columns must be at least 3"))
    return LightGraphs.grid([columns, rows], periodic=true)
end

"""
    Return a random perfect matching
"""
function matching_graph(n::Integer)
    (n%2 == 1) && throw(DomainError("$n is not even"))
    #nodes = [1:100...]
    #left, right = nodes[1:convert(Integer, n/2)], nodes[convert(Integer, n/2)+1:n]
    #nodes = enumerate([1:n...][convert(Integer, n/2) + 1: n])
    nodes = randperm(n)
    left, right = nodes[1:convert(Integer, n/2)], nodes[convert(Integer, n/2)+1:n]
    edgs = map(x -> Edge(x[1], x[2]), zip(left, right))
    return SimpleGraphFromIterator(edgs)
end

"""
    Given a graph g, this method joins it to a random perfect matching
"""
function matching_graph!(g::SimpleGraph)
    n = nv(g)
    (n%2 == 1) && throw(DomainError("The graph has an odd number of vertex"))
    nodes = randperm(n)
    left, right = nodes[1:convert(Integer, n/2)], nodes[convert(Integer, n/2)+1:n]
    edgs = map(x -> Edge(x[1], x[2]), zip(left, right))
    for e in edgs
        add_edge!(g, e)
    end
end

"""
    TO DO: Comment this method
"""
function index2edge(l::Integer)
    i = ceil( (√(1+8l) - 1) / 2 )
    j = ( (i - i^2) / 2 ) + l
    return convert(Integer, (i-j+1)), convert(Integer, j)
end

"""
    Funzione che mappa un indice k nell'intervallo [0, n(n-1)/2 - 1]
    nelle coordinate (i,j) della matrice triangolare superiore (diagonale esclusa)
"""
function index2edge(k::Integer, n::Integer)
    i = n - 2 - floor(√(-8*k + 4*n*(n-1)-7)/2.0 - 0.5)
    j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2
    return convert(Integer, i+1), convert(Integer, j+1)
end

"""
    Implementation of grass hopping method described in
    Section 3.2 of arxiv paper 1709.03438.

    Since only non-direct graphs will be needed,
    I optimized the process by considering only the upper right
    triangular part of the adjacency matrix, and not all of it.

    Returns and Erdos-Renyi random graph G(n, p)
"""
function grass_hop(n::Integer, p::Float64)
    edge_index = - 1
    geometric_distribution = Geometric(p)
    gap = rand(geometric_distribution) + 1  # +1 perché rand ritorna una numero ≥ 0, a noi serve > 0
    G = SimpleGraph(n)
    max_index = (n*(n-1))/2
    while edge_index + gap < max_index
        edge_index += gap
        src::Integer, dst::Integer = index2edge(edge_index, n)
        add_edge!(G, src, dst)
        gap = rand(geometric_distribution) + 1
    end
    return G
end

"""
    Joins a given graph g with an Erdos-Renji random graph G(n, p),
    where n is the number of nodes of g.
"""
function grass_hop!(g::SimpleGraph, p::Float64)
    n::Integer = LightGraphs.nv(g)
    edge_index = - 1
    geometric_distribution = Geometric(p)
    gap = rand(geometric_distribution) + 1  # +1 perché rand ritorna una numero ≥ 0, a noi serve > 0
    max_index = (n*(n-1))/2
    while edge_index + gap < max_index
        edge_index += gap
        src::Integer, dst::Integer = index2edge(edge_index, n)
        add_edge!(g, src, dst)
        gap = rand(geometric_distribution) + 1
    end
end

"""
    This function computes the normalisation factor
    as specified in Definition 1.2 of Isabella. 
"""
function compute_c(alpha::Real, n::Integer)
    c = 0
    for x in 1:n / 2
        c = c + 1 / x^alpha
    end
    return 2 * c
end

"""
    Metodo che genera in maniera efficiente un random graph sopra un toro.
    E' una semplice modifica del codice scritto da Pierluigi Crescenzi
"""
function grass_hop_over_torus(g::SimpleGraph, alpha::Real)
    n = nv(g)
    c = compute_c(alpha, n)
    max_d = (n ÷ 2) * 2     # dimostrare che questo è il diametro di un TORO (grafo)
    for d in 2:max_d
        gd = Geometric(1 / (c * d^alpha))
        src = -1
        gap = rand(gd)
        while (src + gap ≤ n)
            src = src + gap
            dst = (src + d) % n
            add_edge!(g, src + 1, dst + 1)
            gap = rand(gd)
        end
    end
    return g
end

"""
    Implementazione della procedura grass_hop sopra un cycle_graph
    scritta da Pierluigi Crescenzi
"""
function grass_hop_rb(n, alpha)
    c = compute_c(alpha, n)
    g = cycle_graph(n)
    max_d = n ÷ 2
    for d in 2:max_d
        gd::Geometric{Float64} = Geometric(1 / (c * d^alpha))
        src::Int64 = -1
        gap::Int64 = rand(gd)
        while (src + gap <= n)
            src = src + gap
            dst::Int64 = (src + d) % n
            add_edge!(g, src + 1, dst + 1)
            gap = rand(gd)
        end
    end
    return g
end

"""
    Trivial method that returns a Small-World-Graph SWG(n, α)
"""
function cycle_SWG(n::Integer, α::Real)
    c = compute_c(α, n)
    g = cycle_graph(n)
    d(u,v) = min( mod(u-v, n), mod(v-u, n) )
    for i=1:(n-1), j=i:n
        if ( rand() ≤ (1/c)*(1/(d(i,j)^α)) )
            add_edge!(g, i, j)
        end
    end
    return g
end

get_row_index(index::Integer, columns::Integer)::Integer = ((index-1)÷columns) + 1
get_column_index(index::Integer, columns::Integer)::Integer = ((index-1)%columns) + 1
get_coord(index::Integer, columns::Integer) = (get_row_index(index,columns), get_column_index(index,columns))

"""
    Function that return the distance between the nodes u, v 
    in a torus graph (rows × columns), in a constant run-time.
"""
function torus_distances(u::Integer, v::Integer, rows::Integer, columns::Integer)::Integer
    ru::Integer, cu::Integer = get_coord(u, columns)
    rv::Integer, cv::Integer = get_coord(v, columns)
    return min( mod(ru-rv, rows), mod(rv-ru, rows) ) + min( mod(cu-cv, columns), mod(cv-cu, columns) ) 
end 

"""
    Trivial method that returns a Small-World-Graph SWG(n, m, α)
"""
function torus_SWG(rows::Integer, columns::Integer, α::Real)
    n::Integer = rows*columns
    c = compute_c(α, n)
    g = torus_graph(rows, columns)
    d(u,v) = torus_distances(u, v, rows, columns)
    for i=1:(n-1), j=i:n
        if ( rand() ≤ (1/c)*(1/(d(i,j)^α)) )
            add_edge!(g, i, j)
        end
    end
    return g
end

"""
    Returns a stochastic block graph with 2 cluster, with n1, n2 nodes
"""
function two_cluster_sbm(
    n1::Integer,
    n2::Integer,
    inner_p::Float64,
    outer_p::Float64
    )::SimpleGraph

    g::SimpleGraph = SimpleGraph(n1+n2)

    # UP LEFT
    for e in edges( grass_hop(n1, inner_p) )
        add_edge!(g, e.src, e.dst)
    end

    # DOWN RIGHT
    for e in edges( grass_hop(n2, inner_p) )
        add_edge!(g, e.src+n1, e.dst+n1)
    end

    # UP RIGHT
    for e in edges( grass_hop(max(n1, n2), outer_p) )
        if (e.src < n1 && e.dst < n2) 
            add_edge!(g, e.src, e.dst+n1)
        end
    end

    # DOWN LEFT
    for e in edges( grass_hop(max(n1, n2), outer_p) )
        if (e.src < n2 && e.dst < n1) 
            add_edge!(g, e.src+n1, e.dst)
        end
    end

    return g
end

###############################################
###############################################
###############################################

function avg_deg(g::SimpleGraph)::Float64
    n = nv(g)
    sum([degree(g,v) for v in vertices(g)])/n
end

###############################################
###############################################
###############################################

"""
    Function that generate a periodic k-dimensional lattice,
    where k = length(shape).
"""
function lattice(shape::Integer...)
    LightGraphs.grid(shape, periodic=true)
end

lattice(shape::Vector{Integer}) = lattice(shape...)

"""
    Function that generate and return an hypercube
    whit `dim` dimention (at least 3-dimentional)
"""
function hypercube(dim::Int64)
    dim ≥ 3 || throw(DomainError("The dimension must be at least 3."))
    LightGraphs.grid([2 for _=1:dim])
end

function SWG(g::SimpleGraph, α::Real)
    n = nv(g)
    c = compute_c(α, n)
    for u=1:(n-1)
        dists::Vector{Int64} = dijkstra_shortest_paths(g, u).dists
        for v=u:n
            rand() ≤ (1/c)*(1/(dists[v]^α)) && add_edge!(g, u, v)
        end
    end
    return g
end