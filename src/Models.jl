module Models

export
# Abstract type
SirModel,
# Models
NETWORK, ERDOS_RENYI,
TORUS_U_ERDOS_RENYI, TORUS_U_MATCHING, TORUS_U_POWERLAW,
CYCLE_U_ERDOS_RENYI, CYCLE_U_MATCHING, CYCLE_U_POWERLAW,
TWO_CLUSTER_SBM,
LATTICE, LATTICE_U_POWERLAW,
HYPERCUBE, HYPERCUBE_U_POWERLAW,
# Methods
set_infected!, set_recovered!, infected, recovered, is_infected, is_recovered,
add_infected!, add_recovered!, infect!, recover!, is_susceptible,
graph, nv, ne, neighbors, neighborhood, reset_model

using Base: Integer
using LightGraphs
include("generator.jl")

using LightGraphs

abstract type SirModel end

"""
    A simple network based on a given graph
"""
mutable struct NETWORK <: SirModel
    graph::SimpleGraph
    infected::Set{Integer}
    recovered::Set{Integer}
    function NETWORK(graph::SimpleGraph, infected::Set{Integer}, recovered::Set{Integer})
        return new(graph, infected, recovered)
    end
end
NETWORK(graph) = NETWORK(graph, Set{Integer}(), Set{Integer}())


"""
    An Erdos-Renyi random graph
"""
mutable struct ERDOS_RENYI <: SirModel
    graph::SimpleGraph
    p::Float64
    infected::Set{Integer}
    recovered::Set{Integer}
    function ERDOS_RENYI(n::Integer, p::Float64=0.5)
        ( p < 0 || p > 1 ) && throw(DomainError("p must be between 0 and 1. $p given."))
        return new(grass_hop(n, p), p, Set{Integer}(), Set{Integer}())
    end
end

"""
    Torus graph UNION Erdos-Renyi random graph
"""
mutable struct TORUS_U_ERDOS_RENYI <: SirModel
    graph::SimpleGraph
    rows::Integer
    columns::Integer
    p::Real
    infected::Set{Integer}
    recovered::Set{Integer}
    function TORUS_U_ERDOS_RENYI(rows::Integer, columns::Integer, p::Real=0.5)
        ( p < 0 || p > 1 ) && throw(DomainError("p must be between 0 and 1. $p given."))
        t = torus_graph(rows, columns)
        grass_hop!(t, p)
        return new(t, rows, columns, p, Set{Integer}(), Set{Integer}())
    end
end

"""
    Torus graph UNION random perfect matching
"""
mutable struct TORUS_U_MATCHING <: SirModel
    graph::SimpleGraph
    rows::Integer
    columns::Integer
    infected::Set{Integer}
    recovered::Set{Integer}
    function TORUS_U_MATCHING(rows::Integer, columns::Integer)
        t = torus_graph(rows, columns)
        matching_graph!(t)
        return new(t, rows, columns, Set{Integer}(), Set{Integer}())
    end
end

"""
    Torus graph UNION random graph where the probability
    of existence between two nodes u,v is
        C(α,n)^-1 * d(u,v)^-α
    where
    - α > 0 is a given constant value
    - d(u,v) is distance between u and v on the torus graph
"""
mutable struct TORUS_U_POWERLAW <: SirModel
    graph::SimpleGraph
    rows::Integer
    columns::Integer
    α::Real
    infected::Set{Integer}
    recovered::Set{Integer}
    function TORUS_U_POWERLAW(rows::Integer, columns::Integer, α::Real)
        (α ≤ 0) && throw(DomainError("α must be a strictly positive value. $α given."))
        #t = torus_graph(rows, columns)
        t = torus_SWG(rows, columns, α)
        #return new(grass_hop_over_torus(t, α), rows, columns, α, Set{Integer}(), Set{Integer}())
        return new(t, rows, columns, α, Set{Integer}(), Set{Integer}())
    end
end

"""
    Ring graph UNION Erdos-Renyi random graph
"""
mutable struct CYCLE_U_ERDOS_RENYI <: SirModel
    graph::SimpleGraph
    p::Real
    infected::Set{Integer}
    recovered::Set{Integer}
    function CYCLE_U_ERDOS_RENYI(n::Integer, p::Real=0.5)
        (n ≤ 2) && throw(DomainError("n must be greather then 2. $n given."))
        g = cycle_graph(n)
        grass_hop!(g, p)
        return new(g, p, Set{Integer}(), Set{Integer}())
    end
end

"""
    Ring graph UNION Erdos-Renyi random graph
"""
mutable struct CYCLE_U_MATCHING <: SirModel
    graph::SimpleGraph
    infected::Set{Integer}
    recovered::Set{Integer}
    function CYCLE_U_MATCHING(n::Integer)
        (n ≤ 2) && throw(DomainError("n must be greather then 2. $n given."))
        g = cycle_graph(n)
        matching_graph!(g)
        return new(g, Set{Integer}(), Set{Integer}())
    end
end

"""
    Ring graph UNION random graph where the probability
    of existence between two nodes u,v is
        C(α,n)^-1 * d(u,v)^-α
    where
    - α > 0 is a given constant value
    - d(u,v) is distance between u and v on the ring
"""
mutable struct CYCLE_U_POWERLAW <: SirModel
    graph::SimpleGraph
    α::Real
    infected::Set{Integer}
    recovered::Set{Integer}
    function CYCLE_U_POWERLAW(n::Integer, α::Real)
        (n ≤ 2) && throw(DomainError("n must be greather then 2. $n given."))
        (α ≤ 0) && throw(DomainError("α must be a strictly positive value. $α given."))
        return new(grass_hop_rb(n, α), α, Set{Integer}(), Set{Integer}())
        #return new(cycle_SWG(n, α), α, Set{Integer}(), Set{Integer}())
    end
end

"""
    SBM stands for Stochastic Block Model
    This model is composed by the union of two Erdos-Renyi random graph
    with probability 'inner_p', which represent two communities.
    While 'outer_p' is the probability that there are edges between the two communities.
"""
mutable struct TWO_CLUSTER_SBM <: SirModel
    n1::Integer
    n2::Integer
    inner_p::Float64
    outer_p::Float64
    graph::SimpleGraph
    infected::Set{Integer}
    recovered::Set{Integer}
    function TWO_CLUSTER_SBM(n1::Integer, n2::Integer, inner_p::Float64, outer_p::Float64)
        (n1 < 2 || n2 < 2) && throw(DomainError("The cluster must have at least 2 vertices."))
        return new(n1, n2, inner_p, outer_p, two_cluster_sbm(n1, n2, inner_p, outer_p), Set{Integer}(), Set{Integer}())
    end
end

"""
    A generic k-dimensional periodic lattice, where k = length(shape).
"""
mutable struct LATTICE <: SirModel
    graph::SimpleGraph
    shape::NTuple{N,Integer} where N
    infected::Set{Integer}
    recovered::Set{Integer}
    function LATTICE(shape::Integer...)
        return new(lattice(shape...), shape, Set{Integer}(), Set{Integer}())
    end
end

LATTICE(shape::Vector{T}) where T <: Integer = LATTICE(shape...)
LATTICE(shape::Tuple{T}) where T <: Integer = LATTICE(shape...)

mutable struct LATTICE_U_POWERLAW <: SirModel
    graph::SimpleGraph
    shape::NTuple{N,Integer} where N
    α::Real
    infected::Set{Integer}
    recovered::Set{Integer}
    function LATTICE_U_POWERLAW(α::Real, shape::Integer...)
        (α ≤ 0) && throw(DomainError("α must be a strictly positive value. $α given."))
        l = lattice(shape...)
        return new(SWG(l, α), shape, α, Set{Integer}(), Set{Integer}())
    end
end

LATTICE_U_POWERLAW(α::Real, shape::Vector{T}) where T <: Integer = LATTICE_U_POWERLAW(α, shape...)
LATTICE_U_POWERLAW(α::Real, shape::Tuple{T}) where T <: Integer = LATTICE_U_POWERLAW(α, shape...)

"""
    A N-dimension hypercube
"""
mutable struct HYPERCUBE <: SirModel
    graph::SimpleGraph
    dimension::Int64
    infected::Set{Integer}
    recovered::Set{Integer}
    function HYPERCUBE(N::Int64)
        return new(hypercube(N), N, Set{Integer}(), Set{Integer}())
    end
end

mutable struct HYPERCUBE_U_POWERLAW <: SirModel
    graph::SimpleGraph
    dimension::Int64
    α::Real
    infected::Set{Integer}
    recovered::Set{Integer}
    function HYPERCUBE_U_POWERLAW(N::Int64, α::Real)
        (α ≤ 0) && throw(DomainError("α must be a strictly positive value. $α given."))
        h = hypercube(N)
        return new(SWG(h, α), N, α, Set{Integer}(), Set{Integer}())
    end
end

## METHODS

"""
    Sets the infected nodes
"""
function set_infected!(x::SirModel, nodes::Integer...)
    x.infected = Set{Integer}(nodes)
end
set_infected!(x::SirModel, nodes::Vector{Integer}) = set_infected!(x, nodes...)
set_infected!(x::SirModel, nodes::Vector{Int64}) = set_infected!(x, nodes...)

"""
    Sets the recovered nodes
"""
function set_recovered!(x::SirModel, nodes::Integer...)
    x.recovered = Set{Integer}(nodes)
end
set_recovered!(x::SirModel, nodes::Vector{Integer}) = set_infected!(x, nodes...)
set_recovered!(x::SirModel, nodes::Vector{Int64}) = set_infected!(x, nodes...)

"""
    Returns the set of infected nodes
"""
infected(x::SirModel) = x.infected

"""
    Returns the set of recovered nodes
"""
recovered(x::SirModel) = x.recovered

"""
    Returns True if the node v is infected,
    False otherwise
"""
is_infected(x::SirModel, v::Integer)::Bool = v in infected(x)

"""
    Returns True if the node v is recovered,
    False otherwise
"""
is_recovered(x::SirModel, v::Integer)::Bool = v in recovered(x)

"""
    Adds v to the set of infected
"""
add_infected!(x::SirModel, v::Integer) = push!(x.infected, v)

"""
    Adds v to the set of recovered
"""
add_recovered!(x::SirModel, v::Integer) = push!(x.recovered, v)

"""
    Adds v to the set of infected,
    making sure it's not already recovered
"""
function infect!(x::SirModel, v::Integer)
    !(is_recovered(x, v)) && push!(x.infected, v)
end

"""
    Removes v from set of infected,
    and inserts it into set of recovered
"""
function recover!(x::SirModel, v::Integer)
    delete!(x.infected, v)
    push!(x.recovered, v)
end

"""
    Returns True if v is susceptible,
    False otherwise
"""
is_susceptible(x::SirModel, v::Integer)::Bool = !(is_infected(x, v)) && !(is_recovered(x, v))

"""
    Resets the sets of infected and recovered
"""
function reset_model(x::SirModel)
    x.infected = Set{Integer}()
    x.recovered = Set{Integer}()
    return nothing
end

"""
    Return the graph of the model
"""
graph(x::SirModel)::SimpleGraph = x.graph
#number_vertices(x::SirModel) = LightGraphs.nv(graph(x))
#ne(x::SirModel) = LightGraphs.ne(graph(x))
neighbors(x::SirModel, v::Integer) = LightGraphs.neighbors(graph(x), v)
neighborhood(x::SirModel, v::Integer, d) = LightGraphs.neighborhood(graph(x), v, d)

end # module