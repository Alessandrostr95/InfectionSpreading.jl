module Models

export
# Abstract type
SirModel,
# Models
TORUS_U_ERDOS_RENYI, TORUS_U_MATCHING, TORUS_U_RANDOM_GRAPH,
CYCLE_U_ERDOS_RENYI, CYCLE_U_MATCHING, CYCLE_U_RANDOM_GRAPH,
TWO_CLUSTER_SBM,
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
mutable struct TORUS_U_RANDOM_GRAPH <: SirModel
    graph::SimpleGraph
    rows::Integer
    columns::Integer
    α::Real
    infected::Set{Integer}
    recovered::Set{Integer}
    function TORUS_U_RANDOM_GRAPH(rows::Integer, columns::Integer, α::Real)
        (α ≤ 0) &&  throw(DomainError("α must be a strictly positive value. $α given."))
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
mutable struct CYCLE_U_RANDOM_GRAPH <: SirModel
    graph::SimpleGraph
    α::Real
    infected::Set{Integer}
    recovered::Set{Integer}
    function CYCLE_U_RANDOM_GRAPH(n::Integer, α::Real)
        (n ≤ 2) && throw(DomainError("n must be greather then 2. $n given."))
        (α ≤ 0) &&  throw(DomainError("α must be a strictly positive value. $α given."))
        return new(grass_hop_rb(n, α), α, Set{Integer}(), Set{Integer}())
        #return new(trivial_cycle_SWG(n, α), α, Set{Integer}(), Set{Integer}())
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