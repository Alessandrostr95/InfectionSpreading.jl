__precompile__()

module InfectionSpreading

include("Models.jl")
using .Models
include("evolution.jl")

export
# Abstract type
SirModel,
# Models
TORUS_U_ERDOS_RENYI, TORUS_U_MATCHING, TORUS_U_RANDOM_GRAPH,
CYCLE_U_ERDOS_RENYI, CYCLE_U_MATCHING, CYCLE_U_RANDOM_GRAPH,
TWO_CLUSTER_SBM,
models_dict, requirements,
# Methods
set_infected!, set_recovered!, infected, recovered, is_infected, is_recovered,
add_infected!, add_recovered!, infect!, recover!, is_susceptible,
graph, nv, ne, neighbors, neighborhood, reset_model,
reed_frost_process,

plot_process, plot_simulation,
simulation

models_dict = Dict{String, DataType}(
    "TORUS_U_ERDOS_RENYI" => TORUS_U_ERDOS_RENYI,
    "TORUS_U_MATCHING" => TORUS_U_MATCHING,
    "TORUS_U_RANDOM_GRAPH" => TORUS_U_RANDOM_GRAPH,
    "CYCLE_U_ERDOS_RENYI" => CYCLE_U_ERDOS_RENYI,
    "CYCLE_U_MATCHING" => CYCLE_U_MATCHING,
    "CYCLE_U_RANDOM_GRAPH" => CYCLE_U_RANDOM_GRAPH,
    "TWO_CLUSTER_SBM" => TWO_CLUSTER_SBM
)

requirements = Dict{String, Vector{String}}(
    "TORUS_U_ERDOS_RENYI" => ["rows", "cols", "probs"],
    "TORUS_U_MATCHING" => ["rows", "cols"],
    "TORUS_U_RANDOM_GRAPH" => ["rows", "cols", "alpha"],
    "CYCLE_U_ERDOS_RENYI" => ["nodes", "probs"],
    "CYCLE_U_MATCHING" => ["nodes"],
    "CYCLE_U_RANDOM_GRAPH" => ["nodes", "alpha"],
    "TWO_CLUSTER_SBM" => ["nodes", "probs"]
)

using Plots

"""
    Method that plot the number of infected nodes at each
    time t of a given IC process
"""
function plot_process(data::Dict{String, Any}, s::Symbol=:plot)
    # n = data["number of node"]
    series = data["infected series"]
    age = length(series)
    
    if (s == :plot)
        plot(
            1:age,
            series,
            label = "# infected",
            # ylims = (0, n),
            xlabel = "time",
            ylabel = "nodes"
        )
    elseif (s == :bar)
        bar(
            1:age,
            series,
            label = "# infected",
            # ylims = (0, n),
            xlabel = "time",
            ylabel = "nodes"
        )
    end
end

end # module
