"""
    Given a graph, a set of infected and a contagion_probability, this function
    computer a Reed-Frost process and return a dictionary with the following informations:

     - "infected series" => the number of infected over time,
     - "total infected" => the number of total infected,
     - "% infected" => percentage of total infected,
     - "survived" => the number of total survived,
     - "% survived" => percentage of total survived,
     - "age" => the time needed to end the process,
     - "contagion probability" => the given ontagion probability,
     - "number of node" => the number of nodes of the given graph,
     - "number of edges" => the number of edges of the given graph
"""
function reed_frost_process(
    graph::SirModel,
    first_infected::T, 
    contagion_probability::Float64
    ) where T <: Union{Vector, Array, Tuple}

    reset_model(graph)

    set_infected!(graph, first_infected)

    num_infected::Vector{Int64} = Vector{Int64}()
    t::Int64 = 0

    while !isempty(infected(graph))

        append!(num_infected, length(infected(graph)))
        t += 1

        new_infected::Set{Int64} = Set{Int64}()

        for u in infected(graph)
            neig::Vector{Int64} = neighbors(graph, u)
            for v in neig
                rand() ≤ contagion_probability && push!(new_infected, v)
            end
        end

        for u in infected(graph)
            recover!(graph, u)
        end

        for v in new_infected
            infect!(graph, v)
        end

    end

    append!(num_infected, length(infected(graph)))

    n = nv(graph.graph)
    I = length(recovered(graph))
    Dict(
        "infected series" => num_infected,
        "total infected" => I,
        "% infected" => I/n,
        "survived" => n - I,
        "% survived" => (n - I) / n,
        "age" => t,
        "contagion probability" => contagion_probability,
        "number of node" => n,
        "number of edges" => ne(graph.graph)
        )
end

"""
    Method which, given a model and a dictionary with the correct values,
    constructs and returns an instance of the given model
"""
function create(model_type::DataType, params::Dict)::SirModel
    model_type == TORUS_U_ERDOS_RENYI && return TORUS_U_ERDOS_RENYI(params[:rows], params[:columns], params[:p])
    model_type == TORUS_U_MATCHING && return TORUS_U_MATCHING(params[:rows], params[:columns])
    model_type == TORUS_U_RANDOM_GRAPH && return TORUS_U_RANDOM_GRAPH(params[:rows], params[:columns], params[:α])

    model_type == CYCLE_U_ERDOS_RENYI && return CYCLE_U_ERDOS_RENYI(params[:n], params[:p])
    model_type == CYCLE_U_MATCHING && return CYCLE_U_MATCHING(params[:n])
    model_type == CYCLE_U_RANDOM_GRAPH && return CYCLE_U_RANDOM_GRAPH(params[:n], params[:α])

    model_type == TWO_CLUSTER_SBM && return TWO_CLUSTER_SBM(params[:n1], params[:n2], params[:inner_p], params[:outer_p])
end

using Statistics:mean, stdm
import Plots
# using ProgressMeter # un po' buggato, appena risolvono lo inserisco

"""
    This function simulate the RF process over a set of graphs, within a range of probabilities.
    Parameters:
        - samples::Integer => number of graph to samples,
        - first_infected => the set of first infected nodes,
        - model_type::DataType => the model to sample (must be subtype of SirModel),
        - constructor_params::Dict => the parameters for constructing the model (must comply with the requirements the given model)
        - ε::Float64 => the granularity of the probabilities to be considered (default = 0.01),
        - probs_limits::Tuple => the probability limits within which to simulate (default = (1.0e-6, 1 - 1.0e-6)),
        - plot_result => flag for plotting the result graph (dafult = false)
        - log => flag to get information about the simulation status (dafult = true)
    Output: 
        infected => percentage of total infected for each contagion probabilitiy in the probabilities range
        standard_deviation => the standard deviation of the infected for each probabilitiy
        time => the average time to finish the Reed-Frost process, for each probabilitiy
"""
function simulation(
    samples::Integer,
    first_infected::T,
    model_type::DataType,
    constructor_params::Dict;
    #ε::Float64=1.0e-4,
    ε::Float64=0.01,
    probs_limits::Tuple{Float64, Float64}=(1.0e-6, 1 - 1.0e-6),
    plot_result=false,
    log=true
    ) where T <: Union{Vector, Array, Tuple}

    !(model_type <: SirModel) && error("model_type must be a subtype of SirModel. $model_type given.")

    samples_models::Vector{SirModel} = SirModel[]
    for i=1:samples
        push!(samples_models, create(model_type, constructor_params))
        log && begin
            print("\e[2K")
            print("\e[1G")
            printstyled("Generation models at $(round((i*100)/samples; digits=2))%"; color=:yellow)   
        end
    end
    log && println()

    return simulation(
        samples_models,
        first_infected;
        ε=ε,
        probs_limits=probs_limits,
        plot_result=plot_result,
        log=log
        )
end

"""
    This function simulate the RF process over a set of graphs, within a range of probabilities.
    Parameters:
        - samples_models::Vector{<:SirModel} => a set of SirModels,
        - first_infected => the set of first infected nodes,
        - ε::Float64 => the granularity of the probabilities to be considered (default = 0.01),
        - probs_limits::Tuple => the probability limits within which to simulate (default = (1.0e-6, 1 - 1.0e-6)),
        - plot_result => flag for plotting the result graph (dafult = false)
        - log => flag to get information about the simulation status (dafult = true)
    Output: 
        infected => percentage of total infected for each contagion probabilitiy in the probabilities range
        standard_deviation => the standard deviation of the infected for each probabilitiy
        time => the average time to finish the Reed-Frost process, for each probabilitiy
"""
function simulation(
    samples_models::Vector{<:SirModel},
    first_infected::T;
    ε::Float64=0.01,
    probs_limits::Tuple{Float64, Float64}=(1.0e-6, 1 - 1.0e-6),
    plot_result=false,
    log=true
    ) where T <: Union{Vector, Array, Tuple}

    infected = Dict{Float64, Float64}()
    standard_deviation = Dict{Float64, Float64}()
    time = Dict{Float64, Float64}()

    N = length(probs_limits[1]:ε:probs_limits[2]) * length(samples_models)
    i = 0
    for p=probs_limits[1]:ε:probs_limits[2]
        res = Float64[]
        all_times = Int64[]

        for g in samples_models
            reset_model(g)

            rf_process::Dict{String, Any} = reed_frost_process(g, first_infected, p)

            push!(res, rf_process["% infected"])
            push!(all_times, rf_process["age"])

            i += 1
            log && begin
                print("\e[2K")
                print("\e[1G")
                printstyled("Simultaion at $(round((i*100)/N; digits=1))%"; color=:yellow) 
            end
        end

        infected[p] = mean(res)
        standard_deviation[p] = stdm(res, infected[p])
        time[p] = mean(all_times)
    end

    log && println()

    plot_result && plot_simulation([probs_limits[1]:ε:probs_limits[2]...], infected)

    return Dict(
        "infected" => infected,
        "standard_deviation" => standard_deviation,
        "time" => time
        )
end

"""
    Method used in 'simulation' to plot the result
"""
function plot_simulation(
    x_values,
    infected::Dict{Float64, Float64}
    )::Plots.Plot
    
    P = plot(
        x_values,
        sort([values(infected)...]),
        xlabel = "\$ p \$",
        ylabel = "% infected",
        label="\$ p \$",
        show=true
        )

    #=
    thrasholds = [√2-1, .5]
    plot!([thrasholds[1]], seriestype="vline", label="")
    annotate!(thrasholds[1], 0, text("\$ \\sqrt{2} - 1 \$", :bottom, :right, pointsize=12))

    plot!([thrasholds[2]], seriestype="vline", label="")
    annotate!(thrasholds[2], 0, text("\$ 0.5 \$", :bottom, :left, pointsize=12))
    =#

    return P
end