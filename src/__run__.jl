printstyled("[ Info ] "; color=:green)
printstyled("Precompiling InfectionSpreading.jl module ... "; color=:blue)
include("InfectionSpreading.jl")
using .InfectionSpreading
printstyled("Finish ✓\n"; color=:green)

#localARGS = isdefined(Main, :parsed_args) ? parsed_args : error("Invalid argument")
#@show localARGS

model = models_dict[parsed_args["MODEL"]]

function check_params(::Type{T}) where T <: Union{TORUS_U_ERDOS_RENYI, TORUS_U_MATCHING, TORUS_U_RANDOM_GRAPH}
    checked = true
    if parsed_args["rows"] === nothing
        printstyled("Missing number of rows.\n"; color=:red)
        checked = false
    end
    if parsed_args["cols"] === nothing
        printstyled("Missing number of columns.\n"; color=:red)
        checked = false
    end
    if (T <: TORUS_U_ERDOS_RENYI) && (length(parsed_args["probs"]) == 0)
        printstyled("Missing probability parameter.\n"; color=:red)
        checked = false
    end
    checked
end

function check_params(::Type{T}) where T <: Union{CYCLE_U_ERDOS_RENYI, CYCLE_U_MATCHING, CYCLE_U_RANDOM_GRAPH}
    checked = true
    if length(parsed_args["nodes"]) == 0
        printstyled("Missing numbers of nodes.\n"; color=:red)
        checked = false
    end
    if (T<:CYCLE_U_ERDOS_RENYI) && (length(parsed_args["probs"]) == 0)
        printstyled("Missing probability parameter.\n"; color=:red)
        checked = false
    end
    checked
end

function check_params(::Type{TWO_CLUSTER_SBM})
    checked = true
    if length(parsed_args["nodes"]) ≠ 2
        printstyled("Missing the sizes n1,n2 of the two blocks\n"; color=:red)
        checked = false
    end
    if length(parsed_args["probs"]) ≠ 2
        printstyled("Missing the probabilities p,q for inner and outer edges probabilites\n"; color=:red)
        checked = false
    end
    checked
end

get_constructor_dict(::Type{TORUS_U_ERDOS_RENYI}) = Dict(
    :rows => parsed_args["rows"],
    :columns => parsed_args["cols"],
    :p => parsed_args["probs"][1]
    )

get_constructor_dict(::Type{TORUS_U_MATCHING}) = Dict(
    :rows => parsed_args["rows"],
    :columns => parsed_args["cols"]
    )

get_constructor_dict(::Type{TORUS_U_RANDOM_GRAPH}) = Dict(
    :rows => parsed_args["rows"],
    :columns => parsed_args["cols"],
    :α => parsed_args["alpha"]
    )

get_constructor_dict(::Type{CYCLE_U_ERDOS_RENYI}) = Dict(
    :n => parsed_args["nodes"][1],
    :p => parsed_args["probs"][1]
    )

get_constructor_dict(::Type{CYCLE_U_MATCHING}) = Dict(
    :n => parsed_args["nodes"][1]
    )

get_constructor_dict(::Type{CYCLE_U_RANDOM_GRAPH}) = Dict(
    :n => parsed_args["nodes"][1],
    :α => parsed_args["alpha"]
    )

get_constructor_dict(::Type{TWO_CLUSTER_SBM}) = Dict(
    :n1 => parsed_args["nodes"][1],
    :n2 => parsed_args["nodes"][2],
    :inner_p => parsed_args["probs"][1],
    :outer_p => parsed_args["probs"][2]
    )

using PrettyTables
function table_result(result::Dict{Float64, Float64}; precision=0.01)
    n = length(result)
    data = ones(n, 2)
    ks = sort([keys(result)...])
    for i=1:n
        data[i,1] = ks[i]
        data[i,2] = result[ks[i]]
    end

    formatter = (v, i, j) -> round(v, digits = 5);

    hl = (
        Highlighter((data, i, j) -> isapprox(data[i,1], √(2) - 1, atol=precision), crayon"fg:white bold bg:dark_gray"),
        Highlighter((data, i, j) -> isapprox(data[i,1], .5, atol=precision), crayon"fg:white bold bg:dark_gray")
    )

    pretty_table(
        data;
        header = ["p", "% infected"],
        header_crayon = crayon"yellow bold",
        highlighters = hl,
        #formatters = ft_printf("%5.4f", 1:2),
        formatters = formatter,
        tf = tf_unicode_rounded,
        limit_printing = false,
        crop = :none
    )
end

import Plots:png

check_params( model ) && begin
    printstyled("Arguments checked ✓\n"; color=:green)

    result = simulation(
        parsed_args["n-samples"],
        parsed_args["initial-infected"],
        model,
        get_constructor_dict( model );
        ε = parsed_args["granularity"],
        probs_limits = (parsed_args["probs-limits"][1], parsed_args["probs-limits"][2]),
        #plot_result = parsed_args["display"]
    )

    if parsed_args["display"]
        printstyled("Plotting results ... "; color=:yellow)
        P = plot_simulation([parsed_args["probs-limits"][1]:parsed_args["granularity"]:parsed_args["probs-limits"][2]...], result)
        printstyled("Done ✓\n"; color=:green)
    end
    
    parsed_args["show-table"] && table_result(result; precision=parsed_args["granularity"])
    printstyled("Simulation finished ✓\n"; color=:green)
    println("press ENTER to exit and eventully save figure.")
    readline()

    #=
      Purtroppo se salvo il grafico con i metodi png() o savefig()
      essi chiuderanno la finestra che mostra il grafico.
      Aspetto che aggiustino questa cosa.
      =#
    parsed_args["save-plot"] && png(P, parsed_args["output"]) 
end