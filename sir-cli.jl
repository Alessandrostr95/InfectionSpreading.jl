using ArgParse

all_models = [
    "TORUS_U_ERDOS_RENYI", "TORUS_U_MATCHING", "TORUS_U_RANDOM_GRAPH",
    "CYCLE_U_ERDOS_RENYI", "CYCLE_U_MATCHING", "CYCLE_U_RANDOM_GRAPH",
    "TWO_CLUSTER_SBM"
    ]

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin

        "MODEL"
            help = "graph model"
            #required = true
            range_tester = (x -> x in all_models)

        "--nodes", "-n"
            nargs = '+'
            help = "Number of nodes. If the model is TWO_CLUSTER_SBM, two values are needed, one for each cluster."
            arg_type = Int
            range_tester = (x -> x > 0)

        "--rows", "-r"
            help = "number of rows (needed only if the model is TORUS_U_<MATCHING|ERDOS_RENYI|RANDOM_GRAPH>)"
            arg_type = Int
            range_tester = (x -> x > 2)

        "--cols", "-c"
            help = "number of columns (needed only if the model is TORUS_U_<MATCHING|ERDOS_RENYI|RANDOM_GRAPH>)"
            arg_type = Int
            range_tester = (x -> x > 2)

        "-p", "--probs"
            nargs = '+'
            help = "probabilities of edges existances (needed only when erdos-renyi graphs had to be generated). If the model is TWO_CLUSTER_SBM, two values are needed, one for inner cluster prob., and one for outer."
            arg_type = Float64
            range_tester = (x -> 0 < x < 1)
            
        "--alpha", "-a"
            arg_type = Float64
            range_tester = (x -> x > 0)
            default = 1.0
        
        "--n-samples", "-S"
            help = "Number of model to sample."
            arg_type = Int
            default = 5
            range_tester = (x -> x > 0)
        
        "--probs-limits"
            nargs = 2
            help = "Limits of probabilities to simulate"
            arg_type = Float64
            range_tester = (x -> 0 < x < 1)
            default = [1.0e-6, 1 - 1.0e-6]
        
        "--initial-infected", "-I"
            nargs = '+'
            help = "The set of initial infected nodes."
            arg_type = Int
            default = [1]
            range_tester = (x -> x > 0)
        
        "--granularity", "-d"
            help = "the granularity of the probabilities to be considered"
            arg_type = Float64
            default = 0.01
            range_tester = (x -> 0 < x < 0.2)

        "-l"
            help = "describe all models"
            action = :store_true
        
        "--save-plot"
            help = "save plot in given file name"
            action = :store_true
        
        "--output", "-o"
            help = "file name where save plot and/or result data"
            default = "simulation"  

        "--display"
            help = "plots the result of the simulation"
            action = :store_true
        
        "--show-table"
            help = "prints a table containing the simulation's results"
            action = :store_true
        
        "--verbose"
            help = "verbose output"
            action = :store_true
        
        "--version", "-v"
            action = :show_version
    end

    s.version = "1.0"

    s.epilog = """
        examples:\n
        \n
        1. TORUS_U_ERDOS_RENYI\n
        \n
        \ua0\ua0\ua0$(basename(Base.source_path())) TORUS_U_ERDOS_RENYI -r 100 -c 250 -p 0.5 ...\n
        \n
        2. TORUS_U_MATCHING\n
        \n
        \ua0\ua0\ua0$(basename(Base.source_path())) TORUS_U_MATCHING -r 100 -c 250 ...\n
        \n
        3. TORUS_U_RANDOM_GRAPH\n
        \n
        \ua0\ua0\ua0$(basename(Base.source_path())) TORUS_U_RANDOM_GRAPH -r 100 -c 250 --alpha 1.4 ...\n
        \n
        4. CYCLE_U_ERDOS_RENYI\n
        \n
        \ua0\ua0\ua0$(basename(Base.source_path())) CYCLE_U_ERDOS_RENYI -n 1000 -p 0.01 ...\n
        \n
        5. CYCLE_U_MATCHING\n
        \n
        \ua0\ua0\ua0$(basename(Base.source_path())) CYCLE_U_MATCHING -n 1000 ...\n
        \n
        6. CYCLE_U_RANDOM_GRAPH\n
        \n
        \ua0\ua0\ua0$(basename(Base.source_path())) CYCLE_U_RANDOM_GRAPH -n 1000 --alpha 1.4 ...\n
        \n
        7. TWO_CLUSTER_SBM\n
        \n
        \ua0\ua0\ua0$(basename(Base.source_path())) TWO_CLUSTER_SBM -n 1000 700 -p 0.2 0.001 ...\n
        \n
        """

    return parse_args(s)
end


parsed_args = parse_commandline()
#println(ARGS)
#println(parsed_args)

parsed_args["l"] && begin
    printstyled("Avaible models :\n"; color=:green)
    for m in all_models
        println("\t$m")
    end
end

parsed_args["MODEL"] ≠ nothing && include("./src/__run__.jl")