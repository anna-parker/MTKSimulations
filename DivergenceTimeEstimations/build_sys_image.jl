using PackageCompiler

create_sysimage(
    [:MTKTools, :TreeTools, :CSV, :ArgParse,
    :Combinatorics, :TreeAlgs, :Distributions,
    :TreeKnit],
    sysimage_path="arg_simulations.so",
    precompile_execution_file="ARGSimulations.jl" # the new line
)
