using PackageCompiler

create_sysimage(
    [:MTKTools, :TreeTools, :CSV, :ArgParse,
    :Combinatorics, :DataFrames, :StatsBase,
    :TreeKnit],
    sysimage_path="polytomy_simulations.so",
    precompile_execution_file="../PolytomyResolutionSimulations.jl" # the new line
)