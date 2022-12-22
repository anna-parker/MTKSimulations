using PackageCompiler

create_sysimage(
    [:MTKTools, :TreeTools, :CSV, :ArgParse,
    :Combinatorics, :DataFrames, :StatsBase,
    :TreeKnit, :TestRecombTools, :Clustering, :Dagger],
    sysimage_path="topo_compat.so",
    precompile_execution_file="Consistency_TopoCompatibility.jl" # the new line
)