using PackageCompiler

create_sysimage(
    [:MTKTools, :TreeTools, :CSV, :ArgParse,
    :Combinatorics, :DataFrames, :StatsBase,
    :TreeKnit],
    sysimage_path="accuracy_shared_branches.so",
    precompile_execution_file="SharedBranches.jl" # the new line
)