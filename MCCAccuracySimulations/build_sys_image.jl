using PackageCompiler

create_sysimage(
    [:MTKTools, :TreeTools, :CSV, :ArgParse,
    :Combinatorics, :DataFrames, :StatsBase,
    :TreeKnit, :TestRecombTools, :Clustering, :Dagger],
    sysimage_path="mcc_accuracy_simulations.so",
    precompile_execution_file="MCC_accuracy_simulations.jl" # the new line
)
