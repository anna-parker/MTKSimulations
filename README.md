# MTKSimulations
This directory contains multiple folders to run simulations for assessing the performance and accuracy of the new MultiTreeKnit module in [TreeKnit](https://github.com/PierreBarrat/TreeKnit.jl).

To read more about how TreeKnit works please refer to the paper: 

*TreeKnit: Inferring Ancestral Reassortment Graphs of influenza viruses
Pierre Barrat-Charlaix, Timothy G. Vaughan, Richard A. Neher bioRxiv 2021.12.20.473456; doi: https://doi.org/10.1101/2021.12.20.473456*. 


### Installation

I use multiple julia packages which are not registered and additionally use a non-master branch of TreeKnit. Thus, to avoid installation issues I always install julia packages using the following commands:

```
using Pkg
Pkg.add("PackageCompiler")
Pkg.add("TreeTools")
Pkg.add("ArgParse")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Combinatorics")
Pkg.add("StatsBase")
Pkg.add("Distributions")
## branch of TreeKnit with MultiTreeKnit module
Pkg.add(url="https://github.com/PierreBarrat/TreeKnit.jl", rev="MTK_clean")
## modules needed for TestRecombTools - package with MCC metrics designed by Pierre
Pkg.add(url="https://github.com/diegozea/ROC.jl", rev="master")
## needed for sequence simulations
Pkg.add(url="https://github.com/anna-parker/TreeAlgs-Fork")
Pkg.add(url="https://github.com/PierreBarrat/TestRecombTools")
## modules needed to simulate ARGs 
Pkg.add(url="https://github.com/PierreBarrat/ARGTools")
Pkg.add(url="https://github.com/anna-parker/MTKTools", rev="main")
## package to plot results
Pkg.add("Plots")

```

I use `Snakemake` for coordinating my simulations and running `treetime` in the divergence time simulation pipeline. To reconstruct my mamba environment run the command:
```
mamba env create -f environment.yml
```
Incase an error occurs the `feat/multiTK_arg` branch of `treetime` should be used and needs to installed via the command:
```
pip install https://github.com/neherlab/treetime/archive/feat/multiTK_arg.zip
```
### Simulations

#### 1. Accuracy of MCCs

Assessment of MCC accuracy with VI and rand index. 

