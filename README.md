# MTKSimulations
This directory contains multiple folders to run simulations for assessing the performance and accuracy of the new MultiTreeKnit module in [TreeKnit](https://github.com/PierreBarrat/TreeKnit.jl).

To read more about how TreeKnit works please refer to the paper: 

*TreeKnit: Inferring Ancestral Reassortment Graphs of influenza viruses
Pierre Barrat-Charlaix, Timothy G. Vaughan, Richard A. Neher bioRxiv 2021.12.20.473456; doi: https://doi.org/10.1101/2021.12.20.473456*. 


## Installation

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
Pkg.add("Dagger")
Pkg.add("Clustering")
## branch of TreeKnit with MultiTreeKnit module
Pkg.add(url="https://github.com/PierreBarrat/TreeKnit.jl", rev="MTK_no_consistency")
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
## Simulations

As I performed most simulations on a cluster I have included a bash script for launching simulations in SLURM- note that I use parallel TreeKnit which requires all computatiosn to be performed on one node. Simulations are managed by the Snakemake file in each simulation folder - with parameters in the `config.json` file. The general set up is: 
 - Simulate 8 segment ARGs with ARGTools with these parameters. ARGs map to true MCCs and true trees of each segment.
 - Take arbitrary $k \leq 8$ subsamples of these 8 segment trees as input for TreeKnit. Typically internal branches of input trees are additionally removed at random to create unresolved trees.
 - Assess the accuracy of TreeKnit's inference - write a summary of accuracy to a txt file
 - Plot results of simulations
 - remove unnecessary files, zip output files 

### 1. [Accuracy of MCCs](MCCAccuracySimulations/)

Assessment of MCC accuracy with VI and rand index. 

### 2. [Accuracy of Resolved Polytomies](PolytomySimulations/)

Using RF - distance improvement of tree and percentage of missing splits that are corrected or incorrectly added back to the tree.

### 3. [Accuracy of Shared Branches](AccuracySharedBranches/)

Using percentage of branches that are correctly and incorrectly inferred to be shared amongst tree pairs.

### 4. [Accuracy of Divergence Time Estimations](DivergenceTimeEstimations/)

Improvement of divergence time inference with [TreeTime](https://github.com/neherlab/treetime) when TreeKnit is used to infer shared branches and this information is used in ancestral sequence reconstruction and branch length/ divergence time inference. 

### 5. [Consistency of MCCs](Consistency/)

Make sure that the Maximally Compatible Clades found by TreeKnit fulfill all necessary conditions, e.g. they are topologically compatible and are consistent.