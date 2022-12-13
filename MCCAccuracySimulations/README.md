## Pipeline Steps

The MCC Accuracy Simulation pipeline assesses the improvement in MCC inference when MultiTreeKnit is used, as opposed to when standard TreeKnit is used with only two tree segments. Influenza is comprised of 8 segments, MCCs (maximally compatible clades) are clades that are shared between any pair of segment trees. When MultiTreeKnit is applied to $K$ segment trees it calculates the MCCs of all $\frac{K(K-1)}{2}$ tree pairs. As reassortments between tree pairs are not independent and additionally input trees are typically not fully resolved information from other trees and other MCCs can be useful in the calculation of individual MCCs, see the [MultiTreeKnit docs](https://github.com/PierreBarrat/TreeKnit.jl/blob/MTK_clean/docs/src/multitreeknit.md) for more information on how this is done in MultiTreeKnit. 

1. Use the julia package [ARGTools](https://github.com/PierreBarrat/ARGTools) to simulate recombination of eight segments with the flu/kingman coalescence model. See the [docs](https://github.com/PierreBarrat/ARGTools/tree/extended_newick_clean#simulations) for information on how these simulations are completed.  Obtain individual segmnet trees for each segment in the ARG as well as the real MCCs for each tree segment pair.
2. Randomly remove internal branches with probability $$e^{\frac{-\tau}{cN}}.$$ The parameter $N$ corresponds to the population size and $\tau$ is the branch length - branches are scaled to generation number as ARGTools uses a coalescence model. The parameter $c$ is determined in [MTKTools](https://github.com/anna-parker/MTKTools) and has been chosen to result in trees with a desired resolution rate (typically between 0.3 and 0.4 for influenza).
3. Take random subsamples of size $2 \leq K \leq 8$ trees and use MultiTreeKnit to infer the MCCs of all tree pairs in this random subset. (The $K=2$ case corresponds to standard TreeKnit). 
4. Use the [TestRecombTools](https://github.com/PierreBarrat/TestRecombTools) to assess the accuracy of the MCC inference using the VI index and the rand index. 
5. Write a summary of results to a txt file.
6. Plot results.
7. Remove unnecessary files and zip output files.

## MCC measures

#### VI index
#### rand index

## Link to consistency 
