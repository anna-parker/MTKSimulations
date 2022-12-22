# MCC Conditions for ARG reconstruction


In order to be able to reconstruct Ancestral Reassortment Graphs and dendrograms of tree pairs under reassortment the MCCs (Maximally Compatible Clades) that we have inferred must fulfill certain conditions. 

First and foremost they must actually be Maximally Compatible Clades - at least up to resolution. And secondly, whenever we are view more than two trees MCC triplets must fulfill what we call "consistency" conditions - a form of the triangular inequality. For example if a reassortment event occurred between trees $tree_1$ and $tree_2$ above the clade defined by the leaves $A$ and $B$, but no reassortment event occurred between trees $tree_1$ and $tree_3$ then there must also be a reassortment event above the clade $[A, B]$ between the trees $tree_2$ and $tree_3$. 

In this section we shall look closer into these conditions and what can be done to make TreeKnit's results fulfill these conditions.

## Topological Compatibility Issues

Maximally Compatible Clades (MCCs) should fullfil their definition of being clades of two trees which share topology (also share topology upto resolution if desired). An MCC can be defined by its leaves. Given a set of leaves and a tree all leaves not in that set can be masked and the resulting subtree is the maximally compatible clade. If we calculate the split list of each clade in the two trees they should be the same (or if MCCs are defined upto resolution they should be compatible with each other). When run on two trees TreeKnit returns MCCs that fulfill these conditions. However, when MCCs are inferred on multiple trees and these trees are then further resolved in other trees we need to make sure that the first sets of MCCs that we inferred still fulfill these conditions even when trees are further resolved. 

Let us look at two cases where such an issue could arise. 

1. 
<img src="Figures/Example0_pre.png"/> 

Now let us assume that we are using `liberal` resolve and we calculate the following MCC pairs:
$$MCC_{1,2} = [[A, B, C, D, E]]$$
$$MCC_{1,3} = [[D], [A, B, C, E]]$$ 
This leads to Tree1 being resolved using Tree3 and the $(B,C)$ split being added to Tree1
But now assume that the SA does not converge properly when we calculate $MCC_{2,3}$
$$MCC_{2,3} =  [[A], [D], [B, C, E]]$$
Tree2 will not be further resolved. This means that the two clades defined by $MCC_{1,2}$ are no longer the same - however they are still compatible up to resolution. (Note that under strict resolve the Tree1 would not be resolved as the location of $D$ is unclear which means that the two clades would still be compatible up to resolution- even if the MCCs are not consistent.)
<img src="Figures/Example0_post.png"/> 

2. 
<img src="Figures/Example1_pre.png"/> Now assume the same case but with an additional Tree4. 

As before 
$$MCC_{1,2} = [[A, B, C, D, E]]$$
$$MCC_{1,3} = [[D], [A, B, C, E]]$$ 
This leads to Tree1 being resolved using Tree3 and the $(B,C)$ split being added to Tree1
$$MCC_{1,4} = [[C], [D], [A, B, E]]$$
But now assume that the SA does not converge properly when we calculate $MCC_{2,3}$
$$MCC_{2,3} =  [[A], [D], [B, C, E]]$$
$$MCC_{2,4} =  [[D], [A, B, C, E]]$$
This leads to Tree2 being resolved using Tree4 and the $(A, B)$ split being added to Tree2.
This now means that the two clades defined by $MCC_{1,2}$ are not only not the same they are no longer compatible, even up to resolution. This means that $MCC_{1,2}$ is no longer an actual MCC. Which is quite an issue.
<img src="Figures/Example1_post.png"/> 

3. 
<img src="Figures/Example2_pre.png"/> 

Let us now look at an example of topological incompatibility that could arise with `strict` resolve.

$$MCC_{1,2} = [[D], [A, B, C, E]]$$
However, as we are now using `strict` resolve and the location of `D` cannot be determined from the other tree - i.e. we could introduce the split $(B, C)$ or $(B, C, D)$ - Tree1 shall not be further resolved.
$$MCC_{1,3} = [[A, B, C, D, E]]$$ 
Now Tree1 shall be resolved using Tree3 - the split $(A, B)$ shall be introduced.
$$MCC_{2,3} =  [[A], [D], [B, C, E]]$$
This now means that the two clades defined by $MCC_{1,2}$ are not only not the same they are no longer compatible, even up to resolution. This means that $MCC_{1,2}$ is no longer an actual MCC. This is again an issue.
<img src="Figures/Example2_post.png"/> 

## Potential Ways to Solve Topological Incompatibilities

1. Improve SA

a. pre resolve all trees by introducing splits that are compatible with all tree topologies before inferring individual tree pairs. 

b. add a consistent constraint which will push the SA to choose MCCs which are consistent with previously inferred MCCs - leading to less issues with incompatible resolutions.

Note that as SA is a stochastic process this will not fix the issue but could potentially make it less likely.

2. As all these issues arise from incompatible resolution and will potentially occur with `liberal` and `strict` resolve - adding a final round of TreeKnit which nolonger resolves trees but computes MCCs using the trees as resolved in a previous round should prevent any topological incompatibilities. This will however most likely split up MCCs too much. 


## Consistency Conditions

Assume we have $k \geq 3$ trees $\text{tree}_1$, $\text{tree}_2$, ..., $\text{tree}_k$. Recombination occurs when a host has co infections and the multiple influenza strains merge and then split their RNA. This results in each of the 8 influenza segments potentially coming from a different ancestor. As co-infections are relatively unlikely (at least for diverged strains) we further assume that recombination events happen at most between two strains (if multiple were involved this can be represented as consecutive recombination events), leading to each segment in the recombined virus having one of two possible ancestors. 

Let us look now at a 3 tree subset of these $k$ trees. At each recombination event the three segments can either all have the same ancestor or one can have a different ancestor than the other two trees. Using this relationship we can define a 

### Transitive consistency relationship on all triplets of trees

This is necessary for the construction of an ARG from MCCs and trees. For the mcc pairs: $\text{MCC}_{1,2}$, $\text{MCC}_{1,3}$ and $\text{MCC}_{2,3}$, where $\text{MCC}_{1,2}$ is the set of MCCs between the trees  $\text{tree}_1$ and  $\text{tree}_2$ we have the relationship:

$$\forall m_{1,2} \in \text{MCC}_{1,2}, \forall m_{1,3} \in \text{MCC}_{1,3}: \text{   if  } m_{1,2} \cap m_{1,3} \neq \emptyset \text{   then  } \exists m_{2,3} \in \text{MCC}_{2,3} \text{   with  } m_{1,2} \cap m_{1,3} \subseteq m_{2,3}$$

Furthermore, we can distinguish between two cases:

1. $m_{1,2} \cap m_{1,3} = m_{1,2} = m_{1,3}$ In this case a recombination event has occurred between $\text{tree}_1$ and $\text{tree}_2$ as well as between $\text{tree}_1$ and $\text{tree}_3$ above the least common ancestor of the nodes in $m_{1,2} \cap m_{1,3}$. As we assume there are at most two ancestors for each recombination event $\text{tree}_2$ and $\text{tree}_3$ should not have a recombination event between them unless there was another recombination event immediately proceeding this event; therefore both $m_{1,2} \cap m_{1,3} = m_{2,3}$ and $m_{1,2} \cap m_{1,3} \subset m_{2,3}$ are possible.

2.  $m_{1,2} \cap m_{1,3} \subset m_{1,2}$ or $m_{1,2} \cap m_{1,3} \subset m_{1,3}$ In the first of these cases a recombination event has occurred between $\text{tree}_1$ and $\text{tree}_3$ but not between $\text{tree}_1$ and $\text{tree}_2$. This means that at that recombination event $\text{tree}_1$ and $\text{tree}_2$ must have the same ancestor while $\text{tree}_3$ has a different ancestor leading to the condition that $\text{tree}_2$ and $\text{tree}_3$ must have a recombination event at this same location, or in set notation: $m_{1,2} \cap m_{1,3} = m_{2,3}$. 

Furthermore, we can distinguish between two cases:

2.A. $m_{1,2} \cap m_{1,3} = m_{1,2} = m_{1,3}$ 
$\Rightarrow $ both $m_{1,2} \cap m_{1,3} = m_{2,3}$ and $m_{1,2} \cap m_{1,3} \subset m_{2,3}$ possible

2.B. $m_{1,2} \cap m_{1,3} \subset m_{1,2}$ or $m_{1,2} \cap m_{1,3} \subset m_{1,3}$ $\Rightarrow m_{1,2} \cap m_{1,3} = m_{2,3}.$ 


This can be visualized with the following examples:
1.  $$\text{MCC}_{1,2} = [[A, B, C]]$$
    $$\text{MCC}_{1,3} = [[A, B, C]]$$
    $$\Rightarrow \text{MCC}_{2,3} = [[A, B, C]]$$
2.A. $$\text{MCC}_{1,2} = [[A, B], [C]]$$
    $$\text{MCC}_{1,3}  = [[A, B], [C]]$$
    $$\Rightarrow \text{MCC}_{2,3}  = [[A, B], [C]] \text{ or } [[A, B, C]]$$
2.B. $$\text{MCC}_{1,2}  = [[A, B, C]]$$
    $$\text{MCC}_{1,3}  = [[A, B], [C]]$$
    $$\Rightarrow \text{MCC}_{2,3}  = [[A, B], [C]$$

## Pipeline Steps

1. Use the julia package [ARGTools](https://github.com/PierreBarrat/ARGTools) to simulate recombination of eight segments with the flu/kingman coalescence model. See the [docs](https://github.com/PierreBarrat/ARGTools/tree/extended_newick_clean#simulations) for information on how these simulations are completed. Obtain individual segment trees for each segment in the ARG.
2. Randomly remove internal branches with probability $$e^{\frac{-\tau}{cN}}.$$ The parameter $N$ corresponds to the population size and $\tau$ is the branch length - branches are scaled to generation number as ARGTools uses a coalescence model. The parameter $c$ is determined in [MTKTools](https://github.com/anna-parker/MTKTools) and has been chosen to result in trees with a desired resolution rate (typically between 0.3 and 0.4 for influenza).
3. Take random subsamples of size $2 \leq K \leq 8$ trees and use MultiTreeKnit to infer the MCCs of all tree pairs in this random subset, in this process TreeKnit will additionally resolve each tree using every other tree (The $K=2$ case corresponds to standard TreeKnit). 
4. Compute average number of clade pairs defined by an MCC which are not the same but are topologically compatible and the maximal number of clade pairs defined by an MCC which are not  the same and are not topologically compatible - this is an error which should be avoided. Look at cases where multiple rounds are used, no resolution occurs in the final round and where all trees are resolved as much as possible prior to pairwise TreeKnit inference to find the optimal subroutine. 
5. Also compute the total number of inconsistencies - this is calculated by taking a random sample of triplets - three inferred trees and the inferred MCCs of all their pairs. Then we calculate which percentage of branches that should be shared between $tree_i$ and $tree_j$ due to the MCCs seen between $tree_i$ and $tree_k$ and $tree_j$ and $tree_k$ but are not - this is done for every combination of the three trees and then averaged.
6. Write a summary of results to a txt file.
7. Plot results.
8. Remove unnecessary files and zip output files.