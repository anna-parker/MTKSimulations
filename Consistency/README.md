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