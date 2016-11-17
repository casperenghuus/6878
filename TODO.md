# Tasks
## A
- Swap labels (Casper)
- Bi-partite clustering (Casper)
- Multi-layer clustering (JC)

## B
- Higher order consensus

## C
- Vizualization
- Parallelize some of the Fortunato code

# Main goals
## Re-run their algorithm
Stuck because it’s super slow, might reduce network size more

## Implement more clustering algorithms
### Multi-layer clustering
- Feed individual layers to julia code
- Runs with Julia 0.4.0
- Tensor input file similar

### Bi-partite clustering
- Aggregate into bi-partite graph
- Run scikit co-clustering

### Fancier consensus algorithm
- Mess with the cpp code

## Evaluation code
- Functional enrichment, overlap with known functional clusters (can we use GAnet for something
- [modularity](https://en.wikipedia.org/wiki/Modularity_(networks)), [conductance](https://en.wikipedia.org/wiki/Conductance_(graph))
- max(inter) - min(intra)
- Assign labels for differential expression

## Visualization

## Improve existing pipeline
