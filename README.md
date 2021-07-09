# PhyloAcc-interface

## PhyloAcc

[PhyloAcc](https://xyz111131.github.io/PhyloAcc/) is a software package that can estimate substitution rates in non-coding genomic regions of comparative data. PhyloAcc also estimates lineage-specific rate categories for each region to determine elements that may be conserved or evolving at an accelerated rate relative to the background rates. These categories can be used to test for convergent evolution among lineages of interest.

For more info, see the website linked above, and the PhyloAcc paper: 

Hu Z, Sackton TB, Edwards SV, and Liu JS. 2019. Bayesian detection of convergent rate changes of conserved noncoding elements on phylogenetic trees. https://doi.org/10.1093/molbev/msz049.

## A Python based front-end for PhyloAcc

This repository will be used to develop a user-friendly front-end for the various PhyloAcc codebases and models.

The goals of this repository include:

1. Developing and handling command-line options for ease of use.
2. Initial assessment of input data to inform model selection.
3. Batching data for easy parallelization.

