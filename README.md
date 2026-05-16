## Description
This is an R package combines the functions of `asmidx` and `CSBfilter` into a single package, which will make it easier to maintain.

## Citation
If you use this package please cite the original article where we first described asmidx and CSBfilter.

Morinaga, Gen, Darío Balcazar, Athanase Badolo, et al. “From Macro to Micro: De Novo Genomes of Aedes Mosquitoes Enable Comparative Genomics Among Close and Distant Relatives.” Genome Biology and Evolution 17, no. 8 (2025): evaf142. https://doi.org/10.1093/gbe/evaf142.

## Installation
You'll need a way to install R packages from github. One way to do this is:

```
if (!requireNamespace("devtools", quietly = TRUE))
library(devtools)
devtools::install_github('genmor/asmtools')
```
## Purpose
Many programs have been developed that allow the user to assemble genomes from sequenced reads. Alongside these assemblers, a constellation of programs have been developed to assess the quality of these genomes in terms of contiguity, gene content completeness, and error detection.

This R package bundles a handful of functions written to load and compare outputs from programs such as samtools, BUSCO/compleasm, BBMap, inspector, and Juicer. Using outputs from these programs, `asmtools` allows the users to do the following in an interactive R session:
1. load various outputs into the R environment (e.g., `read.busco`, `read.jbat.review`)
2. compare user-defined genome assembly metrics for any number of genomes (`asmidx`)
3. generate filter lists of contigs or scaffolds based on sequence length and BUSCO gene presence (`busco.contig.filter`, `busco.scaffold.filter`)
4. determine if there are any BUSCOs missing from the primary assembly can be found in the secondary/haplitg assembly

##
