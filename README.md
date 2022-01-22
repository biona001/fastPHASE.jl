# fastPHASE.jl

fastPHASE.jl is a Julia wrapper of the popular [fastPHASE](https://stephenslab.uchicago.edu/software.html#fastphase) genetics program, designed for haplotype reconstruction and estimating missing genotypes from population data. 

The methodology is described in the following paper

+ Scheet, P and Stephens, M (2006). A fast and flexible statistical model for large-scale population genotype data: applications to inferring missing genotypes and haplotypic phase. Am J Hum Genet. [https://doi.org/10.1086/502802](https://doi.org/10.1086/502802)

## Supported data files

The original input format for fastPHASE is outdated. See their original [documentation](http://scheet.org/code/fastphase_doc_1.4.pdf).

fastPHASE.jl supports binary PLINK inputs via [SnpArrays.jl](https://github.com/OpenMendel/SnpArrays.jl), by simply decompressing and converting binary PLINK data into the format that is readable by fastPHASE. Thus, it is not advised to run this on large (e.g. > 1GB PLINK) data. 

Adding VCF, BGEN, and PGEN support in this manner via corresponding OpenMendel modules should not be hard. PRs are welcomed! 

## Installation

fastPHASE.jl requires Julia v1.6 or later. It is not yet registered and can be installed, in the Julia Pkg mode, by

```julia
(@v1.6) Pkg> add https://github.com/biona001/fastPHASE.jl
```

The original fastPHASE software thus fastPHASE.jl only support Linux and MacOS.
