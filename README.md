QTLRel
======

###Overview

QTLRel is an R package for quantitative trait mapping in populations
such as advanced intercross lines (AILs) where relatedness among
individuals should not be ignored. QTLRel includes functions to
estimate background genetic variance components, impute missing
genotypes, simulate genotypes, perform a genome scan for quantitative
trait loci (QTLs), and plot the mapping results. QTLRel also includes
functions to efficiently calculate identity coefficients, including
kinship coefficients.

This package implements the methods described in

> Cheng R, Lim J E, Samocha K E, Sokoloff G, Abney M, Skol A D,
> Palmer A A (2010).
> [Genome-wide association studies and the problem of relatedness
among advanced intercross lines and other highly recombinant
populations](http://dx.doi.org/10.1534/genetics.110.116863)
> *Genetics* 185(3): 1033–1044.

If you find QTLRel useful for your own project, we request that you
cite *Genetics* (2010) paper above, and the following paper:

> Cheng R, Abney M, Palmer A A, Skol A D (2011). [QTLRel: an R
package for genome-wide association studies in which relatedness is a
concern](http://dx.doi.org/10.1186/1471-2156-12-66).
> *BMC Genetics* 12(1): 66. 

###License

Put a summary of the GNU Public License (version 3) here.

###Installation

Explain here how to install the R package, either through CRAN (see
Section 2 of the tutorial PDF), or by downloading the source code from
github and installing the package using the local files.

###More information

Point to the documentation for the individual functions, and the
QTLRel tutorial for a detailed introduction.

Also point out that a very small change was made to this version of
the code to remove checks on the pedigree preventing self-mating
(selfing). Since the computation of identity coefficients assume that
the founders are *not* inbred, removing this check allows one to
define (mostly) inbred individual through iterated self-mating.

###Credits

Mention Riyan Cheng, developed at Abe Palmer's lab at the University
of Chicago. Also mention me---I'm the one who this code on Github.
