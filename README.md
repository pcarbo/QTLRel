QTLRel
======

###Overview

QTLRel is an R package for quantitative trait mapping in populations
such as advanced intercross lines where relatedness among individuals
should not be ignored. QTLRel includes functions to estimate
background genetic variance components, impute missing genotypes,
simulate genotypes, perform a genome scan for quantitative trait loci
(QTLs), and plot mapping results. QTLRel also includes functions to
efficiently calculate identity coefficients.

Mention here that anyone who finds this useful should cite a couple
publications.

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
