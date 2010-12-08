\name{qtlMprob-internal}
\alias{qtlMprob-internal}
\alias{em.full}
\alias{em.reverse}
\alias{getmap}
\alias{grad.bc}
\alias{grad.f2}
\alias{loglik.cross.full}
\alias{loglik.cross.rev}
\alias{neg.ll.bc}
\alias{neg.ll.f2}
\alias{my.subset.cross}
\title{Internal functions for qtlMprob}
\description{
Internal functions for binary trait mapping with selective genotyping, as implemented in scanone.full and scanone.reverse.  These functions are not intended to be called by the user.
}
\usage{
em.full(cross.aff, cross.unaff, chr, pos, naff, nunaff, pheno.col, phi.null, maxit, tol)
em.reverse(cross, chr, pos, maxit, tol)
getmap(i, cross)
grad.bc(params, e.aff, e.unaff, phi.null)
grad.f2(params, e.aff, e.unaff, phi.null)
loglik.cross.full(cross.aff, cross.unaff, chr, pos, param = "em", naff = naff, nunaff = nunaff, ...)
loglik.cross.rev(cross, chr, pos, phi = "em", ...)
neg.ll.bc(params, e.aff, e.unaff, phi.null)
neg.ll.f2(params, e.aff, e.unaff, phi.null)
my.subset.cross(x, chr, ind, ...)
}
\author{
Ani Manichaikul \email{amanicha@virginia.edu}
}
\keyword{internal}