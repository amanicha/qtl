\name{scanone.full}
\alias{scanone.full}
\title{ Genome scan for binary trait with selective genotyping, using full likelihood approach}
\description{
This function implements the single QTL genome scan using the full-likelihood approach as described for binary traits with selective genotyping, in Manichaikul and Broman (2009).  
}
\usage{
scanone.full(cross, chr = 1:19, pheno.col = 1, step = 0, maxit = 4000, tol = 0.0001, suppress.warnings = TRUE)
}
\arguments{
  \item{cross}{ An object of class 'cross'. See 'read.cross' for details. }
  \item{chr}{ Optional vector indicating the chromosomes for which LOD scores should be calculated. }
  \item{pheno.col}{ Column number in the phenotype matrix which should be used as the phenotype. }
  \item{step}{ Step size to be passed to 'calc.markprob' if the marker probabilities are not found in the 'cross' object that is passed to the 'scanone.reverse' function. }
  \item{maxit}{ Maximum number of iterations for EM algorithm. }
  \item{tol}{ Tolerance value for determining convergence of EM algorithm. }
  \item{suppress.warnings}{ Flag indicating whether warnings should be displayed }
}
\details{ This approach is implemented imposing a constraint on the likelihood maximization, with segregation proportions set to c(1,1)/2 for a backcross and c(1,2,1)/4 for an intercross.  If the constraint is not necessary, it is equivalent to use the "reverse" approach as implemented in \code{scanone.reverse}.
}
\references{
Manichaikul A, Broman KW (2009) Binary trait mapping in experimental crosses with selective genotyping. Genetics 182: 863-74.
}
\author{
Ani Manichaikul \email{amanicha@virginia.edu}
}
\seealso{\code{\link[qtl]{scanone}},\code{\link{scanone.reverse}}}
\keyword{models}
\keyword{ optimize}