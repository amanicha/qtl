\name{scanone.reverse}
\alias{scanone.reverse}
\title{ Genome scan for binary trait with selective genotyping, conditioning genotypes on phenotypes}
\description{ This function implements the single QTL genome scan using the "reverse" approach as described for binary traits with selective genotyping, in Manichaikul and Broman (2009).  This approach is implemented with the option of imposing a constraint on the likelihood maximization.
}
\usage{
scanone.reverse(cross, chr = 1:19, pheno.col = 1, phi.null = "em", constrain.lik = FALSE, step = 0, maxit = 4000, tol = 0.0001, suppress.warnings = TRUE)
}
\arguments{
  \item{cross}{ An object of class 'cross'. See 'read.cross' for details. }
  \item{chr}{ Optional vector indicating the chromosomes for which LOD scores should be calculated. }
  \item{pheno.col}{ Column number in the phenotype matrix which should be used as the phenotype. }
  \item{phi.null}{ Expected genotype proportions under the null hypothesis.  This can be a numeric vector of length =2 for a backcross, or =3 for an intercross.  The default is to set 'phi.null'="em", in which case we estimate the genotype proportions for each genetic position using the available genotype data.  Alternatively, 'phi.null' can be set to "seg", in which case we set it to c(1,1)/2 for a backcross, and c(1,2,1)/4 for an intercross.   }
  \item{constrain.lik}{ Should the alternative hypothesis likelihood be maximized subject to a constraint of the specified value of 'phi.null'?  This option should always be set to 'FALSE' for affected-only analysis. }
  \item{step}{ Step size to be passed to 'calc.markprob' if the marker probabilities are not found in the 'cross' object that is passed to the 'scanone.reverse' function. }
  \item{maxit}{ Maximum number of iterations for EM algorithm. }
  \item{tol}{ Tolerance value for determining convergence of EM algorithm. }
  \item{suppress.warnings}{ Flag indicating whether warnings should be displayed }
}
\references{
Manichaikul A, Broman KW (2009) Binary trait mapping in experimental crosses with selective genotyping. Genetics 182: 863-74.
}
\author{
Ani Manichaikul \email{amanicha@virginia.edu}
}
\seealso{\code{\link[qtl]{scanone}},\code{\link{scanone.full}}}
\keyword{models}
\keyword{ optimize}