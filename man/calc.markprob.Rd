\name{calc.markprob}
\alias{calc.markprob}
\title{Calculate observed multipoint marker genotype probabilities conditional on specified assumed genotypes at each marker and (optionally) at positions in between markers
}
\description{
This function uses hidden Markov models to calculate 
the probability of observed multipoint marker data for each individual, conditional on assumed genotype values at each marker and (optionally) at points between markers on a chromosome.}
\usage{
calc.markprob(cross, step = 0, off.end = 0, error.prob = 0.0001, map.function = c("haldane", "kosambi", "c-f", "morgan"))
}
\arguments{
  \item{cross}{ An object of class 'cross'. See 'read.cross' for details.}
  \item{step}{  Maximum distance (in cM) between positions at which the genotype probabilities are calculated, though for 'step = 0', probabilities are calculated only at the marker locations. }
  \item{off.end}{ Distance (in cM) past the terminal markers on each chromosome to which the genotype probability calculations will be carried. }
  \item{error.prob}{ Assumed genotyping error rate used in the calculation of the penetrance Pr(observed genotype | true genotype). }
  \item{map.function}{ Indicates whether to use the Haldane, Kosambi or Carter-Falconer map function when converting genetic distances into recombination fractions. }
}
\value{ The input \code{cross} object is returned with a component, \code{markprob}, added to each component of \code{cross$geno}. code{markprob} is an array of size [n.ind x n.pos x n.gen] where n.pos is the number of positions at which the probabilities were calculated and n.gen = 3 for an intercross, and = 2 for a backcross. Attributes \code{"error.prob"}, \code{"step"}, \code{"off.end"}, and code{"map.function"} are set to the values of the corresponding arguments, for later reference.
}
\references{
Manichaikul A, Broman KW (2009) Binary trait mapping in experimental crosses with selective genotyping. Genetics 182: 863-74.
}
\author{
Ani Manichaikul \email{amanicha@virginia.edu}
}
\keyword{models}
\seealso{\code{\link[qtl]{calc.genoprob}}}