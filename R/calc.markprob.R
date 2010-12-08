#####################################################################
#
# calc.markprob.R
#
# copyright (c) 2007-2010, Ani Manichaikul
# last modified August, 2007
# first written August, 2007
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
#
# Contains: calc.markprob
#
######################################################################

######################################################################
#
# calc.markprob: Calculate observed multipoint marker genotype probabilities
# conditional on specified assumed genotypes at each marker and (optionally)
# at positions in between markers
#
######################################################################

calc.markprob <-
function(cross, step=0, off.end=0, error.prob=0.0001,
         map.function=c("haldane","kosambi","c-f","morgan"))
{

  
  # map function
  map.function <- match.arg(map.function)
  if(map.function=="kosambi") mf <- mf.k
  else if(map.function=="c-f") mf <- mf.cf
  else if(map.function=="morgan") mf <- mf.m
  else mf <- mf.h
 
  # don't let error.prob be exactly zero (or >1)
  if(error.prob < 1e-50) error.prob <- 1e-50
  if(error.prob > 1) {
    error.prob <- 1-1e-50
    warning("error.prob shouldn't be > 1!")
  }

  n.ind <- nind(cross)
  n.chr <- nchr(cross)
  n.mar <- nmar(cross)

  type <- class(cross)[1]

  # calculate genotype probabilities one chromosome at a time
  for(i in 1:n.chr) {
    if(n.mar[i]==1) temp.offend <- max(c(off.end,5))
    else temp.offend <- off.end
    
    # which type of cross is this?
    if(type == "f2") {
      one.map <- TRUE
      if(class(cross$geno[[i]]) == "A") { # autosomal
        cfunc <- "calc_markprob_f2"
        n.gen <- 3
        gen.names <- c("A","H","B")
      }
      else {                             # X chromsome 
        cfunc <- "calc_markprob_bc"
        n.gen <- 2
        gen.names <- c("A","H")
      }
    }
    else if(type == "bc") {
      cfunc <- "calc_markprob_bc"
      n.gen <- 2
      gen.names <- c("A","H")
      one.map <- TRUE
    }
    else {
      err <- paste("calc.markprob not available for cross type",
                    type, ".")
      stop(err)
    }

    # genotype data
    gen <- cross$geno[[i]]$data
    gen[is.na(gen)] <- 0
    
    # recombination fractions
    if(one.map) {
      # recombination fractions
      map <- create.map(cross$geno[[i]]$map,step,temp.offend)
      rf <- mf(diff(map))
      rf[rf < 1e-14] <- 1e-14

      # new genotype matrix with pseudomarkers filled in
      newgen <- matrix(ncol=length(map),nrow=nrow(gen))
      dimnames(newgen) <- list(NULL,names(map))
      newgen[,colnames(gen)] <- gen
      newgen[is.na(newgen)] <- 0
      n.pos <- ncol(newgen)
      marnames <- names(map)
    }
    else {
      map <- create.map(cross$geno[[i]]$map,step,temp.offend)
      rf <- mf(diff(map[1,]))
      rf[rf < 1e-14] <- 1e-14
      rf2 <- mf(diff(map[2,]))
      rf2[rf2 < 1e-14] <- 1e-14

      # new genotype matrix with pseudomarkers filled in
      newgen <- matrix(ncol=ncol(map),nrow=nrow(gen))
      dimnames(newgen) <- list(NULL,dimnames(map)[[2]])
      newgen[,colnames(gen)] <- gen
      newgen[is.na(newgen)] <- 0
      n.pos <- ncol(newgen)
      marnames <- colnames(map)
    }

    # call the C function
    if(one.map) {
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(newgen),        # genotype data
              as.double(rf),             # recombination fractions
              as.double(error.prob),     # 
              markprob=as.double(rep(0,n.gen*n.ind*n.pos)))
    }
    else {
      z <- .C(cfunc,
              as.integer(n.ind),         # number of individuals
              as.integer(n.pos),         # number of markers
              as.integer(newgen),        # genotype data
              as.double(rf),             # recombination fractions
              as.double(rf2),            # recombination fractions
              as.double(error.prob),     # 
              markprob=as.double(rep(0,n.gen*n.ind*n.pos)))
    }

    # re-arrange marginal probabilites
    cross$geno[[i]]$markprob <- array(z$markprob,dim=c(n.ind,n.pos,n.gen))
    dimnames(cross$geno[[i]]$markprob) <- list(NULL, marnames, gen.names)
    # attribute set to the error.prob value used, for later
    #     reference, especially by calc.errorlod()
    attr(cross$geno[[i]]$markprob,"error.prob") <- error.prob
    attr(cross$geno[[i]]$markprob,"step") <- step
    attr(cross$geno[[i]]$markprob,"off.end") <- temp.offend
    attr(cross$geno[[i]]$markprob,"map.function") <- map.function
  } # end loop over chromosomes

  cross
}

# end of calc.markprob.R
