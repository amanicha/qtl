#####################################################################
#
# scanone.full.R
#
# copyright (c) 2007-2010, Ani Manichaikul
# last modified Feb, 2010
# first written 2007
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
# Contains: scanone.full, loglik.cross.full, neg.ll.bc, grad.bc, neg.ll.f2, grad.f2, em.full
#
######################################################################

#full likelihood with constraint on segregation probabilities
#assumes we have both affected and unaffecteds, otherwise scanone.reverse should be used
scanone.full <- function(cross, chr=1:19, pheno.col = 1, step=0, maxit = 4000, tol = 1e-04, suppress.warnings=TRUE){
  flag <- FALSE
  if (!("prob" %in% names(cross$geno[[1]])) | !("markprob" %in% names(cross$geno[[1]])))  flag <- TRUE
  else if (dim(cross$geno[[1]]$prob)[2] != dim(cross$geno[[1]]$markprob)[2]) flag <- TRUE

  if(flag){
    if(!suppress.warnings) warning("First calculating genotype and marker probabilities.")
    cross <- calc.genoprob(cross, step=step)
    cross <- calc.markprob(cross, error.prob=0, step=step)
  }

  #determine the appropriate value of phi
  if("bc" %in% class(cross))
    phi.null <- c(1,1)/2
  else if("f2" %in% class(cross))
    phi.null <- c(1,2,1)/4

  #count affecteds and unaffecteds
  na.count <- sum(is.na(cross$pheno[,pheno.col]))
  if(na.count & !suppress.warnings){
    warning(paste("Dropping ",na.count," individuals with phenotype NA."))
    cross <- subset.cross(cross, ind=!is.na(cross$pheno[,pheno.col]))
  }
  
  dis.tab <- table(cross$pheno[,pheno.col])
  nunaff <- dis.tab[1]
  naff <- dis.tab[2]

  if(length(dis.tab)!=2) stop("Full likelihood is not applicable here.  The number of phenotypes should be exactly two.")  
  
  pi.null <- naff / (naff + nunaff)

  param.null <- list(phi.aff=phi.null, phi.unaff=phi.null, pi=pi.null)

  #force phenotypes to be 0 / 1
  cross$pheno[,pheno.col] <- ifelse(cross$pheno[,pheno.col]==names(dis.tab)[2], 1, 0)

  #subset crosses for affected and unaffected
  cross.aff <- subset.cross(cross, ind=(cross$pheno[,pheno.col] == 1))
  cross.unaff <- subset.cross(cross, ind=(cross$pheno[,pheno.col] == 0))

  positions <- sapply(chr, function(i) return(ncol(cross$geno[[i]]$prob)))
    
  loglik.null <- loglik.alt <- matrix(ncol=3, nrow=sum(positions))
  colnames(loglik.null) <- colnames(loglik.alt) <- c("chr", "pos", "lod")
  
  index <- 1
  
  for(i in 1:length(chr))
    for(j in 1:positions[i]){
      mymap <- getmap(chr[i], cross)
      pos.here <- mymap[j]
      
      loglik.null[index,] <- c(chr[i], pos.here, loglik.cross.full(cross.aff=cross.aff, cross.unaff=cross.unaff,
                                            chr=chr[i], pos=j, param=param.null, naff=naff, nunaff=nunaff)$loglik )
      loglik.alt[index,] <- c(chr[i], pos.here, loglik.cross.full(cross.aff=cross.aff, cross.unaff=cross.unaff,
                                           chr=chr[i], pos=j, param="em", naff=naff, nunaff=nunaff, pheno.col, phi.null, maxit, tol)$loglik )

      index <- index+1
    }

  results <- loglik.alt
  results[,"lod"] <- results[,"lod"] - loglik.null[,"lod"]
  
  rownames(results) <- paste("c", rep(chr, positions), ".", unlist(lapply(cross$geno[chr], function(geno) return(colnames(geno$prob)))),sep='')

  results <- as.data.frame(results)
  class(results) <- c("scanone", "data.frame")
  return(results)
}

#################################
# general functions             #
#################################

loglik.cross.full <- function(cross.aff, cross.unaff, chr, pos, param="em", naff=naff, nunaff=nunaff, ...){
  markprob.pos.aff <- cross.aff$geno[[chr]]$markprob[,pos,]
  markprob.pos.unaff <- cross.unaff$geno[[chr]]$markprob[,pos,]

  if(param[1]=="em")
    param <- em.full(cross.aff, cross.unaff, chr, pos, naff, nunaff, ...)

  phi.aff <- param[["phi.aff"]]
  phi.unaff <- param[["phi.unaff"]]
  pi <- param[["pi"]]
  
  ll <- sum(log10(markprob.pos.aff%*%phi.aff) ) + sum(log10(markprob.pos.unaff%*%phi.unaff) ) + log10(pi) * naff + log10(1-pi) * nunaff
  
  return(list(loglik=ll, param=param))
}

##########################
# to be used in optim    #
##########################
neg.ll.bc <- function(params, e.aff, e.unaff, phi.null){
  phi.aa.aff <- params[1]
  pi <- params[2]

  phi.aa.unaff <- (phi.null[1] - phi.aa.aff*pi) / (1-pi)
    
  A <- e.aff[1]; B <- e.aff[2]
  C <- e.unaff[1]; D <- e.unaff[2]

  if(max(phi.aa.aff,phi.aa.unaff,pi)>=1 | min(phi.aa.aff, phi.aa.unaff, pi)<=0)
    ll <- -10^10
  else
    ll <- A * (log(phi.aa.aff) + log(pi)) + B * (log(1-phi.aa.aff) + log(pi)) +
      C * (log(phi.aa.unaff) + log(1-pi)) + D*(log(1-phi.aa.unaff) + log(1-pi))

  return(-ll)
}

grad.bc <- function(params, e.aff, e.unaff, phi.null){
  phi.aa.aff <- params[1]
  pi <- params[2]

  A <- e.aff[1]; B <- e.aff[2]
  C <- e.unaff[1]; D <- e.unaff[2]

  deriv.phi.aa.aff <- A/phi.aa.aff - B/(1-phi.aa.aff) -
    C * pi / (phi.null[1] - phi.aa.aff*pi) + D * pi / (phi.null[2] - (1- phi.aa.aff)*pi)

  deriv.pi <- (A + B)/pi - C * phi.aa.aff / (phi.null[1] - phi.aa.aff*pi) - D*(1-phi.aa.aff) / (phi.null[2] - (1-phi.aa.aff)*pi)
  
  return(c(-deriv.phi.aa.aff, -deriv.pi))  
}

neg.ll.f2 <- function(params, e.aff, e.unaff, phi.null){
  phi.aa.aff <- params[1]
  phi.ab.aff <- params[2]
  pi <- params[3]

  phi.aa.unaff <- (phi.null[1] - phi.aa.aff*pi) / (1-pi)
  phi.ab.unaff <- (phi.null[2] - phi.ab.aff*pi) / (1-pi)

  A <- e.aff[1]; B <- e.aff[2]; C <- e.aff[3]
  D <- e.unaff[1]; E <- e.unaff[2]; F <- e.unaff[3]

  if(max(phi.aa.aff, phi.ab.aff, phi.aa.unaff, phi.ab.unaff,pi) >= 1 | min(phi.aa.aff, phi.ab.aff, phi.aa.unaff,phi.ab.unaff, pi)<=0)
    ll <- -10^10
  else
    ll <- A * (log(phi.aa.aff) + log(pi)) + B * (log(phi.ab.aff) + log(pi)) + C * (log(1-phi.aa.aff-phi.ab.aff) + log(pi)) +
      D * (log(phi.aa.unaff) + log(1-pi)) + E*(log(phi.ab.unaff) + log(1-pi)) + F*(log(1- phi.aa.unaff - phi.ab.unaff) + log(1-pi))
  
  return(-ll)
}

grad.f2 <- function(params, e.aff, e.unaff, phi.null){
  phi.aa.aff <- params[1]
  phi.ab.aff <- params[2]
  pi <- params[3]

  A <- e.aff[1]; B <- e.aff[2]; C <- e.aff[3]
  D <- e.unaff[1]; E <- e.unaff[2]; F <- e.unaff[3]
  
  deriv.phi.aa.aff <- A/phi.aa.aff - C/(1-phi.aa.aff-phi.ab.aff) -
    D * pi / (phi.null[1] - phi.aa.aff*pi) + F * pi / (phi.null[3] - (1- phi.aa.aff-phi.ab.aff)*pi)
  deriv.phi.ab.aff <- B/phi.aa.aff - C/(1-phi.aa.aff-phi.ab.aff) -
    E * pi / (phi.null[2] - phi.ab.aff*pi) + F * pi / (phi.null[3] - (1- phi.aa.aff-phi.ab.aff)*pi)
  
  deriv.pi <- (A + B + C)/pi - D * phi.aa.aff / (phi.null[1] - phi.aa.aff*pi) -
    E * phi.ab.aff / (phi.null[2] - phi.ab.aff*pi) - F* (1-phi.aa.aff-phi.ab.aff) / (phi.null[3] - (1-phi.aa.aff-phi.ab.aff)*pi)
  
  return(c(-deriv.phi.aa.aff, -deriv.phi.ab.aff, -deriv.pi))
    
}

############
em.full <- function(cross.aff, cross.unaff, chr, pos, naff, nunaff, pheno.col, phi.null, maxit, tol){
  markprob.pos.aff <- cross.aff$geno[[chr]]$markprob[,pos,]
  markprob.pos.unaff <- cross.unaff$geno[[chr]]$markprob[,pos,]

  genoprob.pos.aff <- cross.aff$geno[[chr]]$prob[,pos,]
  genoprob.pos.unaff <- cross.unaff$geno[[chr]]$prob[,pos,]

  #initial guess values
  pi.s <- naff / (naff + nunaff)
  phi.s.aff <- apply(genoprob.pos.aff,2,mean)
  phi.s.unaff <- (phi.null - phi.s.aff *pi.s) / (1-pi.s)
  if(max(phi.s.unaff)>=1 | min(phi.s.unaff)<=0)
    phi.s.aff <- phi.s.unaff <- phi.null
  
  flag <- 1
  count <- 1
  
  while(flag){
    #E step
    marginal.prob.aff <- markprob.pos.aff%*%phi.s.aff
    joint.prob.aff <- sweep(markprob.pos.aff,2, phi.s.aff,"*")
    e.s1.aff <- apply(sweep(joint.prob.aff,1, marginal.prob.aff,"/"),2,sum)

    marginal.prob.unaff <- markprob.pos.unaff%*%phi.s.unaff
    joint.prob.unaff <- sweep(markprob.pos.unaff,2, phi.s.unaff,"*")
    e.s1.unaff <- apply(sweep(joint.prob.unaff,1, marginal.prob.unaff,"/"),2,sum)

    #M step
    if("bc" %in% class(cross.aff))
      optim.out <- optim(c(phi.s.aff[1],pi.s), neg.ll.bc, gr=NULL, 
                         e.aff=e.s1.aff, e.unaff=e.s1.unaff, phi.null=phi.null)$par
    else if("f2" %in% class(cross.aff))
      optim.out <- optim(c(phi.s.aff[1:2],pi.s), neg.ll.f2, gr=NULL,
                         e.aff=e.s1.aff, e.unaff=e.s1.unaff, phi.null=phi.null)$par

    pi.s1 <- optim.out[length(optim.out)]
    phi.s1.aff <- optim.out
    phi.s1.aff[length(phi.s1.aff)] <- 1 - sum(phi.s1.aff[-length(phi.s1.aff)])

        
    #flag update
    if(count>=maxit| (max(abs(phi.s1.aff - phi.s.aff))<tol & abs(pi.s1 - pi.s)<tol) )
      flag=0
    else
      count <- count+1

    pi.s <- pi.s1
    phi.s.aff <- phi.s1.aff
    phi.s.unaff <- (phi.null - phi.s.aff *pi.s) / (1-pi.s) 
  }

  return(list(phi.aff=phi.s.aff, phi.unaff=phi.s.unaff, pi=pi.s))
}

