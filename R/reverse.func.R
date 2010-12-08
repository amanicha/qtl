#####################################################################
#
# scanone.reverse.R
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
# Contains: scanone.reverse, loglik.cross.rev
#
######################################################################

scanone.reverse <- function(cross, chr=1:19, pheno.col = 1, phi.null="em", constrain.lik=FALSE, step=0, maxit = 4000, tol = 1e-04, suppress.warnings=TRUE){
  flag <- FALSE
  if (!("prob" %in% names(cross$geno[[1]])) | !("markprob" %in% names(cross$geno[[1]])))  flag <- TRUE
  else if (dim(cross$geno[[1]]$prob)[2] != dim(cross$geno[[1]]$markprob)[2]) flag <- TRUE

  if(flag){
    if(!suppress.warnings) warning("First calculating genotype and marker probabilities.")
    cross <- calc.genoprob(cross, step=step)
    cross <- calc.markprob(cross, error.prob=0, step=step)
  }

  #determine the appropriate value of phi
  if(phi.null=="seg"){ #use the segregation values
    if("bc" %in% class(cross))
      phi.null <- c(1,1)/2
    else if("f2" %in% class(cross))
      phi.null <- c(1,2,1)/4
  }
  
  #count affecteds and unaffecteds
  na.count <- sum(is.na(cross$pheno[,pheno.col]))
  if(na.count & !suppress.warnings){
    warning(paste("Dropping ",na.count," individuals with phenotype NA."))
    cross <- subset.cross(cross, ind=!is.na(cross$pheno[,pheno.col]))
  }

  ##print warning about likelihood constraint mismatch
  if(!suppress.warnings & !is.numeric(phi.null) & constrain.lik==TRUE) stop("Numerator likelihood constrained without denominator constraint.")

  cross.aff <- subset.cross(cross, ind=(cross$pheno[,pheno.col] == names(table(cross$pheno[,pheno.col]))[1]))

  #check for unaffecteds
  if(length(table(cross$pheno[,pheno.col])) > 1){
    flag.unaff <- TRUE
    cross.unaff <- subset.cross(cross, ind=(cross$pheno[,pheno.col] != names(table(cross$pheno[,pheno.col]))[1]))
  }
  else
    flag.unaff <- FALSE
  
  positions <- sapply(chr, function(i) return(ncol(cross$geno[[i]]$prob)))
    
  loglik.null <- loglik.aff <- loglik.unaff <- matrix(ncol=3, nrow=sum(positions))
  colnames(loglik.null) <- colnames(loglik.aff) <- colnames(loglik.unaff) <- c("chr", "pos", "lod")
  
  index <- 1
                                                    
  for(i in 1:length(chr))
    for(j in 1:positions[i]){
      mymap <- getmap(chr[i], cross)
      pos.here <- mymap[j]
      
      loglik.null[index,] <- c(chr[i], pos.here, loglik.cross.rev(cross=cross, chr=chr[i], pos=j, phi=phi.null, maxit, tol)$loglik )
      likaff.out <- loglik.cross.rev(cross=cross.aff, chr=chr[i], pos=j, phi="em", maxit, tol)

      if(flag.unaff){
        likunaff.out <- loglik.cross.rev(cross=cross.unaff, chr=chr[i], pos=j, phi="em", maxit, tol) 

        #decide if we need to constrain the MLEs
        if(constrain.lik==FALSE || #didn't plan to constrain
           prod(sign((likaff.out$phi - phi.null) * (likunaff.out$phi-phi.null))<=0)==1){ #observed MLEs are already in constrained region
          loglik.aff[index,] <- c(chr[i], pos.here, likaff.out$loglik)
          loglik.unaff[index,] <- c(chr[i], pos.here, likunaff.out$loglik)
        }
        else{ #constrain the likelihood
          likaff.constr <- loglik.cross.rev(cross=cross.aff, chr=chr[i], pos=j, phi=phi.null, maxit, tol)
          likunaff.constr <- loglik.cross.rev(cross=cross.unaff, chr=chr[i], pos=j, phi=phi.null, maxit, tol)

          #maximize across the two possible constraints
          if(likaff.out$loglik + likunaff.constr$loglik >= likaff.constr$loglik + likunaff.out$loglik){
            loglik.aff[index,] <- c(chr[i], pos.here, likaff.out$loglik)
            loglik.unaff[index,] <- c(chr[i], pos.here, likunaff.constr$loglik)
          }
          else{
            loglik.aff[index,] <- c(chr[i], pos.here, likaff.constr$loglik)
            loglik.unaff[index,] <- c(chr[i], pos.here, likunaff.out$loglik)
          }
        }
      }
      else #no unaffecteds
        loglik.aff[index,] <- c(chr[i], pos.here, likaff.out$loglik)
      
      index <- index+1
    }

  results <- loglik.aff
  
  if(flag.unaff)
    results[,"lod"] <- results[,"lod"] + loglik.unaff[,"lod"] - loglik.null[,"lod"]
  else
    results[,"lod"] <- results[,"lod"] - loglik.null[,"lod"]
  
  rownames(results) <- paste("c", rep(chr, positions), ".", unlist(lapply(cross$geno[chr], function(geno) return(colnames(geno$prob)))),sep='')

  results <- as.data.frame(results)
  class(results) <- c("scanone", "data.frame")
  return(results)
}

#################################
# general functions             #
#################################

loglik.cross.rev <- function(cross, chr, pos, phi="em", ...){
  markprob.pos <- cross$geno[[chr]]$markprob[,pos,]

  if(phi[1]=="em")
    phi <- em.reverse(cross, chr, pos, ...)
  
  ll <- sum(log10(markprob.pos%*%phi) )
  
  return(list(loglik=ll, phi=phi))
}

em.reverse <- function(cross, chr, pos, maxit, tol){
  markprob.pos <- cross$geno[[chr]]$markprob[,pos,]
  genoprob.pos <- cross$geno[[chr]]$prob[,pos,]
  phi.s <- apply(genoprob.pos,2,mean)

  flag <- 1
  count <- 1
  
  while(flag){
    #E step
    marginal.prob <- markprob.pos%*%phi.s
    joint.prob <- sweep(markprob.pos,2, phi.s,"*")
    e.s1 <- apply(sweep(joint.prob,1, marginal.prob,"/"),2,sum)
    
    #M step
    phi.s1 <- e.s1 / nind(cross)
    
    #flag update
    if(count>=maxit| max(abs(phi.s1 - phi.s))<tol)
      flag=0
    else
      count <- count+1

    phi.s <- phi.s1
  }

  return(phi.s)
}

