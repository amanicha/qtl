/**********************************************************************
 * 
 * hmm_mark.c
 *
 *
 * Based on Karl Broman's hmm_main.c for the R/qtl package
 * Modified by Ani Manichaikul
 *
 * last modified June, 2007
 * first written June, 2007
 *
 *
 * Contains: calc_markprob
 *  
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_main.h"
#include "hmm_mark.h"
#include "util.h"

/**********************************************************************
 * 
 * calc_markprob
 *
 # This function uses hidden Markov models to calculate 
 # the probability of observed multipoint marker data for each individual, 
 # conditional on assumed genotype values at each marker and (optionally) 
 # at points between markers on a chromosome.
 * This function assumes data on a single chromosome
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of positions at which to 
 *              calculate the conditional observed marker probabilities
 *
 * n_gen        Number of different genotypes
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * rf2          A second set of recombination fractions, in case of
 *              sex-specific maps (may be ignored)
 *
 * error_prob   Genotyping error probability
 *
 * markprob     Observed marker probabilities conditional on genotypes 
 *              at the specified positions (the output); a single vector
 *              stored by columns (ind moves fastest, then mar, then
 *              genotype
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void calc_markprob(int n_ind, int n_pos, int n_gen, int *geno, 
		   double *rf, double *rf2, 
		   double error_prob, double *markprob, 
		   double initf(int),
		   double emitf(int, int, double),
		   double stepf(int, int, double, double)) 
{
  int i, j, j2, v, v2;
  double **betal, **betar; /* betas for left and right sides of the chromosome */
  int **Geno;
  double ***Markprob;
  
  /* allocate space for betal and betar and 
     reorganize geno and markprob */
  reorg_geno(n_ind, n_pos, geno, &Geno);
  reorg_genoprob(n_ind, n_pos, n_gen, markprob, &Markprob);
  allocate_alpha(n_pos, n_gen, &betal);
  allocate_alpha(n_pos, n_gen, &betar);

  for(i=0; i<n_ind; i++) { /* i = individual */

    /* initialize betal and betar */
    for(v=0; v<n_gen; v++) {
      betal[v][0] = 0.0;
      betar[v][n_pos-1] = 0.0;
    }

    /* backward equations */
    for(j=1,j2=n_pos-2; j<n_pos; j++, j2--) {
      
      for(v=0; v<n_gen; v++) {
	betal[v][j] = betal[0][j-1] + stepf(v+1, 1, rf[j-1], rf2[j-1]) +
	  emitf(Geno[j-1][i],1,error_prob);
	
	betar[v][j2] = betar[0][j2+1] + stepf(v+1,1,rf[j2], rf2[j2]) + 
	  emitf(Geno[j2+1][i],1,error_prob);

	for(v2=1; v2<n_gen; v2++) {
	  betal[v][j] = addlog(betal[v][j], betal[v2][j-1] + 
			       stepf(v+1,v2+1,rf[j-1],rf2[j-1]) +
			       emitf(Geno[j-1][i],v2+1,error_prob));
	  betar[v][j2] = addlog(betar[v][j2], betar[v2][j2+1] + 
			       stepf(v+1,v2+1,rf[j2],rf2[j2]) +
			       emitf(Geno[j2+1][i],v2+1,error_prob));
	}

      }
    }

    /* calculate genotype probabilities */
    for(j=0; j<n_pos; j++) 
      for(v=0; v<n_gen; v++) 
	Markprob[v][j][i] = exp(betal[v][j] + betar[v][j] + emitf(Geno[j][i], v+1, error_prob)); 
      
  } /* loop over individuals */
  
}



/* end of hmm_mark.c */
