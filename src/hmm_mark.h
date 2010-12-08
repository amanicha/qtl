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


/**********************************************************************
 * 
 * calc_markprob
 *
 * This function uses the hidden Markov model technology to calculate 
 * the observed marker probabilities, conditional on assumed genotype values
 * at each of marker and (optionally) at points between markers on a chromosome.
 * This assumes data on a single chromosome
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
		   double stepf(int, int, double, double)); 
