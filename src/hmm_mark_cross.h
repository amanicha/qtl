/**********************************************************************
 * 
 * hmm_mark_cross.h
 * 
 * Based on Karl Broman functions to call calc_genoprob, found in
 * the files hmm_xx.h
 * Modified by Ani Manichaikul
 *
 * last modified Aug, 2007
 * first written Aug, 2007
 *
 *
 * Contains: calc_markprob_bc, calc_markprob_f2
 *
 * These are hmm wrappers for backcross and intercross.
 *
 * BC Genotype codes:  0=AA; 1=AB
 * F2 Genotype codes:  0=AA; 1=AB; 2=BB
 * Phenotype codes: 0=missing; 1=AA; 2=AB
 *
 **********************************************************************/

void calc_markprob_bc(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *markprob); 


void calc_markprob_f2(int *n_ind, int *n_mar, int *geno, 
                      double *rf, double *error_prob, double *markprob); 


/* end of hmm_mark_cross.h */
