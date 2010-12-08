/**********************************************************************
 * 
 * hmm_mark_cross.c
 * 
 * Based on Karl Broman functions to call calc_genoprob, found in
 * the files hmm_xx.c
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "util.h"
#include "hmm_main.h"
#include "hmm_bc.h"
#include "hmm_f2.h"
#include "hmm_mark.h"
#include "hmm_mark_cross.h"

void calc_markprob_bc(int *n_ind, int *n_mar, int *geno, 
		      double *rf, double *error_prob, double *markprob) 
{
  calc_markprob(*n_ind, *n_mar, 2, geno, rf, rf, *error_prob, markprob,
		init_bc, emit_bc, step_bc);
}

void calc_markprob_f2(int *n_ind, int *n_mar, int *geno, 
                      double *rf, double *error_prob, double *markprob) 
{
  calc_markprob(*n_ind, *n_mar, 3, geno, rf, rf, *error_prob, markprob,
                init_f2, emit_f2, step_f2);
}

/* end of hmm_mark_cross.c */
