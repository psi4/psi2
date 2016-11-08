#include <stdio.h>
#include <libciomr.h>
#define EXTERN
#include "includes.h"
#include "common.h"

/*
** PHASE.C
**
** This routine forms the product: C_new(t)*S*C_old = ~I to check the
** phases of the MOs.  When a phase change occurs, the new MO is
** corrected.  If MO swapping occurs, the phase_check flag is set to zero
** and later written to file30.  This prevents the correlated routines
** from trying to restart from the old wavefunction.
**
** T. Daniel Crawford 6/19/96
**
*/

int phase()
{
  int i,j,k;
  int nn, row_max;
  int phase_chk;
  double maxvalue;
  double **smat, **tmp, **identity;
  double ***cnew;
  struct symm *s;

  phase_chk = 1;
  maxvalue = 0.0;

  /* Make a copy of the new MO vector */
  cnew = (double ***) malloc(num_ir * sizeof(double **));
  for(k=0; k < num_ir; k++) {
      s = &scf_info[k];
      if(nn=s->num_so) {
	  cnew[k] = block_matrix(nn,nn);
	  for(j=0; j < nn; j++)
	      for(i=0; i < nn; i++)
		  cnew[k][j][i] = s->cmat[j][i];
	}
/*      fprintf(outfile, "MOs for Irrep %d\n", k);
      print_mat(cnew[k], nn, nn, outfile); */
    }

  for(k=0; k < num_ir; k++) {
      s = &scf_info[k];
      if(nn=s->num_so) {

	  smat = block_matrix(nn,nn);
	  tmp = block_matrix(nn,nn);
	  identity = block_matrix(nn,nn);

	  /* Unpack the overlap matrix */
	  tri_to_sq(s->smat,smat,nn);

	  /* Form ~I = C^t(new) * S * C(old) */
	  mxmb(cnew[k],nn,1,smat,1,nn,tmp,1,nn,nn,nn,nn);
	  mxmb(tmp,1,nn,s->cmat_orig,1,nn,identity,1,nn,nn,nn,nn);

/*	  fprintf(outfile, "Approximate Identity Matrix for Irrep %d\n", k);
	  print_mat(identity, nn, nn, outfile); */

	  /* Check for MO swapping */
	  for(j=0; j < nn; j++) {
	      maxvalue = 0.0;
	      for(i=0; i < nn; i++) {
		  if(fabs(identity[j][i]) > maxvalue) {
		      maxvalue = fabs(identity[j][i]);
		      row_max = i;
		    }
		}
	      if(row_max != j) phase_chk = 0;
	    }

	  /* Now correct the MO phases, if necessary */
	  if(phase_chk) {
	      for(i=0; i < nn; i++) {
		  if(identity[i][i] < 0.0) {
		      for(j=0; j < nn; j++) cnew[k][j][i] = -(cnew[k][j][i]);
		    }
		}
/*	      fprintf(outfile, "Corrected MOs for irrep %d\n", k);
	      print_mat(cnew[k], nn, nn, outfile); */
	    }
	  free_block(smat);
	  free_block(tmp);
	  free_block(identity);
	}
    }

  /* Finally, put the corrected MOs back into s->cmat if no swapping
     occurred */
  if(phase_chk) {
      fprintf(outfile, "\n Correcting phases of orbitals.\n");
      for(k=0; k < num_ir; k++) {
	  s = &scf_info[k];
	  if(nn=s->num_so) {
	      for(j=0; j < nn; j++)
		  for(i=0; i < nn; i++) s->cmat[j][i] = cnew[k][j][i];
	      free_block(cnew[k]);
	    }
	}
    }
  else fprintf(outfile, "\n No phase correction possible.\n");

  free(cnew);

  return(phase_chk);

}
