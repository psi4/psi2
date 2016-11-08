#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_blk_scf():  Reads in a symmetry block of the SCF eigenvector 
**   (or whatever is stored in its place).
**  
**   arguments: int irrep   designates the particular irrep to which the block
**     belongs.
**  
**   returns: double **scf_vector    This square matrix has orbspi[irrep] 
**     rows.
*/


double **file30_rd_blk_scf(int irrep)
{
  int i,j;
  double **scf_vector;
  int nmo, mxcoef, nirreps, count, num_orbs;
  int sym, first, last, column, offset;
  int *opi;
  PSI_FPTR junk, mo_coeff_ptr, scf_ptr;
  int tmp;
  double *tmp_vector;

  scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 - 1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);
  mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);

  nirreps = file30_rd_nirreps();
  mxcoef = file30_rd_mxcoef();
  opi = file30_rd_orbspi();

  scf_vector = NULL;

  if(opi[irrep])
   {
    tmp_vector = init_array(opi[irrep]*opi[irrep]);

    offset = 0;
    for(i=0; i < irrep; i++) {
        offset += opi[i]*opi[i];
      }
    mo_coeff_ptr += (PSI_FPTR) offset*sizeof(double);
  
    wreadw(info30_.filenum, (char *) tmp_vector, 
           (opi[irrep]*opi[irrep])*sizeof(double), mo_coeff_ptr, &junk);
  
    scf_vector = init_matrix(opi[irrep],opi[irrep]);
    
    count = 0;
    for(i=0; i < opi[irrep] ; i++)
       for(j=0; j < opi[irrep] ; j++, count++) {
          scf_vector[j][i] = tmp_vector[count];
         }
    free(tmp_vector);
    }

  free(opi);

  return scf_vector;
}
