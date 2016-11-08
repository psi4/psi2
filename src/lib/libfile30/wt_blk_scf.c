#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_wt_blk_scf():  Writes in the SCF eigenvector (or whatever is to 
**     be stored in its place).
**
**   arguments: 
**
**     int irrep   the number of the irrep to which the symmetry block 
**         belongs (this includes irreps with orbspi[irrep] == 0)
**         n.b. this routine assumes that the first irrep will have
**         irrep == 0.
**
**     double **scf_vector    This should be a single symmetry
**         block of the SCF eigenvector.  Its dimension should be 
**         orbspi[irrep]*orbspi[irrep];
**
**   returns: nothing.
*/

void file30_wt_blk_scf(double **scf_vector, int irrep)
{
  int i,j,k;
  int nmo, mxcoef, nirreps, count, offset;
  int *opi;
  PSI_FPTR junk, mo_coeff_ptr, scf_ptr;
  int tmp;
  double *wt_vector;

  scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 - 1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);
  mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);

  opi = file30_rd_orbspi();

  if(opi[irrep]) 
   {
    offset = 0;
    for(i=0; i < irrep; i++) {
        offset += opi[i]*opi[i];
      } 
    mo_coeff_ptr += (PSI_FPTR) offset * sizeof(double);

    count = 0;
    wt_vector = init_array(opi[irrep]*opi[irrep]);
    for(j=0;j<opi[irrep];j++)
       {
        for(k=0;k<opi[irrep];k++,count++)
          {
           wt_vector[count] = scf_vector[k][j];
          }
       }
  
    wwritw(info30_.filenum, (char *) wt_vector,
          opi[irrep]*opi[irrep]*sizeof(double), mo_coeff_ptr, &junk);

    free(wt_vector);
   }

  free(opi);
} 
  
  

