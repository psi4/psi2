#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_evals():  Reads in the SCF eigenvalues: the SCF orbital energies.
**
**  takes no arguments.
**
**  returns: double *evals   an array of _all_ of the SCF eigenvalues,
**      ordered by irrep, and by increasing energy within each irrep.  
**      (i.e. for sto water, the four a1 eigenvalues all come first, and 
**      those four are ordered from lowest energy to highest energy,
**      followed by the single b1 eigenvalue, etc.)
*/

double *file30_rd_evals(void)
{
  PSI_FPTR scf_ptr, evalsptr, mo_coeff_ptr, junk;
  int tmp;
  double *energies;
  
  scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 -1)*sizeof(int);

  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);

  mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);
  evalsptr = mo_coeff_ptr + (PSI_FPTR) file30_rd_mxcoef()*sizeof(double);

  energies = init_array(file30_rd_nmo());
  wreadw(info30_.filenum, (char *) energies, file30_rd_nmo()*sizeof(double),
	 evalsptr, &junk);

  return energies;

}
