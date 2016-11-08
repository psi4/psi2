#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_usotao():  Reads in the untary SO to AO transformation matrix.
**
**   takes no arguments.
**
**   returns: double **usotao   The unitary SO to AO matrix.  It is a matrix
**     of dimentions nmo*nao (there are nmo rows of length nao).
*/

double **file30_rd_usotao(void)
{
  int i, natom, nmo, nao;
  PSI_FPTR usotao_ptr;
  double **usotao;

  nao = file30_rd_nao();
  nmo = file30_rd_nmo();
  natom = file30_rd_natom();
  usotao_ptr = (PSI_FPTR) 
        ((info30_.mcalcs[0]+60+20+natom*6-1)*sizeof(int)+10*sizeof(double));

  usotao = init_matrix(nmo,nao);

  for(i=0;i<nmo;i++)
   {
    wreadw(info30_.filenum, (char *) usotao[i], (int)sizeof(double)*nao,
	 usotao_ptr, &usotao_ptr);
   }

  return usotao;
}
