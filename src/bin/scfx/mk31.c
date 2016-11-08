
/*
 * this is a quick hack to create file31
 */

#include <stdio.h>

#if !defined(AIX)
#  define swrit swrit_
#  define prints prints_
#endif

#define ioff(i) ((i)*((i)+1)/2)

void
#if defined(AIX)
mk31(double *cdens, double *odens, int *nsym, int *nlamda, int *dc, int *ddo)
#else
mk31_(double *cdens, double *odens, int *nsym, int *nlamda, int *dc, int *ddo)
#endif
{
  double *tempc,*tempo;
  int nbasis=0,i,j,k,nn,loff;
  int ntri;

  for (i=0; i < *nsym; i++) nbasis+=nlamda[i];

  ntri = nbasis*(nbasis+1)/2;

  tempc = (double*) malloc(sizeof(double)*ntri);
  tempo = (double*) malloc(sizeof(double)*ntri);

  if(!tempc || !tempo) {
    fprintf(stderr,"mk31: malloc failed for density matrices\n");
    fprintf(stderr,"tempc = 0x%x, tempo = 0x%x\n",tempc,tempo);
    exit(1);
    }

  bzero(tempc,sizeof(double)*nbasis);
  bzero(tempo,sizeof(double)*nbasis);

  for (i=loff=0; i < *nsym; i++) {
    if(nn=nlamda[i]) {
      for(j=0; j < nn; j++)
        for (k=0; k <= j; k++) {
          tempc[ioff(j+loff)+k+loff]  = cdens[ioff(j)+k+loff];
          tempo[ioff(j+loff)+k+loff]  = odens[ioff(j)+k+loff];
          }
      loff += nn;
      }
    }

#if 0
  nn=3;
  prints(tempc,&ntri,&nbasis,&nn);
  prints(tempo,&ntri,&nbasis,&nn);
#endif

  nn=31; i=sizeof(double)/sizeof(int)*ntri;
  swrit(&nn,tempc,&i);
  swrit(&nn,tempo,&i);


  free(tempc);
  free(tempo);
  }
