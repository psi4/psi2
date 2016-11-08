/* $Log: diis.c,v $
/* Revision 1.1  1991/06/15 20:22:20  seidl
/* Initial revision
/* */

static char *rcsid = "$Id: diis.c,v 1.1 1991/06/15 20:22:20 seidl Exp $";

#define EXTERN
#include "includes.h"
#include "common.h"

extern double delta;

static double *btemp, **bold, **bmat;
static struct diis_mats {
  double ***fock_c;
  double ***fock_o;
  double ***error;
  } *diism,dtemp;

void diis(scr1,scr2,scr3,c1,c2,cim,newci)
   double cim;
   double **scr1, **scr2, **scr3,*c1,*c2;
   int newci;
{
   int i,j,k,ij;
   int m,nn;
   double occi, occj;
   int try = 0;
   int last = iter-1;
   int col = iter+1;
   double etemp, dotp, norm, determ, etempo;
   double scale;
   struct symm *s;
   struct diis_mats *d;

   if(diism == NULL) {
      bmat = (double **) init_matrix(ndiis+1,ndiis+1);
      bold = (double **) init_matrix(ndiis,ndiis);
      btemp = (double *) init_array(ndiis+1);

      diism = (struct diis_mats *) malloc(sizeof(struct diis_mats)*ndiis);

      for(i=0; i < ndiis ; i++) {
         diism[i].fock_c = (double ***) malloc(sizeof(double **)*num_ir);
         if(iopen) 
            diism[i].fock_o = (double ***) malloc(sizeof(double **)*num_ir);
         diism[i].error = (double ***) malloc(sizeof(double **)*num_ir);
         for(j=0; j < num_ir ; j++) {
            if(nn=scf_info[j].num_so) {
               diism[i].fock_c[j] = (double **) init_matrix(nn,nn);
               if(iopen) diism[i].fock_o[j] = (double **) init_matrix(nn,nn);
               diism[i].error[j] = (double **) init_matrix(nn,nn);
               }
            }
         }
      }

   scale = 1.0 + dampsv;

   if (iter > ndiis) {
      last = ndiis-1;
      col = ndiis+1;
      dtemp = diism[0];
      for (i=0; i < last ; i++) {
         diism[i] = diism[i+1];
         }
      diism[last] = dtemp;
      }
      
 /* save ao fock matrices in fock_save */
   
   d = &diism[last];
   for (m=0; m < num_ir ; m++) {
      s = &scf_info[m];
      if(nn=s->num_so) {
         tri_to_sq(s->fock_pac,d->fock_c[m],nn);
         if(iopen) tri_to_sq(s->fock_open,d->fock_o[m],nn);

      /* form error matrix in mo basis */
         mxmb(s->cmat,nn,1,d->fock_c[m],1,nn,scr1,1,nn,nn,nn,nn);
         mxmb(scr1,1,nn,s->cmat,1,nn,scr2,1,nn,nn,nn,nn);
         if(iopen) {
            mxmb(s->cmat,nn,1,d->fock_o[m],1,nn,scr1,1,nn,nn,nn,nn);
            mxmb(scr1,1,nn,s->cmat,1,nn,scr3,1,nn,nn,nn,nn);
            }

         for (i=0; i < nn; i++) {
            occi = s->occ_num[i];
            for (j=0; j <= i ; j++ ) {
               occj = s->occ_num[j];
               if (!iopen) {
                  if (occi == 0.0 && occj != 0.0 ) {
                     scr1[i][j]= scr2[i][j];
                     scr1[j][i]= scr2[i][j];
                     }
                  else {
                     scr1[i][j]=scr1[j][i]=0.0;
                     }
                  }
               else if(!twocon) {
                  if ((occi == 1.0 || occi == 0.5) && occj == 2.0 )
                     etemp = scr2[i][j]-0.5*scr3[i][j];
                  else if(occi == 0.0) {
                     if (occj == 2.0) etemp = scr2[i][j];
                     else if (occj != 0.0) etemp = 0.5*scr3[i][j];
                     else etemp = 0.0;
                     }
                  else etemp = 0.0;
                  scr1[i][j] = scr1[j][i] = etemp;
                  }
               else {
                  if (occi != 2.0 && occi != 0.0 && occj == 2.0 )
                     etemp = cim*(c1[m]*scr2[i][j]-c2[m]*scr3[i][j]);
                  else if(occi == 0.0) {
                     if (occj == 2.0) etemp = scr2[i][j];
                     else if (occj != 0.0) etemp = cim*scr3[i][j];
                     else etemp = 0.0;
                     }
                  else etemp = 0.0;
                  scr1[i][j] = scr1[j][i] = etemp;
                  }
               }
            }

      /* transform error matrix into ao basis */
         mxmb(s->cmat,1,nn,scr1,1,nn,scr2,1,nn,nn,nn,nn);
         mxmb(scr2,1,nn,s->cmat,nn,1,d->error[m],1,nn,nn,nn,nn);

         for(i=0; i < nn ; i++) {
            for(j=0; j <= i ; j++) {
               etemp=fabs(scr1[i][j]);
               diiser = MAX0(diiser,etemp);
               }
            }
         }
      }
               
  /* then set up b matrix */

   if (iter > ndiis) {
      for (i=0; i < last ; i++) {
         for (j=0; j <= i ; j++) {
            bold[i][j]=bold[j][i]=bold[i+1][j+1];
            }
         }
      }
   for (i=0; i <= last ; i++) {
      etemp=0.0;
      for (m=0; m < num_ir ; m++) {
         s = &scf_info[m];
         if(nn=s->num_so) {
            sdot(diism[i].error[m],diism[last].error[m],nn,&dotp);
            etemp += dotp;
            }
         }
      bold[i][last]=bold[last][i] = etemp;
      }

   bmat[0][0] = 0.0;
   btemp[0] = -1.0;
   norm = 1.0/bold[0][0];
   for (i=1; i <= last+1 ; i++) {
      bmat[i][0]=bmat[0][i] = -1.0;
      btemp[i] = 0.0;
      for (j=1; j <= i ; j++) {
         bmat[i][j]=bmat[j][i] = bold[i-1][j-1]*norm;
         if(i==j) bmat[i][j] *= scale;
         }
      }

   if (iter-1) {
      flin(bmat,btemp,col,1,&determ);

    /* test for poorly conditioned equations */
      while (fabs(determ) < 1.0e-19 && try < last) {

         try++;
         col--;

         bmat[0][0] = 0.0;
         btemp[0] = -1.0;
         norm=1.0/bold[try][try];
         for (i=1; i <= ndiis-try ; i++) {
            bmat[i][0]=bmat[0][i] = -1.0;
            for (j=1; j <= i ; j++) {
               bmat[i][j]=bmat[j][i]=bold[i+try-1][j+try-1]*norm;
               if(i==j) bmat[i][j] *= scale;
               }
            btemp[i] = 0.0;
            }

         flin(bmat,btemp,col,1,&determ);
         }

      if(fabs(determ) < 10.0e-20) {
          printf(" try %d no good\n",try);
          return;
          }

      if(iter >= it_diis) {
         for (m=0; m < num_ir ; m++) {
            s = &scf_info[m];
            if(nn=s->num_so) {
               for (i=ij=0; i < nn ; i++) {
                  for (j=0; j <= i ; j++,ij++) {
                     int kk=1;
                     etemp=0.0;
                     etempo=0.0;
                     for (k=try; k < last+1 ; k++) {
                        if(iopen) etempo += btemp[kk]*diism[k].fock_o[m][i][j];
                        etemp += btemp[kk]*diism[k].fock_c[m][i][j];
                        kk++;
                        }
                     if(iopen) s->fock_open[ij] = etempo;
                     s->fock_pac[ij] = etemp;
                     }
                  }
               }
            }
         }
      }
   }
