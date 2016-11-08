
/* $Log: tqli.c,v $
/* Revision 2.2  1991/07/30 04:13:57  seidl
/* do not exit if 30 iterations exceeded
/*
 * Revision 2.1  1991/06/15  18:30:12  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id: tqli.c,v 2.2 1991/07/30 04:13:57 seidl Exp $";

# include "includes.h"

/* diagonalizes tridiagonal matrix output by tred2 */
/* gives only eigenvalues if matz = 0, both eigenvalues and eigenvectors */
/* if matz = 1 */

int tqli(n,d,z,e,matz,toler)
   int n, matz;
   double *d, **z, *e;
   double toler;

   {
      register int k;
      int i,j,l,m,iter;
      double dd,g,r,s,c,p,f,b,h;
      double azi;

      f=0.0;
      if (n == 1) {
         d[0]=z[0][0];
         z[0][0] = 1.0;
         return;
         }

      for (i=1; i < n ; i++) {
         e[i-1] = e[i];
         }
      e[n-1] = 0.0;
      for (l=0; l < n; l++) {
         iter = 0;
L1:
         for (m=l; m < n-1;m++) {
            dd = fabs(d[m]) + fabs(d[m+1]);
#if 0
            if (fabs(e[m])+dd == dd) goto L2;
#else
            if (fabs(e[m]) < toler) goto L2;
#endif
            }
         m=n-1;
L2:
         if (m != l) {
            if (iter++ == 30) {
               fprintf (stderr,"tqli not converging\n");
                continue;
#if 0
               exit(30);
#endif
               }

            g = (d[l+1]-d[l])/(2.0*e[l]);
            r = sqrt(g*g + 1.0);
            g = d[m] - d[l] + e[l]/((g + DSIGN(r,g)));
            s=1.0;
            c=1.0;
            p=0.0;
            for (i=m-1; i >= l; i--) {
               f = s*e[i];
               b = c*e[i];
               if (fabs(f) >= fabs(g)) {
                  c = g/f;
                  r = sqrt(c*c + 1.0);
                  e[i+1] = f*r;
                  s=1.0/r;
                  c *= s;
                  }
               else {
                  s = f/g;
                  r = sqrt(s*s + 1.0);
                  e[i+1] = g*r;
                  c = 1.0/r;
                  s *= c;
                  }
               g = d[i+1] - p;
               r = (d[i]-g)*s + 2.0*c*b;
               p = s*r;
               d[i+1] = g+p;
               g = c*r-b;

               if (matz) {
                  double *zi = z[i];
                  double *zi1 = z[i+1];
                  for (k=n; k ; k--,zi++,zi1++) {
                     azi = *zi;
                     f = *zi1;
                     *zi1 = azi*s + c*f;
                     *zi = azi*c - s*f;
                     }
                  }
               }

            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
            goto L1;
            }
         }
   }
