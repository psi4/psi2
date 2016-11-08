
/* $Id: structs.h,v 2.1 1991/06/15 18:32:29 seidl Exp $ */
/* $Log: structs.h,v $
/* Revision 2.1  1991/06/15 18:32:29  seidl
/* *** empty log message ***
/* */


struct integrals {
     double *smat;
     double *hmat;
     double *twomat;
     double *pmat;
     double **sahalf;
     double **fock_mat;
     double **cmat;
     double **fock_trans;
     double *fock_trns_pac;
     double *fock_eig_vals;
     double **ctrans;
     double *gmat;
     double *fock_pac;
     double *pmat_old;
     } ;

