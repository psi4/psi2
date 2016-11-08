/* $Id: common.h,v 1.9 1998/06/30 14:11:02 sbrown Exp $ */
/* $Log: common.h,v $
/* Revision 1.9  1998/06/30 14:11:02  sbrown
/* *************************************************************
/* *Program Modification                                       *
/* *By: Shawn Brown                                            *
/* *Date: June 30, 1998                                        *
/* *Altered program to make a guess at occupations from the    *
/* *diagonalized core hamiltonian matrix.  Program can now     *
/* *make a guess at the beginning of the calculation or at     *
/* *or at every iteration.  Use the latter at your own risk.   *
/* *See man pages for details on new keywords.                 *
/* *************************************************************
/*
 * Revision 1.8  1997/09/12  13:54:27  crawdad
 * Changing marco name from ULL to PSI_FPTR.
 *
 * Revision 1.7  1997/08/25  21:51:22  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.6  1996/07/03  12:25:44  crawdad
 * Added new phase-checking routines to cscf.  These have been tested more
 * extensively than the last, and are far more intelligently written.
 *
 * Revision 1.5  1995/07/26  20:32:59  sherrill
 * Added support for getting second root of TCSCF and for checking the
 * orbital rotation to make sure orbitals are not swapped.
 *
 * Revision 1.4  1995/07/21  17:37:02  psi
 * Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
 * unidentified changes made after that date were causing problems.
 *
 * Revision 1.1  1991/06/15  20:22:55  seidl
 * Initial revision
 * */

#include <libciomr.h>

#define MAX_BASIS 200
#define MAX_STRING 512

#define SMAT 0
#define TMAT 1
#define VMAT 2

#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

EXTERN FILE *infile, *outfile;

EXTERN double dampsv;           /* scale factor in diis */
EXTERN double repnuc;           /* nuclear repulsion */
EXTERN double etot;             /* electronic and total energies */
EXTERN double lshift;           /* levelshift */
EXTERN double diiser;           /* max off-diag. element in MO fock mat. */
EXTERN double save_ci1,save_ci2; /* ci coefficients for tcscf */
EXTERN double dampd;
EXTERN double dampo;

EXTERN int diisflg;             /* 0 for diis, 1 disables diis */
EXTERN int iconv;               /* dmat convg. criterion */
EXTERN int iopen;               /* 0 for closed, 1 for open, 2 for twocon */
EXTERN int inflg;               /* 0 default, 1 use old guess, 2 use core H */
EXTERN int print;               /* print flag */
EXTERN int fock_typ;            /* 0 for default, 1 for simpler op sh fock m */
EXTERN int ndiis;               /* # of error matrices to keep in diis */
EXTERN int it_diis;             /* iteration to begin diis extrapolation */
EXTERN int itmax;               /* max iterations */

EXTERN int multp;		/* multiplicity of the molecule */
EXTERN int charge;		/* charge of the molecule */
EXTERN int natom;		/* number of atoms in the molecule */
EXTERN int nelec;		/* number of electrons in the molecule */
EXTERN int nbfso;		/* total number of molecular orbitals */
EXTERN int occ_fix;             /* 1 if guess occcupations after*/
                                /*   first iteration,0 if not */
EXTERN int mo_out;              /* 1 if display orbitals in new format at end*/
EXTERN int n_so_typs;           /* number of irreps w/ non-zero num of so's */
EXTERN int nbasis;              /* # basis functions */
EXTERN int nsfmax;              /* max # of so's per irrep */
EXTERN int n_open;              /* # open shells */
EXTERN int num_ir;              /* # of symmetry types */
EXTERN int mxcoef;              /* sum of (# so's per irrep)**2 */
EXTERN int mxcoef2;             /* sum of ioff[# so's per irrep] */
EXTERN int readflg;             /* 1 if using buffered io */
EXTERN int maxbuf;              /* number of integrals per buffer */
EXTERN int num_bufs;            /* number of buffers used */
EXTERN int num_ints;            /* total integrals written to supermatrix */
EXTERN int iter;                /* iteration */
EXTERN int converged;           /* 1 if converged */
EXTERN int hsos;                /* 1 if high spin open shell */
EXTERN int singlet;             /* 1 if open shell singlet */
EXTERN int special;             /* 1 if OPENTYPE=special */
EXTERN int twocon;              /* 1 if tcscf */
EXTERN PSI_FPTR pos34;               /* pointer to location in file34 */
EXTERN int nint;                /* number of pki ints in present batch */
EXTERN int opshl1,opshl2;
EXTERN int opblk1,opblk2;
EXTERN int second_root;         /* get the second root of the MCSCF */
EXTERN int icheck_rot;          /* check orbital rotations? */

EXTERN int itap30,itap34,itap92,itap93;
EXTERN double alpha1,alpha2,alpha3;  /* two configuration things */

EXTERN int ioff[1024];          /* matrix offsets */

EXTERN double *alpha, *beta;    /* arrays with energy coupling coeffs */
EXTERN double *zvals;		/* array for nuclear charges */
EXTERN int *symm_tot;     /* array containing the orbital symmetries in order of energy */
EXTERN double *ener_tot;     /* array containing the orbital energies in order */

EXTERN int *i10;  

EXTERN union psi_buffer {
          int *lbli;
          double *pki;
          double **pki_p;
          } oubuf;

EXTERN struct symm {
  double *smat;
  double *tmat;
  double *hmat;
  double *fock_pac;
  double *fock_open;
  double *fock_eff;
  double *fock_evals;
  double *gmat;
  double *gmato;
  double *pmat;
  double *pmato;
  double *pmat2;
  double *pmato2;
  double *dpmat;
  double *dpmato;
  double **cmat;
  /* STB(4/1/98) - Added array for saving evalues of core H */
  double *hevals;
  /* TDC(6/19/96) - Added array for saving original MO vector */
  double **cmat_orig;
  double **sahalf;
  double *occ_num;
  int nclosed;
  int nopen;
  int nhalf;
  int degeneracy;
  int num_so;
  int os_num;
  int ideg;
  char irrep_label[5];
  } *scf_info;

/* TDC(6/19/96) - Added flag for success or failure of phase checking
   routine */
EXTERN int phase_check;

