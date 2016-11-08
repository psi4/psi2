#ifndef FILE30_H
#define FILE30_H

#define MPOINT 200
#define MCONST 200
#define MCALCS 100
#define MAXL 5    /* Up to g-functions (L = 4) currently supported */

int file30_init(void);
int file30_close(void);
double **file30_rd_blk_scf(int);
double **file30_rd_ccvecs(void);
double *file30_rd_contr(void);
double file30_rd_ecorr(void);
double file30_rd_enuc(void);
double file30_rd_escf(void);
double *file30_rd_evals(void);
double *file30_rd_exps(void);
double **file30_rd_geom(void);
double **file30_rd_scf(void);
double **file30_rd_usotao(void);
double *file30_rd_zvals(void);
char *file30_rd_corr_lab(void);
char **file30_rd_hfsym_labs(void);
char **file30_rd_irr_labs(void);
char *file30_rd_label(void);
char *file30_rd_title(void);
int *file30_rd_clsdpi(void);
int **file30_rd_ict(void);
int file30_rd_iopen(void);
int file30_rd_mxcoef(void);
int file30_rd_nao(void);
int file30_rd_natom(void);
int file30_rd_ncalcs(void);
int file30_rd_nirreps(void);
int file30_rd_nmo(void);
int file30_rd_nprim(void);
int file30_rd_nshell(void);
int file30_rd_nsymhf(void);
int *file30_rd_openpi(void);
int *file30_rd_orbspi(void);
int file30_rd_phase_check(void);
int *file30_rd_sloc(void);
int *file30_rd_smax(void);
int *file30_rd_smin(void);
int *file30_rd_snuc(void);
int *file30_rd_snumg(void);
int *file30_rd_sprim();
int *file30_rd_stype();
void file30_wt_ecorr(double);
void file30_wt_escf(double);
void file30_wt_clsdpi(int *);
void file30_wt_corr_lab(char *);
void file30_wt_geom(double **);
void file30_wt_iopen(int);
void file30_wt_openpi(int *);
void file30_wt_scf(double **);
void file30_wt_evals(double *);
void file30_wt_zvals(double *);
void file30_wt_blk_scf(double **,int);

#endif
