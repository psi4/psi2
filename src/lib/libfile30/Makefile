# Generated automatically from Makefile.in by configure.

srcdir = /home/usrb/b3/longhorn/psi/src/lib/libfile30
VPATH = /home/usrb/b3/longhorn/psi/src/lib/libfile30

include ../MakeVars

TARGET = $(top_objdir)/lib/libfile30.a

TRUESRC = \
init.c \
close.c \
rd_ccvecs.c \
rd_snumg.c \
rd_stype.c \
rd_sloc.c \
rd_sprim.c \
rd_smax.c \
rd_smin.c \
rd_snuc.c \
rd_exps.c \
rd_nshell.c \
rd_nprim.c \
rd_contr.c \
rd_geom.c \
rd_nirreps.c \
rd_nsymhf.c \
rd_nmo.c \
rd_nao.c \
rd_mxcoef.c \
rd_iopen.c \
rd_label.c \
rd_title.c \
rd_phase_check.c \
rd_natom.c \
rd_evals.c \
rd_scf.c \
rd_blk_scf.c \
rd_escf.c \
rd_enuc.c \
rd_ecorr.c \
rd_irr_labs.c \
rd_hfsym_labs.c \
rd_corr_lab.c \
rd_orbspi.c \
rd_clsdpi.c \
rd_openpi.c \
rd_usotao.c \
rd_zvals.c \
rd_ict.c \
wt_escf.c \
wt_ecorr.c \
wt_clsdpi.c \
wt_corr_lab.c \
wt_geom.c \
wt_iopen.c \
wt_openpi.c \
wt_scf.c \
wt_evals.c \
wt_blk_scf.c \
wt_zvals.c \

DEPENDINCLUDE = file30.h

LIBOBJ = $(TRUESRC:%.c=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(incdir)
	$(INSTALL_INCLUDE) $^ $(incdir)

