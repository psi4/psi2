# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/cphf
VPATH = /home/usrb/b3/longhorn/psi/src/bin/cphf

include ../MakeVars

CDEF += -DTIME=0

PSILIBS = libciomr.a libipv1.a

LD = $(CC)
LDFLAGS = $(CFLAGS)
LDLIBS = $(CLIBS)

ifneq (,$(findstring ibm-aix,$(target)))
  LIBS += libmalloc.a
  CC = xlc
endif

COMMON_SRC = ao_to_mo.c dipole_derivs.c init_master.c mwrit.c

CPHF_SRC = bafmat.c bamat.c bamat_o.c bfmat.c cphf.c cphf_iter.c densmatd.c \
           densmatd_o.c derivs.c dipmo.c famat.c famat_o.c gprgid.c hamat.c \
           hemat.c linear.c makeab.c makeab_o.c make_vec.c orthog.c polar.c \
           scf2nd.c scf2nd_o.c uxmat.c wamat.c wamat_o.c

CDIP_SRC = cdip.c link.c gprgid_cdip.c

ALL_SRC = $(COMMON_SRC) $(CPHF_SRC) $(CDIP_SRC)

BINOBJ = $(CPHF_SRC:%.c=%.o) $(COMMON_SRC:%.c=%.o)

CDIPOBJ = $(CDIP_SRC:%.c=%.o) $(COMMON_SRC:%.c=%.o)

ALLOC = # main is provided so get rid of allocate.o
include ../MakeRules

depend:: includes.h common.h dipder.h

cdip: $(CDIPOBJ) $(LIBS) $(SYSLIBS)
	$(CC) $(LDFLAGS) $^ -o $@

ifneq (,$(findstring linux,$(target)))
cphf_iter.o: cphf_iter.c
	$(CC) $(CDBG) $(CDEF) $(COTH) $(CINC) -c $< $(OUTPUT_OPTION)
endif


ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(ALL_SRC:%.c=%.d)
endif

install_man:: cphf.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/cphf.1
