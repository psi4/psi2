# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/makepk
VPATH = /home/usrb/b3/longhorn/psi/src/bin/makepk

include ../MakeVars

CDEF += -DTIME=0

PSILIBS = libipv1.a libciomr.a


ifneq (,$(findstring ibm-aix,$(target)))
  LIBS += libmalloc.a
endif

LD = $(CC)
LDFLAGS = $(CFLAGS)
LDLIBS = $(CLIBS)

TRUESRC = abmat.c clscf.c cmake37.c density.c findit.c gprgid.c init_master.c \
          make_fock.c make_zeta.c mwrit.c packit.c rdone.c rdtwo.c \
          scf_stuff.c vec_sort.c zetmat.c

BINOBJ = $(TRUESRC:%.c=%.o)

ALLOC =
include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

install_man:: makepk.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/makepk.1

