# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/stgrmo
VPATH = /home/usrb/b3/longhorn/psi/src/bin/stgrmo

include ../MakeVars

ifeq ($(target)-f77,mips-dec-ultrix3.0-f77)
  FFLAGS = $(FDBG) $(FOTH)
endif

TRUESRC = alp.F amat.F bet.F blkdta000.F cimul.F davliu.F davlix.F dotab.F \
e0cal.F fentry.F gprgid.F i21d.F pqmat.F regist.F scale.F sort.F stabil.F \
stabix.F stdiag.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: stgrmo.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/stgrmo.1

