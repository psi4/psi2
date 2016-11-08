# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/mastux
VPATH = /home/usrb/b3/longhorn/psi/src/bin/mastux

include ../MakeVars

ifeq ($(target)-f77,mips-dec-ultrix3.0-f77)
  FFLAGS = $(FDBG) $(FOTH)
endif

TRUESRC = abmat.F fentry.F formms.F gprgid.F moconv.F onei.F opordr.F ordr.F \
rdone.F rdtwo.F rotate.F scfcl.F scfgr.F scfhs.F scfinp.F sort.F zetmat.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: mastux.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/mastux.1

