# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/opdm_so
VPATH = /home/usrb/b3/longhorn/psi/src/bin/opdms_so

include ../MakeVars

ifeq ($(target)-f77,mips-dec-ultrix3.0-f77)
  FFLAGS = $(FDBG) $(FOTH)
endif

TRUESRC = fentry.F getdrt.F getit.F gprgid.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: opdm_so.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/opdm_so.1

