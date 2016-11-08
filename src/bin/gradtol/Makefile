# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/gradtol
VPATH = /home/usrb/b3/longhorn/psi/src/bin/gradtol

include ../MakeVars

TRUESRC = \
gprgid.F gradtol.F

BINOBJ = $(TRUESRC:%.F=%.o)
ALLOC =

include ../MakeRules

install_man:: gradtol.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
          $(mandir)/man1/gradtol.1

