# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/dipder
VPATH = /home/usrb/b3/longhorn/psi/src/bin/dipder

include ../MakeVars

TRUESRC = fentry.F gprgid.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: dipder.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/dipder.1
