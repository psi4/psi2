# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/intdif
VPATH = /home/usrb/b3/longhorn/psi/src/bin/intdif

include ../MakeVars

TRUESRC = \
difa.F difb.F dotx.F fentry.F flin.F gprgid.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: intdif.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/intdif.1

