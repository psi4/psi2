# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/lag
VPATH = /home/usrb/b3/longhorn/psi/src/bin/lag

include ../MakeVars

TRUESRC = \
fentry.F getdrt.F getit.F gordr.F gprgid.F initg.F liordr.F makex.F out1.F out2.F squr.F srtone.F tumult.F vecs.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: lag.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/lag.1

