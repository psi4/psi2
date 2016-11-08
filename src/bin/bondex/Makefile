# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/bondex
VPATH = /home/usrb/b3/longhorn/psi/src/bin/bondex

include ../MakeVars

TRUESRC = \
basset.F blkdta000.F densty.F diporb.F ecal.F eivec.F fentry.F gcal.F gprgid.F ligen3.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: bondex.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/bondex.1
