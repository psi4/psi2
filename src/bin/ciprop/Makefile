# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/ciprop
VPATH = /home/usrb/b3/longhorn/psi/src/bin/ciprop

include ../MakeVars

TRUESRC = \
basset.F blkdta000.F center.F densty.F dipole.F ecal.F eivec.F fentry.F gcal.F gprgid.F onehgt.F popul.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: ciprop.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/ciprop.1
