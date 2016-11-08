# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/backtr
VPATH = /home/usrb/b3/longhorn/psi/src/bin/backtr

include ../MakeVars

ifeq ($(ARCH),DEC)
  FOPTIONS = -O
endif

TRUESRC = \
back1.F back2.F ctsqc.F dime.F eivt.F fentry.F gprgid.F init55.F read30.F setup1.F sort2.F squr.F trandm.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: backtr.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
          $(mandir)/man1/backtr.1
