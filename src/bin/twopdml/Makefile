# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/twopdml
VPATH = /home/usrb/b3/longhorn/psi/src/bin/twopdml

include ../MakeVars

ifeq ($(target)-$(FC),mips-dec-ultrix3.0-f77)
  FFLAGS = -O1 -static -Nl2000
endif

TRUESRC = \
dump.F fentry.F fold.F foldw.F getdrt.F gprgid.F ieorx1.F ieorx2.F ieorx3.F init4x.F initdg.F initex.F initlp.F initpg.F initx2.F initxx.F loopyd.F lpdata.F main.F square.F squarw.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: twopdml.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/twopdml.1

