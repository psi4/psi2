# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/inputth
VPATH = /home/usrb/b3/longhorn/psi/src/bin/inputth

include ../MakeVars

TRUESRC = abort.F atomd.F atoms.F berror.F blkatm.F bonda.F charges.F \
          cipa.F count.F enuc.F fentry.F find.F geom.F gprgid.F inttyp.F \
          local.F mole.F out.F patoms.F print2.F printbasis.F ptgrp.F \
          readbasis.F rot.F salc.F shmidt.F start.F symbol.F symorb.F \
          trans.F trmat6.F xmain.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: input.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/input.1
