# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/drt
VPATH = /home/usrb/b3/longhorn/psi/src/bin/drt

include ../MakeVars

TRUESRC = bkspac.F bloksd.F constd.F drtcdd.F drtind.F fentry.F getc.F getcnt.F getcod.F getkey.F getlin.F getspe.F getsym.F gprgid.F intgrl.F mndrt.F nextc.F paldus.F prin.F putdrt.F ref.F renumb.F reordr.F scan.F scancj.F skipbl.F squeez.F tapesd.F weight.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: drt.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/drt.1
