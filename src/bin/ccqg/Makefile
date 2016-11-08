# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/ccqg
VPATH = /home/usrb/b3/longhorn/psi/src/bin/ccqg

include ../MakeVars

TRUESRC = \
consym.F fentry.F formga.F formgb.F formgc.F formgd.F formge.F formgf.F formq.F frmtau.F frozen.F gprgid.F ncount.F offset.F read30.F symarr.F wrtem.F wrtga.F wrtgc.F wrtgd.F wrtge.F wrtqq.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: ccqg.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/ccqg.1
