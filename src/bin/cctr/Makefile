# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/cctr
VPATH = /home/usrb/b3/longhorn/psi/src/bin/cctr

include ../MakeVars

TRUESRC = \
bckchn.F bcktdmp.F consym.F fentry.F gprgid.F htrans.F onei.F pak.F rdone.F rdtwo.F read30.F read34.F reints.F setmbf.F sort.F sort2.F squr.F wrints.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: cctr.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/cctr.1
