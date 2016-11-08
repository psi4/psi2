# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/lagzbak
VPATH = /home/usrb/b3/longhorn/psi/src/bin/lagzbak

include ../MakeVars

TRUESRC = \
addz.F back1.F back2.F consym.F cphf.F dsinv.F fentry.F frozen.F gprgid.F indata.F lagr.F nred30.F rdcore.F rdins.F set55.F sort.F sort1.F sort2.F symcnt.F zlag.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: lagzbak.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/lagzbak.1

