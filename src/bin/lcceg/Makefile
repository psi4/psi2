# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/lcceg
VPATH = /home/usrb/b3/longhorn/psi/src/bin/lcceg

include ../MakeVars

TRUESRC = \
addd.F atau.F btau.F ccsdt1.F cdt1.F cisd9.F clus9.F consym.F diisd.F divide.F ecorr.F et2.F fentry.F fint.F fock2.F frozen.F fterm.F goint.F got2.F gprgid.F gvt2.F hiaint.F hoint.F ht1.F hvint.F initt.F itype.F jkt2.F ncount.F pcorr.F rdinpt.F rdints.F read30.F symarr.F tdiff.F trnps.F twrit.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: lcceg.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/lcceg.1

