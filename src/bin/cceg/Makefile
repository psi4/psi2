# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/cceg
VPATH = /home/usrb/b3/longhorn/psi/src/bin/cceg

include ../MakeVars
FDEF =  -DHAS_IAND -DHAS_IOR -DHAS_IEOR -DHAS_ISHFT -DLinux 
FOPT = -O 

TRUESRC = \
addd.F atau.F btau.F ccsdt1.F cct1.F cdt1.F cisd9.F clus9.F consym.F diisd.F divide.F ecorr.F et2.F fentry.F fock2.F frozen.F ft2.F goint.F got2.F gprgid.F gvint.F gvt2.F hiaint.F hoint.F ht1.F hvint.F initt.F itype.F jkt2.F matout.F ncount.F pcorr.F rdinpt.F rdints.F read30.F sdt.F sdtg.F symarr.F tdiff.F trnps.F twrit.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: cceg.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/cceg.1
