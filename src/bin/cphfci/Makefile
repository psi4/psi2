# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/cphfci
VPATH = /home/usrb/b3/longhorn/psi/src/bin/cphfci

include ../MakeVars

TRUESRC = \
aborts.F alloc.F b0mat.F cidip.F dermat.F dot.F drmat2.F famat.F famat2.F fentry.F flinx.F geta.F getb.F getcv.F getdrt.F getint.F gprgid.F gtint2.F gtints.F hamat.F mainin.F moconv.F order.F order1.F order2.F prn11.F read40.F regist.F samat.F squr.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: cphfci.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/cphfci.1
