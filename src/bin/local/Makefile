# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/local
VPATH = /home/usrb/b3/longhorn/psi/src/bin/local

include ../MakeVars

#
# there is no memory allocation for this program, linking allocate.o
#   messes up the DEC3100 compilation
#
ALLOC =

TRUESRC = \
basset.F cmoao.F derao.F ecal.F gprgid.F locapm.F pm.F readmo.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: local.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/local.1

