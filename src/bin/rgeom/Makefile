# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/rgeom
VPATH = /home/usrb/b3/longhorn/psi/src/bin/rgeom

include ../MakeVars

BINOBJ = gprgid.o rgeom.o

ALLOC =

include ../MakeRules

install_man:: rgeom.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/rgeom.1

