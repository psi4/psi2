# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/ugeom
VPATH = /home/usrb/b3/longhorn/psi/src/bin/ugeom

include ../MakeVars

ALLOC =

PSILIBS = io.a libparse.a libcoallc.a iomr.a

BINOBJ = ugeom.o gprgid.o rotate.o

include ../MakeRules

install_man:: ugeom.1 geom.2
	$(MKDIRS) $(mandir)/man1
	$(MKDIRS) $(mandir)/man2
	cat $(top_srcdir)/lib/macro.psi $(top_srcdir)/src/bin/ugeom/ugeom.1 \
             | sed "s OS_NAME_HERE UNIX "> $(mandir)/man1/ugeom.1
	cat $(top_srcdir)/lib/macro.psi $(top_srcdir)/src/bin/ugeom/geom.2 \
             | sed "s OS_NAME_HERE UNIX " > $(mandir)/man2/geom.2

