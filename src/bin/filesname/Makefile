# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/filesname
VPATH = /home/usrb/b3/longhorn/psi/src/bin/filesname

include ../MakeVars

ifeq ($(ARCH),RS6000)
  FFLAGS = $(FFLAGS_NOOPT)
endif

BINOBJ = filesname.o gfname.o
ALLOC =

include ../MakeRules

install:: psiclean
	$(MKDIRS) $(bindir)
	$(INSTALL_PROGRAM) $< $(bindir)

install_man:: psiclean.1 files.2
	$(MKDIRS) $(mandir)/man1
	$(MKDIRS) $(mandir)/man2
	cat $(top_srcdir)/lib/macro.psi $(top_srcdir)/src/bin/filesname/psiclean.1 | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/psiclean.1
	cat $(top_srcdir)/lib/macro.psi $(top_srcdir)/src/bin/filesname/files.2 | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man2/files.2
