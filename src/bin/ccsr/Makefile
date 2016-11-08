# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/ccsr
VPATH = /home/usrb/b3/longhorn/psi/src/bin/ccsr

include ../MakeVars

ifeq ($(ARCH),DEC)
  FOPTIONS = -O1 -static
endif

TRUESRC = \
bktdmp.F consym.F fentry.F gprgid.F read30.F sortb.F sorti.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: ccsr.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/ccsr.1
