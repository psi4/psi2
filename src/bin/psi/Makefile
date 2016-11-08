# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/psi
VPATH = /home/usrb/b3/longhorn/psi/src/bin/psi

include ../MakeVars

ALLOC =

TRUESRC = execut.F gprgid.F lexec.F psi.F runcmd.c

BINOBJ = $(TRUESRC:%.F=%.o)
BINOBJ := $(BINOBJ:%.c=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
include runcmd.d
endif

install_man:: psi.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/psi.1

