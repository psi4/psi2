# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/extr3D_mp0703
VPATH = /home/usrb/b3/longhorn/psi/src/bin/extr3D_mp0703

include ../MakeVars

ifeq ($(target)-f77,mips-dec-ultrix3.0-f77)
  FFLAGS = $(FDBG) $(FOTH)
endif

TRUESRC = gprgid.F extr3D_mp0703.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: extr3D_mp0703.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/extr3D_mp0703.1

