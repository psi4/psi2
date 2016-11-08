# Generated automatically from Makefile.in by configure.

srcdir = /home/usrb/b3/longhorn/psi/src/lib/libcoallc
VPATH = /home/usrb/b3/longhorn/psi/src/lib/libcoallc

include ../MakeVars

TARGET = $(top_objdir)/lib/libcoallc.a

TRUESRC = \
  coallc.F cochec.F codump.F cofree.F coinit.F corllc.F cosize.F sizeof.F

LIBOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: memory_allocation.5
	$(MKDIRS) $(mandir)/man5
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man5/memory_allocation.5

