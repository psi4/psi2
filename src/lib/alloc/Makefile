# Generated automatically from Makefile.in by configure.

srcdir = /home/usrb/b3/longhorn/psi/src/lib/alloc
VPATH = /home/usrb/b3/longhorn/psi/src/lib/alloc

include ../MakeVars

TARGET = $(top_objdir)/lib/alloc.a

CSRC = allocate.c
FSRC = fndcor.F fndcro.F

LIBOBJ = $(CSRC:%.c=%.o) $(FSRC:%.F=%.o)

include ../MakeRules

default:: allocate.o
	$(INSTALL) $< $(top_objdir)/lib

install_target:: allocate.o
	$(INSTALL) $< $(libdir)
