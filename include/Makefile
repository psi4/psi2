# Generated automatically from Makefile.in by configure.

.PHONY:	default all install depend clean targetclean

default: all

all:

install_man:
depend:
clean:
targetclean:

top_srcdir = /home/usrb/b3/longhorn/psi
prefix = /usr/local/psi
srcdir = /home/usrb/b3/longhorn/psi/include
incdir = $(prefix)/include
VPATH = /home/usrb/b3/longhorn/psi/include
top_objdir = /home/usrb/b3/longhorn/psi
MKDIRS = $(top_srcdir)/mkdirs.sh
INSTALL = /usr/bin/install -c
INSTALL_INCLUDE = $(top_srcdir)/installinc.sh

install: install_inc

install_inc: cfuncs.h ffuncs.h psiconfig.h sizeof.h
	$(MKDIRS) $(incdir)
	$(INSTALL_INCLUDE) $^ $(incdir)

$(top_srcdir)/configure: $(top_srcdir)/configure.in $(top_srcdir)/aclocal.m4
	cd $(top_srcdir) && autoconf

$(top_objdir)/config.status: $(top_srcdir)/configure
	cd $(top_objdir) && ./config.status --recheck

Makefile: $(srcdir)/Makefile.in $(top_objdir)/config.status
	cd $(top_objdir) && CONFIG_FILES=include/Makefile ./config.status
