# Generated automatically from Makefile.in by configure.

.PHONY:	default all install depend clean targetclean

default: all

all:

install_inc:
depend:
clean:
targetclean:

top_srcdir = /home/usrb/b3/longhorn/psi
prefix = /usr/local/psi
srcdir = /home/usrb/b3/longhorn/psi/lib
datadir = ${prefix}/share
VPATH = /home/usrb/b3/longhorn/psi/lib
top_objdir = /home/usrb/b3/longhorn/psi
MKDIRS = $(top_srcdir)/mkdirs.sh
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644

install: pbasis.dat psi.dat tmpdisks.dat
	$(MKDIRS) $(datadir)
	$(INSTALL_DATA) $(srcdir)/pbasis.dat $(datadir)
	$(INSTALL_DATA) $(srcdir)/psi.dat $(datadir)
	$(INSTALL_DATA) $(srcdir)/tmpdisks.dat $(datadir)

install_man: macro.psi
	$(MKDIRS) $(datadir)
	$(INSTALL_DATA) $< $(datadir)

$(top_srcdir)/configure: $(top_srcdir)/configure.in $(top_srcdir)/aclocal.m4
	cd $(top_srcdir) && autoconf

$(top_objdir)/config.status: $(top_srcdir)/configure
	cd $(top_objdir) && ./config.status --recheck

Makefile: $(srcdir)/Makefile.in $(top_objdir)/config.status
	cd $(top_objdir) && CONFIG_FILES=lib/Makefile ./config.status
