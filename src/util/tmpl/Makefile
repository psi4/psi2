# Generated automatically from Makefile.in by configure.

srcdir = /home/usrb/b3/longhorn/psi/src/util/tmpl
top_srcdir = /home/usrb/b3/longhorn/psi
top_objdir = /home/usrb/b3/longhorn/psi
VPATH = /home/usrb/b3/longhorn/psi/src/util/tmpl

prefix = /usr/local/psi
incdir = $(prefix)/include
objincdir = /home/usrb/b3/longhorn/psi/include
target = i686-generic-linux2.4
exec_prefix = $(prefix)/$(target)
bindir = ${exec_prefix}/bin

YACC = bison -y
LEX = flex
LEXLIB = -lfl
INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL}
INSTALL_INCLUDE = $(top_srcdir)/installinc.sh

COPT = -O
CINC = -I. -I$(objincdir) -I$(incdir)
CDBG = 
COTH = 
CDEF =  -DMAIN_FUNC=MAIN__ -DFCLINK=1 -DLinux
CFLAGS = $(CDBG) $(COPT) $(CDEF) $(COTH) $(CINC)

CC = gcc

DODEPEND = no

MKDIRS = $(top_srcdir)/mkdirs.sh

BINOBJ = main.o tmpllex.o

.PHONY:	all default install intall_target install_inc clean oclean dclean targetclean

all: default
default: tmpl $(objincdir)/tmpl.h

tmpl: $(BINOBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LEXLIB)

$(objincdir)/tmpl.h: tmpl.h
	$(MKDIRS) $(objincdir)
	$(INSTALL_INCLUDE) $< $(objincdir)

install: install_target install_inc

install_target: tmpl
	$(MKDIRS) $(bindir)
	$(INSTALL_PROGRAM) tmpl $(bindir)

install_inc: tmpl.h
	$(MKDIRS) $(incdir)
	$(INSTALL_INCLUDE) $< $(incdir)

tmpllex.c: tmpllex.l
	$(LEX) -t $^ > tmpllex.c

tmplparse.c: tmplparse.y
	$(YACC) -d -v $^
	mv -f y.tab.c tmplparse.c
	mv -f y.tab.h tmplparse.h

clean: oclean dclean targetclean
	rm -f tmplparse.c tmplparse.h y.tab.c y.tab.h

oclean:
	rm -f *.o

dclean:
	rm -f *.d

targetclean:
	rm -f tmpl

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

