# Generated automatically from Makefile.in by configure.

srcdir = /home/usrb/b3/longhorn/psi/src/lib/libipv1
VPATH = /home/usrb/b3/longhorn/psi/src/lib/libipv1
YACC = bison -y
LEX = flex
LEXLIB = -lfl

include ../MakeVars

CINC += -I$(srcdir)

TARGET = $(top_objdir)/lib/libipv1.a

DEFINES += -DNO_TEMPLATES

TRUESRC = ip_alloc.c ip_cwk.c ip_data.c ip_error.c ip_karray.c ip_print.c \
            ip_read.c scan.l parse.y

CSRC = ip_alloc.c ip_cwk.c ip_data.c ip_error.c ip_karray.c ip_print.c ip_read.c
GENSRC = scan.c y.tab.c
ALLSRC = $(CSRC) $(GENSRC)

GENINC = $(ALLSRC:%.c=%.gbl)

ALLINCLUDE = ip_lib.h ip_error.h ip_types.h ip_global.h $(GENINC)

LIBOBJ = $(ALLSRC:%.c=%.o)

include ../MakeRules

# make sure that the include files get generated
$(LIBOBJ): $(GENINC)

default:: ip_lib.h
	$(MKDIRS) $(objincdir)
	$(INSTALL_INCLUDE) $< $(objincdir)/ip_libv1.h

default:: ip_error.h ip_types.h $(GENINC)
	$(MKDIRS) $(objincdir)/libipv1
	$(INSTALL_INCLUDE) $^ $(objincdir)/libipv1

y.tab.c y.tab.h: parse.y
	$(YACC) -v -d $<

scan.c: scan.l
	$(LEX) -t $^ > scan.c

install_inc:: ip_lib.h
	$(MKDIRS) $(incdir)
	$(INSTALL_INCLUDE) $< $(incdir)/ip_libv1.h

install_inc:: ip_error.h ip_types.h $(GENINC)
	$(MKDIRS) $(incdir)/libipv1
	$(INSTALL_INCLUDE) $^ $(incdir)/libipv1

clean::
	-/bin/rm -f y.tab.* y.output scan.c

depend:: $(ALLINCLUDE) y.tab.h

ifneq ($(DODEPEND),no)
include $(ALLSRC:%.c=%.d)
endif
