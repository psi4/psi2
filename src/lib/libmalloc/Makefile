# Generated automatically from Makefile.in by configure.
#
#	SQ sccsid @(#)Makefile	1.4 88/08/25
#	NOT copyright by SoftQuad.  - msb, 1988
#
# This is configured to make the traced version of malloc; see mallck(1)
# and ./TRACE_README.  To make the regular version, delete -DMALLOCTRACE.
#							Mark Brader
#

srcdir = /home/usrb/b3/longhorn/psi/src/lib/libmalloc
VPATH = /home/usrb/b3/longhorn/psi/src/lib/libmalloc

include ../MakeVars

TARGET = $(top_objdir)/lib/libmalloc.a

CDEF += -Ddebug

MALLOC_C = malloc.c free.c realloc.c forget.c write_trace.c init_trace.c
LIBOBJ = $(MALLOC_C:%.c=%.o)

include ../MakeRules

depend:: malloc.h

mtests:	mtest1 mtest2 mtest3

mtest1:	mtest1.c

mtest2:	mtest2.c

mtest3:	mtest3.c

mtest1 mtest2 mtest3:	malloctrace.a
	${CC} $(CARGS) $@.c malloctrace.a -o $@

clean::
	-/bin/rm -f tstmalloc

$(LIBOBJ:%.o=%.d): $(DEPENDINCLUDE)

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif
