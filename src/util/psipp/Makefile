# Generated automatically from Makefile.in by configure.

srcdir = /home/usrb/b3/longhorn/psi/src/util/psipp
top_srcdir = /home/usrb/b3/longhorn/psi
top_objdir = /home/usrb/b3/longhorn/psi
MKDIRS = $(top_srcdir)/mkdirs.sh
VPATH = /home/usrb/b3/longhorn/psi/src/util/psipp

prefix = /usr/local/psi
incdir = $(prefix)/include
target = i686-generic-linux2.4
objincdir = /home/usrb/b3/longhorn/psi/include
exec_prefix = $(prefix)/$(target)
bindir = ${exec_prefix}/bin

YACC = bison -y
INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL}

COPT = -O
CINC = -I. -I$(objincdir) -I$(incdir)
CDBG = 
COTH = 
CDEF =  -DMAIN_FUNC=MAIN__ -DFCLINK=1 -DLinux
CFLAGS = $(CDBG) $(COPT) $(CDEF) $(COTH) $(CINC)

CC = gcc

DODEPEND = no

OBJ = \
error.o \
if.o \
include.o \
main.o \
name.o \
operator.o \
parse.o \
scan.o \
print.o

IOBJ = \
ipsipp.o

.PHONY:	all default install intall_target install_inc clean oclean dclean targetclean

all: default
default: psipp ipsipp

psipp: $(OBJ)
	$(CC) $(CFLAGS) -o $@  $(OBJ)

ipsipp: $(IOBJ)
	$(CC) $(CFLAGS) -o ipsipp $(IOBJ)

parse.c: parse.y
	$(YACC) -v -d $<
	mv y.tab.c parse.c
	mv y.tab.h parse.h

parse.h: parse.c

install: default
	$(MKDIRS) $(bindir) 
	$(INSTALL_PROGRAM) psipp $(bindir)
	$(INSTALL_PROGRAM) ipsipp $(bindir)

install_inc:

clean:
	rm -f *.o

targetclean:
	rm -f psipp ipsipp

ifneq ($(DODEPEND),no)
$(OBJ:.o=.d): parse.h
include $(OBJ:.o=.d)
else
$(OBJ): parse.h
endif

