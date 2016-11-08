# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/cscf
VPATH = /home/usrb/b3/longhorn/psi/src/bin/cscf

include ../MakeVars

CDEF += -DTIME=0 -DNO_TEMPLATES

LD = $(CC)
LDFLAGS = $(CFLAGS)
LDLIBS = $(CLIBS)

PSILIBS = libfile30.a libciomr.a libipv1.a

TRUESRC = cscf.c cleanup.c diis.c dmat.c dmat_2.c ecalc.c errchk.c findit.c \
          formg2.c formgc.c formgo.c form_vec.c gprgid.c init_scf.c make_47.c \
          packit_c.c packit_o.c rdone.c rdtwo.c rotate_vector.c scf_input.c \
          scf_iter.c scf_iter_2.c schmit.c sdot.c shalf.c check_rot.c phases.c \
          guess.c sortev.c \

BINOBJ = $(TRUESRC:%.c=%.o)
ALLOC = 

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

install_man:: cscf.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/cscf.1

