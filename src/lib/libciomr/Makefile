# Generated automatically from Makefile.in by configure.

srcdir = /home/usrb/b3/longhorn/psi/src/lib/libciomr
VPATH = /home/usrb/b3/longhorn/psi/src/lib/libciomr

include ../MakeVars

CDEF += -DBUFF=0 -DSITEDIR=\"$(datadir)\"

TARGET = $(top_objdir)/lib/libciomr.a

TRUESRC = \
add_arr.c \
add_mat.c \
block_alloc.c \
block_matrix.c \
dot.c \
eigout.c \
eigsort.c \
eivout.c \
errors.c \
ffile.c \
flen.c \
flin.c \
fndcor.c \
get_file_info.c \
get_file_ptr.c \
get_param.c \
i2sec.c \
init_array.c \
init_matrix.c \
init_ptrs.c \
int_array.c \
int_pac.c \
io_getline.c \
io_locate.c \
ioopen.c \
lubksb.c \
ludcmp.c \
mat_to_arr.c \
mmult.c \
mxmb.c \
mxmbol.c \
old_inp.c \
print_array.c \
print_mat.c \
r_async.c \
ram.c \
rclose.c \
resource.c \
rfile.c \
rgetsa.c \
rread.c \
rsetsa.c \
rsp.c \
rwrit.c \
s_async.c \
sec2i.c \
sequential.c \
sq_rsp.c \
sq_to_tri.c \
sread.c \
srew.c \
swrit.c \
tqli.c \
tred2.c \
tri_to_block.c \
tri_to_sq.c \
tstart.c \
wreadw.c \
wwritw.c \
zero.c 

DEPENDINCLUDE = includes.h iomrparam.h common.h pointers.h structs.h types.h

LIBOBJ = $(TRUESRC:%.c=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

# this is needed for compilation of other libs without installation of libciomr
default:: libciomr.h iomrparam.h
	$(MKDIRS) $(objincdir)
	$(INSTALL_INCLUDE) $^ $(objincdir)

install_inc:: libciomr.h iomrparam.h
	$(MKDIRS) $(objincdir)
	$(INSTALL_INCLUDE) $^ $(objincdir)
