# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/contour3d
VPATH = /home/usrb/b3/longhorn/psi/src/bin/contour3d

include ../MakeVars

#
# there is no memory allocation for this program, linking allocate.o
#   messes up the DEC3100 compilation
#
ALLOC =

TRUESRC = \
basset.F contour.F drawps.F fioflu2.F gprgid.F intplt.F psi.F spline.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

clean::
	/bin/rm -f *.f

install_man:: contour3d.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/contour3d.1
