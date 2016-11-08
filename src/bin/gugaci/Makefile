# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/gugaci
VPATH = /home/usrb/b3/longhorn/psi/src/bin/gugaci

include ../MakeVars

ifeq ($(ARCH),DEC)
   FOPTIONS = -O2 -automatic
endif

TRUESRC = \
allext.F cimain.F davids.F diagnl.F dot.F fdiag.F fentry.F fold.F foldw.F fourx.F getdrt.F getint.F gprgid.F hi.F hii.F hii2.F hiis.F hiis2.F hij.F hs.F initdg.F initdv.F initpg.F initx2.F loopdt.F loopex.F loopin.F loopy.F loopyd.F lpydti.F normal.F nxtblk.F orthog.F roots.F shape4.F shapes.F square.F squarw.F start.F vecout.F vewpsx.F wfile0.F ww.F ww4x.F wx.F wy.F xw.F xx.F xx4xs.F xy.F yw.F yx.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: gugaci.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/gugaci.1
