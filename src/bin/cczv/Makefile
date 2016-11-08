# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/cczv
VPATH = /home/usrb/b3/longhorn/psi/src/bin/cczv

include ../MakeVars

TRUESRC = \
abterm.F adde.F aints.F bdents.F bdints.F bints.F cints.F consym.F diisd.F dintss.F ecorr.F eints.F eterm.F fentry.F fints.F frozen.F fterm.F gprgid.F gterm.F hdints.F hfints.F hjterm.F initz.F iterm.F mterm.F ncount.F newz.F nfints.F normz.F offset.F rcints.F rdints.F read30.F read69.F rqdent.F rqdint.F rqeint.F rqfint.F rqterm.F sortt.F symarr.F t2z2.F twrit.F y3.F ydterm.F yeterm.F yfterm.F ysdt1.F zdiff.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: cczv.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/cczv.1
