# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/bmat
VPATH = /home/usrb/b3/longhorn/psi/src/bin/bmat

include ../MakeVars

TRUESRC = aclear.F arcos.F bottom.F checkd.F cntrct.F convef.F convgd.F \
          dmcopy.F dmdotv.F dmprnt.F dmsymm.F dnegv.F dvaddv.F dvcopy.F \
          dvnegv.F dvskal.F dvtv.F efc.F estime.F expand.F fentry.F fifth.F \
          fndneg.F formd.F formnr.F gdiis.F gprgid.F initef.F machb.F \
          magchk.F nom.F normal.F optefc.F osinv1.F ovrlap.F readdisp.F \
          readsimple.F readsymm.F s2.F scalar.F update.F vektor.F wrtgeo.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: bmat.1 intco.2
	$(MKDIRS) $(mandir)/man1
	$(MKDIRS) $(mandir)/man2
	cat $(top_srcdir)/lib/macro.psi $(top_srcdir)/src/bin/bmat/bmat.1 | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man1/bmat.1
	cat $(top_srcdir)/lib/macro.psi $(top_srcdir)/src/bin/bmat/intco.2 | sed "s OS_NAME_HERE UNIX " > \
           $(mandir)/man2/intco.2
