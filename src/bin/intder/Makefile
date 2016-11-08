# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/intder
VPATH = /home/usrb/b3/longhorn/psi/src/bin/intder

include ../MakeVars

TRUESRC = \
ahx2.F ahx3.F ahx4.F binvrt.F bottom.F brow.F fcout.F fentry.F fill3a.F fill3b.F fill4a.F gfmat.F gprgid.F hijks1.F hijks2.F hijks3.F hijks6.F hijs1.F hijs2.F hijs3.F hijs6.F inpfkm.F intder.F lintr.F lprt.F machb.F machb2.F machf2.F machf3.F machf4.F machx.F machy.F machz.F mat1.F mat2.F normco.F numx.F numy.F numz.F readsimple.F readsymm.F scapro.F srtst1.F srtst2.F table1.F table2.F table3.F vcder1.F vcxki.F vecpro.F vect1.F vect2.F vect3.F vect4.F vect5.F vect6.F xf2.F xf3.F xin.F xout.F xrow.F yf2.F yin.F yout.F yrow.F zin.F zout.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

install_man:: intder.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/intder.1

