# Generated automatically from Makefile.in by configure.
srcdir = /home/usrb/b3/longhorn/psi/src/bin/intsth
VPATH = /home/usrb/b3/longhorn/psi/src/bin/intsth

include ../MakeVars

PSILIBS = io.a alloc.a libparse.a libcoallc.a iomrth.a iomr.a

TRUESRC = aoso2.F aoso2e.F cecob.F fentry.F forms.F genral.F gprgid.F in.F \
          intr.F jandk.F oneint.F s0000.F sf3eq1.F sf3eq2.F sf3eq3.F \
          sf3eq4.F sf3eq5.F sf3eq6.F sf3eq7.F sf3eq8.F sints.F socfpd.F \
          sotoao.F stvcz.F stvint.F syminp.F tints.F togcs.F topitz.F \
          twoint.F vints.F wtint2.F wzero.F xyzint.F

BINOBJ = $(TRUESRC:%.F=%.o)

include ../MakeRules

ifeq ($(target),mips-dec-ultrix3.0)
twoint.o: twoint.f
	$(FC) -O1 -c $< -o $@
endif

install_man:: intsth.1
	$(MKDIRS) $(mandir)/man1
	cat $(top_srcdir)/lib/macro.psi $< | sed "s OS_NAME_HERE UNIX " > \
            $(mandir)/man1/intsth.1

