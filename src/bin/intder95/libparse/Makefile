
CODE = libparse.a

F77 = f77
FFLAGS = -O

AR = ar
ARFLAGS = rv
RANLIB = ranlib

PSIPP = ../psipp/psipp

INCLUDES = ./

FSRC = charlib.F error.F fcount.F fio.F frdboo.F frdc.F frdc1.F frdc2.F \
frdc3.F frdc4.F frdca1.F frdca2.F frdcm.F frdct2.F frdf.F frdf1.F \
frdf2.F frdf3.F frdf4.F frdfa1.F frdfa2.F frdfm.F frdft2.F frdi.F \
frdi1.F frdi2.F frdi3.F frdi4.F frdia1.F frdia2.F frdim.F frdit2.F \
fwk.F imat.F initialize.F ioabort.F jabort.F mabort.F pabort.F \
parse.F parsea.F print.F token.F trunit.F

LIBOBJ = $(FSRC:%.F=%.o)

default: $(CODE)

$(CODE): $(LIBOBJ)
	$(AR) $(ARFLAGS) $(CODE) $(LIBOBJ)
	$(RANLIB) $(CODE)

%.o: %.F

%.f: %.F
	$(PSIPP) -I$(INCLUDES) -o $@ $<

%.o: %.f ../psipp/stamp-psipp
	$(F77) $(FFLAGS) -c $<

fio.f: fio.F
	$(PSIPP) -I$(INCLUDES) -o fio.i $<
	sed "s PSI_LIB \"/usr/local/psi/\" " < fio.i > fio.f
	/bin/rm -f fio.i

clean:
	/bin/rm -f $(CODE) *.o core *.i fio.f
