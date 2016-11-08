
CODE = intder_io.a

F77 = f77
FFLAGS = -O

AR = ar
ARFLAGS = rv
RANLIB = ranlib

PSIPP = ../psipp/psipp

INCLUDES = ./

FSRC = \
bioflu.F digit.F flin.F intowp.F ioerr.F locate.F mabort.F rclose.F \
rfile.F rsp.F srew.F tql2.F tqlrat.F trbak3.F tred3.F tstart.F wreadw.F \
wwritw.F

LIBOBJ = $(FSRC:%.F=%.o)

default: ../psipp/psipp $(CODE)

$(CODE): $(LIBOBJ)
	$(AR) $(ARFLAGS) $(CODE) $(LIBOBJ)
	$(RANLIB) $(CODE)

%.o: %.F

%.f: %.F
	$(PSIPP) -I$(INCLUDES) -o $@ $<

%.o: %.f ../psipp/stamp-psipp
	$(F77) $(FFLAGS) -c $<

clean:
	/bin/rm -f $(CODE) *.o core

