# Generated automatically from Makefile.in by configure.

srcdir = /home/usrb/b3/longhorn/psi/src/lib/iomr
VPATH = /home/usrb/b3/longhorn/psi/src/lib/iomr

include ../MakeVars

TARGET = $(top_objdir)/lib/iomr.a

FSRC = \
ambc.F ambct.F ambtc.F ambtct.F apbc.F apbct.F apbctt.F apbtc.F apbtct.F \
atebc.F atembc.F atmbc.F atpbc.F atpbct.F bioflu.F bit.F copy.F dfunc.F \
digit.F dmpdat.F dmprfl.F dmpunt.F dnode.F dot.F dott.F droot.F drum.F \
dsmit.F ebc.F ebct.F ebtc.F ebtct.F eigout.F eivout.F embc.F embct.F embtc.F \
embtct.F ffile.F fioflu.F flin.F flinq.F frqout.F givens.F i2sec.F iadtwp.F \
init47.F initmf.F intowp.F iobuf.F ioerr.F iround.F isoout.F itrunc.F izero.F \
lnkerr.F locate.F mabort.F matinv.F matout.F mread.F mwrit.F mxmb_br.F \
mxmmips_c_cd.F mxmtnc.F mxmtnd.F mxmtrn.F mxmttc.F mxmttd.F ncdlbl.F nrec.F \
pack.F pack32.F pack64.F pak.F print.F rclose.F rfile.F rgetsa.F root4.F \
root5.F rread.F rsetsa.F rsp.F rt123.F rwrit.F sec2i.F sread.F srew.F swrit.F \
tql2.F tqlrat.F trbak3.F tred3.F tstart.F vecout.F wfile.F wind.F \
wpadti.F wreadw.F wt47.F wwritw.F zero.F zfile.F btest.F

CSRC = init_matrix.c min_mmult.c mmult.c ftime.c tstartc.c

LIBOBJ = $(FSRC:%.F=%.o) $(CSRC:%.c=%.o)

include ../MakeRules

# These routines can be called by different threads at the same time
# so local data must be on the stack.
zero.$(OBJSUF): zero.$(F77SUF)
	$(FC) $(FFLAGS_NOSTATIC) -c $< $(OUTPUT_OPTION)
bit.$(OBJSUF): bit.$(F77SUF)
	$(FC) $(FFLAGS_NOSTATIC) -c $< $(OUTPUT_OPTION)

ifneq (,$(findstring dec-ultrix,$(target)))
ifeq ($(FC),f77)
tstart.o: tstart.f
	$(FC) -g -static -c $< -o $@
endif
endif

ifneq (,$(findstring power2-ibm-aix,$(target)))
ifeq ($(FC),xlf)
rt123.o: rt123.f
	$(FC) -qxflag=xalias -qdpc -qarch=pwrx -qtune=pwrx -c $< -o $@

root4.o: root4.f
	$(FC) -qxflag=xalias -qdpc -qarch=pwrx -qtune=pwrx -c $< -o $@

root5.o: root5.f
	$(FC) -qxflag=xalias -qdpc -qarch=pwrx -qtune=pwrx -c $< -o $@
endif
endif

ifneq (,$(findstring linux,$(target)))
givens.o: givens.f
	$(FC) $(FDBG) $(FDEF) $(FOTH) $(FINC) -c $< $(OUTPUT_OPTION)
endif

ifneq ($(DODEPEND),no)
include $(CSRC:%.c=%.d)
endif
