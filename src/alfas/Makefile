CC=gcc

CFLAGS= -O2

OBJS=  attilabda.o selfattilabda.o layerdata.o bldatomd.o relrates.o\
      jmpedg.o omega.o emutable.o characwl.o init_tables.o\
      define_smpl.o set_spectrom.o pilabda.o  lzero.o secfluij.o\
      enhanc.o intrasecflu.o gen_spec.o linesms.o mkbulklay.o\
      compdata.o convcomp.o edges.o nkls.o costerkr.o\
      mulabfil.o  conciter.o\
      rouss.o  mkbinbulk.o  criss.o func_el.o sum_el.o wilabdak.o\
      bscatter.o pij_labda.o delta_ij.o k_j.o rho_ij.o alfa_ij.o\
      copysample.o

OBJSROUS= xrfrous2.o alphrho.o respfuncnv.o  phase2.o fitconc.o\
            prtmltic.o simplex.o

OBJSSMPLX= respfunc.o xrfsmplx.o prtmltic.o simplex.o

OBJSTSTFIT= tstfitcomp.o fitconc.o compdata.o convcomp.o\
            simplex.o bldatomd.o init_tables.o jmpedg.o omega.o\
            emutable.o edges.o costerkr.o characwl.o relrates.o\
            layerdata.o nkls.o prtmltic.o respfuncnv.o copysample.o

OBJSCONVERSION= convertml.o layerdata.o convcomp.o compdata.o bldatomd.o


all:  laguit xrfalfas xrfrous xrfsmplx  massabs conversion

install:
	cp xrfsmplx ../../bin
	cp laguit ../../bin
	cp xrfalfas ../../bin
	cp xrfrous2 ../../bin/xrfrous

xrfsmplx: $(OBJS) $(OBJSSMPLX)  xrfluor.h xrfdefs.h 
	$(CC) $(CFLAGS) -o xrfsmplx  $(OBJS) $(OBJSSMPLX)  -lm


laguit: $(OBJS) laguit.c  xrfluor.h xrfdefs.h
	$(CC) $(CFLAGS) -o laguit laguit.c  $(OBJS) -lm

xrfalfas: $(OBJS) xrfalfas.c  xrfluor.h xrfdefs.h
	$(CC) $(CFLAGS) -o xrfalfas xrfalfas.c  $(OBJS) -lm

xrfcrude: $(OBJS) xrfcrude.c  xrfluor.h xrfdefs.h
	$(CC) $(CFLAGS) -o xrfcrude xrfcrude.c  $(OBJS) -lm

xrfrous: $(OBJS)  $(OBJSROUS) xrfluor.h xrfdefs.h 
	$(CC) $(CFLAGS) -o xrfrous2 $(OBJSROUS) $(OBJS) -lm

massabs:  $(OBJS) massabs.c xrfluor.h xrfdefs.h
	$(CC) $(CFLAGS) -o massabs massabs.c $(OBJS) -lm

tstfitcomp: $(OBJSTSTFIT) xrfluor.h xrfdefs.h
	$(CC) $(CFLAGS) -o tstfitcomp $(OBJSTSTFIT) -lm

conversion: $(OBJSCONVERSION) xrfluor.h xrfdefs.h convertml.c
	$(CC) $(CFLAGS) -o conversion $(OBJSCONVERSION) -lm

attilabda.o: xrfluor.h xrfdefs.h attilabda.c
	$(CC) -c $(CFLAGS) attilabda.c

selfattilabda.o: xrfluor.h xrfdefs.h selfattilabda.c
	$(CC) -c $(CFLAGS) selfattilabda.c

layerdata.o: xrfluor.h xrfdefs.h layerdata.c
	$(CC) -c $(CFLAGS) layerdata.c

bldatomd.o: xrfluor.h xrfdefs.h bldatomd.c
	$(CC) -c $(CFLAGS) bldatomd.c

relrates.o: xrfluor.h xrfdefs.h relrates.c
	$(CC) -c $(CFLAGS) relrates.c

jmpedg.o: xrfluor.h xrfdefs.h jmpedg.c
	$(CC) -c $(CFLAGS) jmpedg.c

emutable.o: xrfluor.h xrfdefs.h emutable.c
	$(CC) -c $(CFLAGS) emutable.c

characwl.o: xrfluor.h xrfdefs.h characwl.c
	$(CC) -c $(CFLAGS) characwl.c

init_tables.o: xrfluor.h xrfdefs.h init_tables.c
	$(CC) -c $(CFLAGS) init_tables.c

define_smpl.o: xrfluor.h xrfdefs.h define_smpl.c
	$(CC) -c $(CFLAGS) define_smpl.c

set_spectrom.o: xrfluor.h xrfdefs.h set_spectrom.c
	$(CC) -c $(CFLAGS) set_spectrom.c

pilabda.o: xrfluor.h xrfdefs.h pilabda.c
	$(CC) -c $(CFLAGS) pilabda.c

lzero.o: xrfluor.h xrfdefs.h lzero.c
	$(CC) -c $(CFLAGS) lzero.c

secfluij.o: xrfluor.h xrfdefs.h secfluij.c
	$(CC) -c $(CFLAGS) secfluij.c

enhanc.o:	xrfluor.h xrfdefs.h enhanc.c
	$(CC) -c $(CFLAGS) enhanc.c

intrasecflu.o: 	xrfluor.h xrfdefs.h intrasecflu.c
	$(CC) -c $(CFLAGS) intrasecflu.c

gen_spec.o:	xrfluor.h xrfdefs.h gen_spec.c
	$(CC) -c $(CFLAGS) gen_spec.c

linesms.o:	xrfluor.h xrfdefs.h linesms.c
	$(CC) -c $(CFLAGS) linesms.c

mkbulklay.o:	xrfluor.h xrfdefs.h mkbulklay.c
	$(CC) -c $(CFLAGS) mkbulklay.c 

compdata.o:	xrfluor.h xrfdefs.h compdata.c
	$(CC) -c $(CFLAGS) compdata.c

convcomp.o:	xrfluor.h xrfdefs.h convcomp.c
	$(CC) -c $(CFLAGS) convcomp.c

edges.o:	xrfluor.h xrfdefs.h edges.c
	$(CC) -c $(CFLAGS) edges.c

nkls.o:		xrfluor.h xrfdefs.h nkls.c
	$(CC)  -c $(CFLAGS) nkls.c

costerkr.o:	xrfluor.h xrfdefs.h costerkr.c
	$(CC) -c $(CFLAGS) costerkr.c

omega.o:	xrfluor.h xrfdefs.h omega.c
	$(CC) -c $(CFLAGS) omega.c

respfunc.o:	xrfluor.h xrfdefs.h respfunc.c
	$(CC) -c $(CFLAGS) respfunc.c

simplex.o:	xrfluor.h xrfdefs.h simplex.c
	$(CC) -c $(CFLAGS) simplex.c

mulabfil.o:	xrfluor.h xrfdefs.h mulabfil.c
	$(CC) -c $(CFLAGS) mulabfil.c

conciter.o:	xrfluor.h conciter.c
	$(CC) -c $(CFLAGS) conciter.c

newton.o:	xrfluor.h newton.c
	$(CC) -c $(CFLAGS) newton.c

copysample.o:	xrfluor.h copysample.c
	$(CC) -c $(CFLAGS) copysample.c

criss.o:	xrfluor.h criss.c
	$(CC) -c $(CFLAGS) criss.c

func_el.o:	xrfluor.h func_el.c
	$(CC) -c $(CFLAGS) func_el.c

sum_el.o:	xrfluor.h  sum_el.c
	$(CC)  -c  $(CFLAGS) sum_el.c

wilabdak.o:	xrfluor.h wilabdak.c
	$(CC)  -c $(CFLAGS) wilabdak.c

alfa_ij.o:	xrfluor.h alfa_ij.c
	$(CC)  -c $(CFLAGS) alfa_ij.c

rho_ij.o:	xrfluor.h rho_ij.c
	$(CC) -c $(CFLAGS) rho_ij.c

k_j.o:	xrfluor.h k_j.c
	$(CC) -c $(CFLAGS) k_j.c

delta_ij.o:	xrfluor.h delta_ij.c
	$(CC) -c $(CFLAGS) delta_ij.c

pij_labda.o:	xrfluor.h pij_labda.c
	$(CC) -c $(CFLAGS) pij_labda.c

alphrho.o:	xrfluor.h alphrho.c
	$(CC) -c $(CFLAGS) alphrho.c

phase2.o:	xrfluor.h phase2.c
	$(CC) -c $(CFLAGS) phase2.c

bscatter.o:	xrfluor.h bscatter.c
	$(CC) -c $(CFLAGS) bscatter.c

respfuncnv.o:	xrfluor.h respfuncnv.c
	$(CC) -c $(CFLAGS) respfuncnv.c

fitconc.o:	xrfluor.h fitconc.c
	$(CC) -c $(CFLAGS) fitconc.c
