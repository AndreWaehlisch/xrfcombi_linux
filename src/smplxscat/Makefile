CC=gcc

CFLAGS= -O3 -march=pentium4 -msse2 -ffast-math

OBJS=  attilabda.o selfattilabda.o layerdata.o bldatomd.o relrates.o\
      jmpedg.o omega.o emutable.o characwl.o init_tables.o\
      define_smpl.o set_spectrom.o pilabda.o  lzero.o secfluij.o\
      enhanc.o intrasecflu.o gen_spec.o linesms.o mkbulklay.o\
      compdata.o convcomp.o edges.o nkls.o costerkr.o\
      comptcrossmix.o comptoninteg.o comptisesr.o calcsxz.o\
      adapq.o raylcrossmix.o ffstoresm.o raylisesr.o raylfotinteg.o\
      raylispfr.o fotraylinteg.o defsmplscat.o pilabdascat.o\
      secfluijscat.o intrasecfluscat.o convcompscat.o mulabfil.o\
      prtmltic.o respfunc.o simplex.o xrfsim.o fitconc.o respfuncnv.o\
      copysample.o respfuncscat.o


all: $(OBJS) xrfsmplx.c xrfluor.h
	$(CC) $(CFLAGS) -o smplxscat xrfsmplx.c  $(OBJS) -lm

install:
	cp smplxscat ../../bin


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

comptcrossmix.o:	xrfluor.h xrfdefs.h comptcrossmix.c
	$(CC) -c $(CFLAGS) comptcrossmix.c

comptoninteg.o:	xrfluor.h xrfdefs.h comptoninteg.c
	$(CC) -c $(CFLAGS) comptoninteg.c

comptisesr.o:	xrfluor.h xrfdefs.h comptisesr.c
	$(CC) -c $(CFLAGS) comptisesr.c

calcsxz.o:	xrfluor.h xrfdefs.h calcsxz.c
	$(CC) -c $(CFLAGS) calcsxz.c

adapq.o:	adapq.h adapq.c
	$(CC) -c $(CFLAGS) adapq.c

raylcrossmix.o:	xrfluor.h xrfdefs.h raylcrossmix.c
	$(CC) -c $(CFLAGS) raylcrossmix.c

ffstoresm.o:	xrfluor.h xrfdefs.h ffstoresm.c
	$(CC) -c $(CFLAGS) ffstoresm.c

raylisesr.o:	xrfluor.h xrfdefs.h raylisesr.c
	$(CC) -c $(CFLAGS) raylisesr.c

raylfotinteg.o:	xrfluor.h xrfdefs.h raylfotinteg.c
	$(CC) -c $(CFLAGS) raylfotinteg.c

raylispfr.o:	xrfluor.h xrfdefs.h raylispfr.c
	$(CC) -c $(CFLAGS) raylispfr.c

fotraylinteg.o:	xrfluor.h xrfdefs.h fotraylinteg.c
	$(CC) -c $(CFLAGS) fotraylinteg.c

defsmplscat.o:	xrfluor.h xrfdefs.h defsmplscat.c
	$(CC) -c $(CFLAGS) defsmplscat.c

pilabdascat.o:	xrfluor.h xrfdefs.h pilabdascat.c
	$(CC) -c $(CFLAGS) pilabdascat.c

secfluijscat.o:	xrfluor.h xrfdefs.h secfluijscat.c
	$(CC) -c $(CFLAGS) secfluijscat.c

intrasecfluscat.o:	xrfluor.h xrfdefs.h intrasecfluscat.c
	$(CC) -c $(CFLAGS) intrasecfluscat.c

convcompscat.o:	xrfluor.h xrfdefs.h convcompscat.c
	$(CC) -c $(CFLAGS) convcompscat.c

mulabfil.o:	xrfluor.h xrfdefs.h mulabfil.c
	$(CC) -c $(CFLAGS) mulabfil.c

prtmltic.o:	xrfluor.h xrfdefs.h prtmltic.c
	$(CC) -c $(CFLAGS) prtmltic.c

respfunc.o:	xrfluor.h xrfdefs.h respfunc.c
	$(CC) -c $(CFLAGS) respfunc.c

respfuncscat.o:	xrfluor.h xrfdefs.h respfunc.c
	$(CC) -c $(CFLAGS) respfuncscat.c

simplex.o:	xrfluor.h xrfdefs.h simplex.c
	$(CC) -c $(CFLAGS) simplex.c


xrfsim.o:	xrfluor.h xrfdefs.h xrfsim.c
	$(CC) -c $(CFLAGS) xrfsim.c

fitconc.o:	xrfluor.h xrfdefs.h fitconc.c
	$(CC) -c $(CFLAGS) fitconc.c

respfuncnv.o:	xrfluor.h xrfdefs.h respfuncnv.c
	$(CC) -c $(CFLAGS) respfuncnv.c

copysample.o:	xrfluor.h xrfdefs.h copysample.c
	$(CC) -c $(CFLAGS) copysample.c

