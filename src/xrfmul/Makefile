CC=gcc

CFLAGS=-O2

OBJS=  attilabda.o selfattilabda.o layerdata.o bldatomd.o relrates.o\
      jmpedg.o omega.o emutable.o characwl.o init_tables.o\
      define_smpl.o set_spectrom.o pilabda.o  lzero.o secfluij.o\
      enhanc.o intrasecflu.o gen_spec.o linesms.o mkbulklay.o\
      compdata.o convcomp.o edges.o nkls.o costerkr.o mupdown.o\
      bigx.o intermediate.o secflufact.o intersecflu.o respfunc.o\
      simplex.o prtmltic.o mulabfil.o



xrfmul:	$(OBJS) xrfsmplx.c respfunc.c xrfluor.h xrfdefs.h
	$(CC) -o xrfmul xrfsmplx.c  $(OBJS) -lm

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

mupdown.o:	xrfluor.h xrfdefs.h mupdown.c
	$(CC) -c $(CFLAGS) mupdown.c

bigx.o: xrfluor.h xrfdefs.h bigx.c
	$(CC) -c $(CFLAGS) bigx.c

intermediate.o:	xrfluor.h xrfdefs.h intermediate.c
	$(CC) -c $(CFLAGS) intermediate.c

secflufact.o:	xrfluor.h xrfdefs.h secflufact.c
	$(CC) -c $(CFLAGS) secflufact.c

intersecflu.o:	xrfluor.h xrfdefs.h intersecflu.c
	$(CC) -c $(CFLAGS) intersecflu.c

respfunc.o:	xrfluor.h xrfdefs.h respfunc.c
	$(CC) -c $(CFLAGS) respfunc.c

simplex.o:	xrfluor.h xrfdefs.h simplex.c
	$(CC) -c $(CFLAGS) simplex.c

prtmltic.o:	xrfluor.h xrfdefs.h prtmltic.c
	$(CC) -c $(CFLAGS) prtmltic.c

mulabfil.o:	xrfluor.h xrfdefs.h mulabfil.c
	$(CC) -c $(CFALGS) mulabfil.c


