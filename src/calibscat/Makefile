CC=gcc

CFLAGS= -O2

OBJS=  attilabda.o selfattilabda.o layerdata.o bldatomd.o relrates.o\
      jmpedg.o omega.o emutable.o characwl.o init_tables.o\
      define_smpl.o set_spectrom.o pilabda.o  lzero.o secfluij.o\
      enhanc.o intrasecflu.o gen_spec.o linesms.o mkbulklay.o\
      compdata.o convcomp.o edges.o nkls.o costerkr.o bscatter.o\
	adapq.o calcsxz.o comptcrossmix.o comptoninteg.o\
	ffstoresm.o fotraylinteg.o raylcrossmix.o raylfotinteg.o\
	raylisesr.o raylispfr.o comptisesr.o

calib: $(OBJS) calscat.c xrfluor.h xrfdefs.h
	$(CC)  $(CFLAGS) -o calscat calscat.c $(OBJS) -lm

install:	
	cp calscat ../../bin/calscat

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


bscatter.o:	xrfluor.h xrfdefs.h bscatter.c
	$(CC)  -c $(CFLAGS) bscatter.c

adapq.o:	adapq.h adapq.c
	$(CC)  -c $(CFLAGS) adapq.c

calcsxz.o:	xrfluor.h xrfdefs.h calcsxz.c
	$(CC)  -c $(CFLAGS) calcsxz.c

comptcrossmix.o:	xrfluor.h xrfdefs.h comptcrossmix.c
	$(CC)  -c $(CFLAGS) comptcrossmix.c

comptisesr.o:	xrfluor.h xrfdefs.h comptisesr.c
	$(CC)  -c $(CFLAGS) comptisesr.c

comptoninteg.o:	xrfluor.h xrfdefs.h comptoninteg.c
	$(CC)  -c $(CFLAGS) comptoninteg.c

ffstoresm.o:	xrfluor.h xrfdefs.h ffstoresm.c
	$(CC)  -c $(CFLAGS) ffstoresm.c

fotraylinteg.o:	xrfluor.h xrfdefs.h fotraylinteg.c
	$(CC)  -c $(CFLAGS) fotraylinteg.c

raylcrossmix.o:	xrfluor.h xrfdefs.h raylcrossmix.c
	$(CC)  -c $(CFLAGS) raylcrossmix.c

raylfotinteg.o:	xrfluor.h xrfdefs.h raylfotinteg.c
	$(CC)  -c $(CFLAGS) raylfotinteg.c


raylisesr.o:	xrfluor.h xrfdefs.h raylisesr.c
	$(CC)  -c $(CFLAGS) raylisesr.c

raylispfr.o:	xrfluor.h xrfdefs.h raylispfr.c
	$(CC)  -c $(CFLAGS) raylispfr.c

