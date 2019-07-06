/***********************************************************************
*
*                    COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*                    xrfluor.c
*                    Program for quantitative 
*                    X-ray fluorescence analysis
*
*                    M.Bos
*                    University of Twente
*                    Faculty of Chem. Technology
*                    Department of Chemical Analysis
*                    september 1997
*
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"
#include "xrfdefs.h"

int  main(int argc, char *argv[])
{
   multilayer sample ;
   multilayer bulksample ;
   layer *ptop, *pbulk;
   double primint, secint, totint, primint_bulk, totint_bulk;
   double totsecint ;
   struct spectrum_pair *pntrspectrum ;
   double labda_in ;
   int i, j ;
   char el[3];
   char lin[6];
   int z_filt;
   double mth_filt, kev;
   FILE *fpresult;

    
   fpresult=fopen("uitvoer","w");
   init_tables();
   /* debugging
   printf("na init tables \n"); /* */

   set_spectrom();
   /* debugging
   printf("na set spectrom \n"); 
   for (i=0; i < nr_of_meas_lines; i++)
   {
       printf("%d Z: %d lijn %s HT %f mthickn filt %g  \n",
               i, (lines_to_meas[i]).z, (lines_to_meas[i]).line,
                (lines_to_meas[i]).kev, (lines_to_meas[i]).mth_filt);
   } /* */
   def_sample(&sample);
   ptop=get_layer(&sample, 0);
   for (j=0; j < nr_of_meas_lines; j++)
   {
     strcpy(el, get_asymb((lines_to_meas[j]).z, &lijst_els));
     strcpy(lin, (lines_to_meas[j]).line);
  /*   printf("in main lijn %s element Z %d %s el %s\n", lin,
           (lines_to_meas[j]).z,
             get_asymb((lines_to_meas[j]).z, &lijst_els),el); /* */
     z_filt=(lines_to_meas[j]).filt_z;
     mth_filt=(lines_to_meas[j]).mth_filt;
     /* debugging
     printf("Filter Z %d, filter massthickness %g \n",
            z_filt, mth_filt); /* */
     kev=(lines_to_meas[j]).kev;
     pntrspectrum=gen_spec(anode_el, take_off, berryliumd, kev,
                 PNTS_SPECTRUM, z_filt, mth_filt);
     pbulk=mkbulklay(el);
     /* now set up complete sample for bulk element */
      bulksample.nr_layers= 1;
      bulksample.lagen=&pbulk; 
     totint=0.0;
     totsecint= 0.0;
     totint_bulk = 0.0;
     for (i=0; i < PNTS_SPECTRUM; i++)
      {
        labda_in=(pntrspectrum[i]).lambda;
/* printf("Voor pilamba i is %d en intens %g \n",i, (pntrspectrum[i]).intens);  /* debugging */
        primint = pilabda(el,lin,labda_in, ptop)*(pntrspectrum[i]).intens; 
        primint_bulk = pilabda(el,lin,labda_in, pbulk)*(pntrspectrum[i]).intens;
        secint = intrasecflu(el,lin, labda_in , &sample, 0, primint);
    /*    printf("primint %g primbulk %g secint %g\n",
                primint, primint_bulk, secint);  /* */
      /*  printf("secint %g \n",secint); /* */
        totint += primint+secint;
        totsecint += secint ;
        totint_bulk += primint_bulk;
     /*   printf(" totint %g totint_bulk %g \n", totint, totint_bulk); /* */
	 }
      /* debugging 	 */
     printf("totint %s %s is %e bulkint %e secint %e\n",
          el, lin, totint,totint_bulk, totsecint); /* */
     printf("Rel. int %s %s is %e \n",el,lin, totint/totint_bulk);
     fprintf(fpresult,"%s %s %e\n", el, lin, totint/totint_bulk);
   }
   fclose(fpresult);
  return 0;
}
