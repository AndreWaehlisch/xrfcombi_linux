/**********************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*     criss.c
*
*     M.Bos
*     feb 1998
*
**********************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern int nr_of_meas_lines ;
extern multicomplayer compsample ;
extern multilayer sample ;
extern struct spectrum_pair **pntrspectra;
extern meas_line lines_to_meas[MAX_LINES];
extern elem_list lijst_els;
extern double *rxi_meas;
extern double *rxi_calc;
extern double *intprep_calc;
extern double *intbulk_calc;
extern int nr_free_vars;
extern char anode_el[4];
extern double take_off, berryliumd ;
extern double **aij;
extern double **aijj;
extern double ***aijk;



void criss(multilayer *sample)
{
  int i,j;
  layer *ptop ;
  int nr_elem ;
  double delta;
  elem *pelements ;
  char meas_elem[3] ;
  int z_meas_elem ;
  int *seq_elms_in_layer ;
  double newconc, oldconc;
  double ssq, ssqold;
  int no_oxygen=0;    /* default, suppose there is oxygen involved */
  int fix_oxygen ;

  fix_oxygen=1; /*default fixed */
  ptop=(*sample).lagen[0];
  nr_elem= (*ptop).nr_elements ;
  pelements= (*ptop).elementen ;
  seq_elms_in_layer= (int *) malloc(nr_of_meas_lines*sizeof(int));

  /* set up correspondence of elements to measured line */

  for (i=0; i < nr_of_meas_lines ; i++)
  {
     z_meas_elem= (lines_to_meas[i]).z ;
     strcpy(meas_elem, get_asymb(z_meas_elem, &lijst_els));
     j=0;
     while ( strcmp(meas_elem, (pelements[j]).name) !=0 ) j++;
     if (j> nr_elem)
     { 
      printf("error in finding element %s  belonging to line %d\n",
              meas_elem, i);
      exit(1);
     }
     seq_elms_in_layer[i]=j;     
   }
  ssq=1.0e3;
  do
  {
    ssqold=ssq;
    delta =0.0;
   ssq= func_el(sample);
    for (i=0; i < nr_of_meas_lines; i++)
    {
      oldconc=(pelements[seq_elms_in_layer[i]]).wfract ;
      newconc=rxi_meas[i]*oldconc*(1.0 -rxi_calc[i])/
         (rxi_meas[i]*(oldconc-rxi_calc[i])+rxi_calc[i]*(1.0-oldconc));
      delta+= fabs((newconc- oldconc)/newconc);
      (pelements[seq_elms_in_layer[i]]).wfract=newconc;
      /* debugging 
      printf("newconc %d %g \n", i, newconc); /* */
    }
    } while ((delta > 1.0E-3)&& (ssq < ssqold)  )  ; /* */
/*    } while ((delta > 1.0E-3))  ; /* */

       /* now print result  */

/*       printf("Final results weightfractions of elements\n");  */

       printf("element  weightfraction\n");
       printf("-----------------------\n");
       for (j=0; j < nr_of_meas_lines; j++)
       {
         printf("%s       %g\n",pelements[seq_elms_in_layer[j]].name,
                 pelements[seq_elms_in_layer[j]].wfract);
       } 
      printf("Ssq is %g\n", ssq); 
}

