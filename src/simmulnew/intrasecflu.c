/***********************************************************
*
*  intrasecflu.c
*
*  calculates total secundary fluoescence 
*  from the same layer as the primary fluorescence
*
*  M.Bos
*  oct 1997
*
**********************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

extern elem_line *enhances [MAX_ENH_LINES] ;
extern  elem_list lijst_els ;
double intrasecflu(char *meas_elem, char *meas_line, double lambda,
                   multilayer *sample, int nr_layer, double primint)
{
   int i, nr_enh, z_enh;
   char el_enh[3];
   elem_line *pline;
   char line[6];
   layer *player;
   double secint, tmpsecint, lambda_meas;

   if (primint == 0.0) return(0.0) ; /* quickly done */
   /* first check if line is fluoresced, if not no enhancement */

   lambda_meas=characwl(meas_line, meas_elem, &lijst_els);
   if ((lambda_meas <= lambda)||(primint==0.0)) return(0.0);


   player=get_layer(sample, nr_layer);
   nr_enh=find_enhanc(meas_elem, meas_line, player,lambda);
 /*  printf("Na find enhanc in intrasecflu nr_enh is %d for %s %s\n",
           nr_enh, meas_elem, meas_line);  */
   if (nr_enh==0) return(0.0); /* no enhancements */
   secint= 0.0;
 
   for (i=0; i < nr_enh; i++)
   {
     pline=enhances[i];
 /*    printf("pline is %p\n",pline); */
     z_enh=(*pline).z;
     
     strcpy(line, (*pline).line);
     strcpy(el_enh, get_asymb(z_enh, &lijst_els));
     tmpsecint=sijlambda(meas_elem, el_enh, meas_line, line, player, lambda,
                         primint);
/*    printf("in intrasecflu enhancmnt of  %s %s by line %s %s at labda %f is %e \n",
             meas_elem, meas_line, el_enh, line, lambda, tmpsecint); */ 
     secint += tmpsecint;
   }
   clean_up_enh(nr_enh);
   return(secint);
}
                
