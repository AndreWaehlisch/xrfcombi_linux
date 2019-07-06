/*********************************************************
*
*    COPYRIGHT (C) M.BOS 1999 , for details see copying
*
*
*   respfuncnv.c
*
*   contains definition for funcnv()
*   that calculates response ( sum of squared errors)
*   for set of input parameters in compounds structure 
*   by calcn of elemental composition and determination
*   of squared error with found elemental composition
*
*   M.Bos
*   September 1999
*
*  included sum conc =1.0 in condx 22 Sep 1999
*  included rel squared error April 11th 2004
*************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

int funcnv(struct pstruct *parin)
{
   int i, j, k  ;
   double  ssq ;
   double sumfracts;  /* sum of fractions of compounds */
   extern multicomplayer compsample ;
   extern multilayer sample; 
   multilayer tmpsample ;
   extern elem_list lijst_els;
   extern int nr_free_vars;
   extern meas_line lines_to_meas[MAX_LINES] ;
   extern int nr_of_meas_lines ;
   int nr_layers, nr_elements ;
   layer *calc_layer, *found_layer ;
   elem *list_elem_calc, *list_elem_found ;
   double wfract_calc, wfract_found ;
   int elem_is_measured ;


  /* test for negative params -- if so return very large ssq */
   for (i=0; i <   nr_free_vars ; i++)
   {
      if ((*parin).parm[i] < 0.0 )     /* || ((*parin).parm[i] > 1.0)) */
      {
         /* debugging
          printf("Negative parameters found \n"); /* */
         (*parin).val=1.0e4;
         return(1);
      }
   }
 /* test for param s > 1.0 -- if so return very large ssq */

   for (i=0; i <   nr_free_vars ; i++)
   {
      if ((*parin).parm[i] > 1.0 )     /* || ((*parin).parm[i] > 1.0)) */
      {
         /* debugging
         printf("parameters > 1.0  found \n");  /* */
         (*parin).val=1.0e4;
         return(1);
      }
   }
  
  
   
/*   printf("Na testen op neg params in func\n");  */

   put_params(&compsample, (*parin).parm);   /* put params in datastruct smpl */

   ssq=0.0 ; 
   copysample(&sample, &tmpsample);
   sumfracts=convcomp(&compsample, &tmpsample);
   /* debugging  */
   printf("sumfracts is %g \n",sumfracts); /* */
   
   nr_layers=sample.nr_layers;
   for (i=0; i < nr_layers; i++)
   {
      calc_layer= tmpsample.lagen[i];

      nr_elements=calc_layer->nr_elements;
      found_layer= sample.lagen[i];
      list_elem_calc = calc_layer->elementen;
      list_elem_found = found_layer->elementen;

      for (j=0; j < nr_elements; j++)
      { 
       elem_is_measured=0;    /* start with not measured */
       for (k=0; k < nr_of_meas_lines; k++)
       {
         if ( (lines_to_meas[k]).z==get_z(list_elem_found[j].name, &lijst_els))
         { 
           elem_is_measured=1;
          /* debugging */
/*         printf("found %s is measured element \n", list_elem_found[j].name); /* */
         break ;
         
         }
       }

       if (elem_is_measured)
       {     
        wfract_calc=list_elem_calc[j].wfract;
        wfract_found= list_elem_found[j].wfract;
 /*       printf("debug wfraccalc %g wfracfound %g \n", wfract_calc,
                wfract_found); /*  */
       ssq += (wfract_calc-wfract_found)*(wfract_calc-wfract_found);
      /* added pa0mbo april 11th 2004 */
       ssq +=  (1.0-wfract_calc/wfract_found)*(1.0-wfract_calc/wfract_found); 
      }
      }
    }
    ssq += (1.0-sumfracts)*(1.0-sumfracts);  /* pa0mbo added 22-9-99 */
   (*parin).val=ssq ;
   /* debugging  
   printf("Debugging funcnv net voor return ssq %g \n",ssq); /* */
   free_multil(&tmpsample);
   return(0);  /* voorlopig altijd OK later ook ERROR */
}
