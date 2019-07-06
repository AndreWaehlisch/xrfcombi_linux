/***********************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*   sum_el.c
*
*   calculates sum of element fractions in a layer
*
*   M.Bos
*   April 1998
*
*************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

double sum_el(multilayer *sample, int layer_nr)
{
   int i ;
   layer *curr_layer;
   int curr_nr_el ;
   double sum_el_fracts;
   elem  *pelems;

   curr_layer= sample->lagen[layer_nr];  /* address curr layer */
   curr_nr_el= curr_layer->nr_elements;
   pelems=curr_layer->elementen;
   /* debugging 
   printf("aantal elementen in laag %d is %d \n", layer_nr, curr_nr_el); /* */
   sum_el_fracts=0.0;  
   for (i=0; i < curr_nr_el; i++)
   {
     sum_el_fracts += (pelems[i]).wfract;
     /* debugging  
     printf("wfract element %d is %g \n",i, (pelems[i]).wfract); /* */
   } 
   return(sum_el_fracts);
}

    
