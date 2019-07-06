/***********************************************************
*
*  Copyright (c) M.Bos 1998, for details see file copying  
*
*  intersecflu.c
*
*  calculates total secundary fluoescence 
*  from one other layer 
*
*  M.Bos
*  sep 1998
*  
* see D.K.G. de Boer, X-Ray Spectrom. 19 (1990) 145
*
**********************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

extern elem_line *enhances [MAX_ENH_LINES] ;
extern  elem_list lijst_els ;
double intersecflu(char *meas_elem, char *meas_line, double lambda,
                   multilayer *sample, int nr_curr_layer,
                   int nr_other_layer,  double primint)
{
   int i, nr_enh, z_enh;
   char el_enh[3];
   elem_line *pline;
   char line[6];
   layer *pcurrlayer, *potherlayer;
   double secint, tmpsecint, lambda_meas;

   /* first check if line is fluoresced, if not no enhancement */

   lambda_meas=characwl(meas_line, meas_elem, &lijst_els);
   if ((lambda_meas <= lambda)||(primint==0.0)) return(0.0);


   pcurrlayer=get_layer(sample, nr_curr_layer);
   potherlayer=get_layer(sample, nr_other_layer);
   nr_enh=find_enhanc(meas_elem, meas_line, potherlayer, lambda);

   /* debugging  
   printf("Na find enhanc in intersecflu nr_enh is %d for %s %s\n",
           nr_enh, meas_elem, meas_line);  /* */

   if (nr_enh==0) return(0.0); /* no enhancements */
   secint= 0.0;
 
   for (i=0; i < nr_enh; i++)
   {
     pline=enhances[i];

     /* debugging 
     printf("pline is %p\n",pline); /*  */

     z_enh=(*pline).z;
     
     strcpy(line, (*pline).line);
     strcpy(el_enh, get_asymb(z_enh, &lijst_els));

     if (nr_curr_layer < nr_other_layer)
     {
       tmpsecint=sij_lambda_up(meas_elem, el_enh, meas_line, line, 
                               sample, nr_curr_layer, nr_other_layer,
                               lambda, primint);

       /* debugging  
       printf("in intersecflu up sec enh is %g \n", tmpsecint ); /* */

     }
     else
     { 
       tmpsecint=sij_lambda_down(meas_elem, el_enh, meas_line, line, 
                               sample, nr_curr_layer, nr_other_layer,
                               lambda, primint);
       /* debugging 
       printf("in intersecflu down sec enh is %g \n", tmpsecint ); /* */
     }

    /* debugging 
    printf("in intersecflu enhancmnt factor of  %s %s by line %s %s at labda %f \
              in layer %d is %g \n",
             meas_elem, meas_line, el_enh, line, lambda, nr_other_layer,
              tmpsecint); /* */ 

     secint += tmpsecint;
   }
   /* debugging 
   printf(" in intersecflu secint is % g \n", secint);  /* */

   clean_up_enh(nr_enh);
   return(secint*primint);
}
                
