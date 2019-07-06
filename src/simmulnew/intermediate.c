/***********************************************************
*
* COPYRIGHT M.BOS (c) 1998, for details see file copying
*
* intermediate.c
*
* routines to calculate thickness of intermediate layer
* and mean absorption coefficient at lambda
*
* M.Bos
* UNIVERSITY of TWENTE
*
* sept 1998
*
*************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern emudata abstabel ;
extern elem_list lijst_els ;


double calc_th_b(multilayer *sample, int nr_curr_layer, int nr_other_layer)
{
    int i ;
    double  tot_thickn ;

    tot_thickn =0.0;
    if (nr_curr_layer < nr_other_layer)
    {
      for (i=nr_curr_layer+1; i < nr_other_layer; i++)
      {
         tot_thickn += (*(sample->lagen[i])).massthickness ;
      }
    }
    else 
    {
      if (nr_other_layer < nr_curr_layer)
      {
        for (i=nr_other_layer+1; i < nr_curr_layer; i++)
        {
           tot_thickn += (*(sample->lagen[i])).massthickness ; 
        }
      }
      else
      {
         printf("Error in calc_th_b with layer indices \n");
         exit(1);
      }
   }
	return(tot_thickn);
}
        

double calc_mean_mu( multilayer *sample,  int nr_curr_layer,
                     int nr_other_layer, double lambda)
{
  int i;
  double energy, dikte, mudtot;
  layer *player ;
  int start, end ;
  
  energy= 12.396/lambda ;

  if (nr_curr_layer < nr_other_layer)
  {
     start = nr_curr_layer+1;
     end= nr_other_layer;
  }
  else
  {
    start = nr_other_layer+1 ;
    end = nr_curr_layer;
  }

  mudtot=0.0 ;
  dikte = 0.0 ;

  for (i=start; i < end ; i++)
  {
    player= sample->lagen[i];
    mudtot+= calc_mumix(energy, player, &lijst_els, &abstabel)* 
             (*player).massthickness;
    dikte += (*player).massthickness ;
    /* debugging 
    printf("in calc_mu_mean mudtot %g dikte %g \n", mudtot, dikte); /* */
  }
   if (dikte==0.0)
     return(0.0);
   else
     return (mudtot/dikte);
}

