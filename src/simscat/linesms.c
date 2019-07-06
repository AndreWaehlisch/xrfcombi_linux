/********************************************
*
*   COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*   linesms.c
*
*   builds database of measured lines
*   from data file lines.dat
*
*   M.Bos
*   oct 1997
*
*****************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "xrfluor.h"

extern meas_line lines_to_meas[MAX_LINES] ;
extern elem_list lijst_els ;

int get_mslines(FILE *fp)
{
   int  z ,nr_mslines, act_length, z_filt ;
   char regel[LINELENGTH];
   char *position;
   double mth_filt, kev;
  
   nr_mslines=0;
   while (fgets(regel,LINELENGTH,fp) != NULL)
   { 
       act_length= strlen(regel);
       /* debugging
       printf("actuele lengte regel van lines.dat %s is %d\n",
               regel, act_length); /* */
       position=strtok(regel," ");        /* el symbol */
       /* debugging
       printf("in readin lines.dat %s\n",position); /8 */
       z=get_z(position, &lijst_els);
       /* debugging
       printf("Atomic number is %d \n",z); /* */
       position=strtok(NULL, " \n");       /* line symbol */
       /* debugging
       printf("in readin lines.dat %s\n",position); /* */
       (lines_to_meas[nr_mslines]).z = z;
       strcpy((lines_to_meas[nr_mslines]).line,position);
       position=strtok(NULL, " ");       /* HV in keV */
       sscanf(position,"%lf", &kev);
       /* debugging
       printf("HT is %g keV \n",kev); /* */
       (lines_to_meas[nr_mslines]).kev=kev;
       position=strtok(NULL," \n");
       if ( ((position-regel) <= act_length)&& (position != NULL) )
       {
         /* debugging
         printf("position %p  regel %p \n",position, regel);
         printf("Filter is %s \n",position); /* */
         z_filt=get_z(position,&lijst_els);
         (lines_to_meas[nr_mslines]).filt_z=z_filt;
         position=strtok(NULL, " \n");  /* to thickness */
         sscanf(position,"%lf", &mth_filt);
         /* debugging
         printf("FILTER Z: %d, massthickness %g \n",z_filt, mth_filt); /* */
         (lines_to_meas[nr_mslines]).mth_filt=mth_filt;
       }
       else
         (lines_to_meas[nr_mslines]).mth_filt=0.0;
       nr_mslines++ ;
    }
    return(nr_mslines);
}
      
