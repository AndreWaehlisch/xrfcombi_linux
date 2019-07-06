/************************************************************
*  
*  Copyright (c) 1999 M.Bos, for details see file copying
*
*  file getgifct.c   
*
*  reads file ./data/alllines.dat
*  and extracts the Pi-values according to Wang,
*  X-Ray Spectrom. 25 (1996) 245-257
*
*  M.Bos
*  Sep 1999
*
*************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "xrfluor.h"

extern meas_line lines_to_meas[MAX_LINES] ;
extern elem_list lijst_els;
extern char anode_el[3];
extern double berryliumd, take_off;
extern double *pi_facts, *rik_meas;
extern meas_line lines_to_meas[MAX_LINES];

int get_nr_lines_meas(FILE *fp)
{
   int  line_nr ;
   char regel[LINELENGTH] ;
   line_nr=0;   
   while ( fgets(regel,LINELENGTH,fp) != NULL)
   {
       /* debugging
        printf("%s\n",regel); /* */
        line_nr++;
   }
  /* debugging 
   printf("Found in get_nr_lines_meas %d lines\n", line_nr);  /* */
   return(line_nr);
}

int  get_gifactors(FILE *fp)
{

    char *position ;
    char regel[LINELENGTH];
    int z, line_nr;
    char el_name[3];
    char line_symb[6];
    char filt_el[5];
    char xtal_name[7];
    double filt_mth;
    double tube_volts ;
    double   pi, kcps;
    
   line_nr=0;   
   while ( fgets(regel,LINELENGTH,fp) != NULL)
   {
       /* debugging  
       fputs(regel, stdout);  /* */
       position=strtok(regel," ");
       strcpy(el_name, position);
       /* debugging  
       printf("Element %s \n",position); /* */
       z=get_z(position, &lijst_els); 
      
        /* get line_name  */
        position=strtok(NULL," ");
        /* debugging  
        printf("line symbol %s \n",position); /* */
       
        strcpy(line_symb, position);
       
        position=strtok(NULL," ");
        sscanf(position,"%lf",&tube_volts);
        
        position=strtok(NULL," ");
        strcpy(filt_el,position); 
        position=strtok(NULL," ");
        sscanf(position,"%lf",&filt_mth);

        position=strtok(NULL," ");
        strcpy(xtal_name,position);

        position=strtok(NULL," ");  /* mask name */

        position=strtok(NULL," ");
        sscanf(position,"%lf",&kcps);
         
        /* debugging 
        printf("kcps is %e \n", kcps); /* */

        position=strtok(NULL," \n");
        sscanf(position,"%lf",&pi);

       /*  print and store pi-factors and meas intensities */
       
       /* debugging
       printf("results for line nr %d pi=%e kcps=%e\n",
               line_nr, pi, kcps);  /* */
       pi_facts[line_nr]=pi;
       rik_meas[line_nr]=kcps;

       /* store all relevant results in lines_to_meas structure */
       (lines_to_meas[line_nr]).z=z;
       strcpy((lines_to_meas[line_nr]).line, line_symb);
       (lines_to_meas[line_nr]).kev=tube_volts ;
       if ( (strcmp("None", filt_el)!=0))
         {
          /* debugging 
           printf("filt_el is %s van line_nr %d mth %f Z filt %d\n", 
                  filt_el, line_nr, filt_mth, get_z(filt_el,&lijst_els)); /* */
           (lines_to_meas[line_nr]).mth_filt=filt_mth;
           (lines_to_meas[line_nr]).filt_z=get_z(filt_el,&lijst_els);
         }
       else
         {
           /* no filter */
           (lines_to_meas[line_nr]).mth_filt=0.0;
         }
       /* debugging 
       printf("filt el Z %f filt mth %d \n",
        (lines_to_meas[line_nr]).mth_filt, 
          (lines_to_meas[line_nr]).filt_z); /* */
       line_nr++;
   }
      --line_nr;
   /* debugging 
   printf("ref-line from get_gifactors is %d\n", line_nr); /* */
   return(0 ); /* for now first line is the reference-line */
}

void calc_rikmeas(double *ins,  int nr_lines, int ref_line)
{
    int i ;
    
    for (i=0; i < nr_lines; i++)
      if ( i != ref_line)
       {
        ins[i]= ins[i]/ins[ref_line] ;
       }
    ins[ref_line]=1.0 ;

}

