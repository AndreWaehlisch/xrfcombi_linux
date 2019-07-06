/**********************************************
*
*    set_spectrom.c
*    enter spectrometer settings
*
*    M.Bos
*    sept 1997
*
***********************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "xrfluor.h"

FILE *own_fopen(char *);

extern double Gfactor;
extern double psi1;
extern double psi2;
extern char anode_el[3];
extern double  berryliumd ;
extern double take_off ;
extern int nr_of_meas_lines ;

void set_spectrom(void)
{
   FILE *fp;

   fp=own_fopen("./data/pwsettings.dat");
   fscanf(fp,"%s %lf %lf %lf %lf %lf %lf", anode_el,  &berryliumd,
             &take_off, &Gfactor, &psi1, &psi2);
   psi1=(psi1/360.0)*2.0*PI;
   psi2=(psi2/360.0)*2.0*PI;
   fclose(fp);
   fp=own_fopen("./data/deflines.dat");
   nr_of_meas_lines=get_mslines(fp);
   fclose(fp); 
}



