/********************************************************
*
*    COPYRIGHT (C) M.BOS 1999, FOR DETAILS SEE COPYING
*
*
*    File tstfitcomp.c
*
*    test program for routine
*    to calculate compound composition
*    for a mixture with known compounds
*    from given elemental composition
*
*    M.Bos
*    September 1999
*
*********************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"
#include "xrfdefs.h"
  
double *free_pars;
double *pars_min;
double *pars_max ;


int main(int argc, char *argv[])
{
  extern multicomplayer compsample;
  extern multilayer sample;
  multicomplayer nep;
  int i,j, k ;
  FILE *fp, *fpconv ;
  char filename[80]="ksample.dat" ;
  char filename2[80]="usample.dat" ;

  init_tables();
  free_pars= malloc((nr_free_vars+1)*sizeof(double));
  pars_min= malloc((nr_free_vars+1)*sizeof(double));
  pars_max=malloc((nr_free_vars+1)*sizeof(double));

  fp=fopen(filename,"r");
  bldcompos(fp,&nep);
  fclose(fp);
  convcomp(&nep,&sample); 
  fp=fopen(filename2,"r");
  bldcompos(fp,&compsample); 
  fclose(fp);
  fitconc(&compsample, &sample); 
  exit(0);
}

