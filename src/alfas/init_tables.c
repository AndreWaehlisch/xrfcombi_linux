/*****************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*      init_tables.c
*
*      routine to load all needed tables
*      for XRFLUOR from files
*
*      data is stored in global vars that are
*      defined in xrfdefs.h
*
*      M.Bos
*      sep 1997
*
*      
********************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"


extern elem_list lijst_els ;
extern double omegas[TOTAL][NR_OMEGAS];
extern double jmpdata[ALL_ELMS][2];
extern double relrates[TOTAL][TOT_RELRATES];
extern double klines[98][6];
extern double llines[98][14];
extern double mlines[98][7];
extern emudata abstabel;
extern double edges[TOTAL][10];
extern double nkl123[TOTAL] [3];
extern double fij[TOTAL][13];

FILE *own_fopen(char *);
 
void init_tables(void)
{
  FILE *fp1, *fp2, *fp3 ;
  
  /* element gegevens */
  fp1=own_fopen("periotbl.dat");
  bld_atlist(fp1,&lijst_els);
  fclose(fp1);

  /* fluorescence yields omega */
  fp1=own_fopen("omkl.txt");
  fp2=own_fopen("zm45.txt");
  fp3=own_fopen("zm23.txt");
  bldomega(fp1,fp2,fp3,omegas);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fp1=own_fopen("edges.dat");
  bldedges(fp1,edges);
  fclose(fp1);
  fp1=own_fopen("nkl123.dat");
  bldnkls(fp1, nkl123);
  fclose(fp1);
  fp1=own_fopen("costerkr.dat");
  fp2=own_fopen("costerkm.dat");
  bldfij(fp1,fp2, fij);
  fclose(fp1);
  fclose(fp2);
  /* jump factors at edges */
  /* debugging
  printf("Voor openen appendix8.dat \n"); /* */
  fp1=own_fopen("appendix8.txt");
  bldjmpedg(fp1,jmpdata);
  fclose(fp1);
  /* debugging
  printf("Na bldjmpedg \n");  /* */
 
  /* relative emission rates */
  bldrelrates(relrates);

  /* characteristic wavelengths */
  /*  in angstrom */
  fp1=own_fopen("kser.dat");
  bld_k_lines(fp1,klines);
  fclose(fp1);
  fp1=own_fopen("lser.dat");
  bld_l_lines(fp1,llines);
  fclose(fp1);
  fp1=own_fopen("mser.dat");
  bld_m_lines(fp1,mlines);
  fclose(fp1);
  /* debugging
  printf("golflengte leta uranium is %g \n", llines[91][13]);  /* */
  /* emu data */
 fp1=own_fopen("emutable.dat");
 bldemutable(fp1,&abstabel); 

}



FILE *own_fopen(char *filenaam) 
{
    FILE *fp;
    if ((fp=fopen(filenaam,"r"))==NULL)
    {
         printf("Cannot open file %s \n",filenaam);
         exit(1);
    }
    return(fp);
}


