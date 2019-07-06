/*************************************************
*
*
*   COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*   costerkr.c
*  
*   routine to build datastructure for
*   Coster Kornig yields
*
*   data from
*   M.O. Krause
*   J.Phys.ref.Data 8 (1979) 307
*
*   M.Bos
*   november 1997
*
******************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"


void bldfij(FILE *fp1, FILE *fp2, double fij[][13])
{
    int i, j ;
    char regel[80];
    for (i=12; i < 95; i++)  /* all elements from Mg - Pu */
    {
       fscanf(fp1,"%s",regel);  /* element name */
       for (j=0; j < 3 ; j++)
         fscanf(fp1,"%lf",&fij[i][j]);
    }
   for (i=20; i < 91; i++)  /* elements Ca- Th */
   {
      fscanf(fp2,"%s",regel); /* element name */
      for (j=3; j < 13; j++)
        fscanf(fp2,"%lf",&fij[i][j]);
   }

   /* debugging  
   printf(" bld fij \n");
   for (i=12; i <95 ; i++)
   {
     for (j=0; j < 13 ; j++)
          printf("%g ",fij[i][j]);
     printf("\n");
   } /* */
}

