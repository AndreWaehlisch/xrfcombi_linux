/************************************************
*
*
*   COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*   nkls.c
*
*   routine to build datastructure
*   with nkl1 nkl2 nkl3 values
*   from D.K.G. de Boer
*   SpectroChim Acta, 44B (1989) 1171
*
*   M.Bos
*   november 1997
*
************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"

void bldnkls(FILE *fp, double nkl123[][3])
{
    int i, j ;

    for (i=20; i < 95 ; i++)  /* from Ca - Pu */
    {
       for (j=0; j < 3 ; j++)
       {
         fscanf(fp,"%lf",&nkl123[i][j]);
       }
     }
   /* debugging 
   printf("In nkls \n"); 
   for (i=20; i < 95; i++)
   {
    for (j=0; j < 3 ; j++)
      printf("%lf ",nkl123[i][j]);
    printf("\n");
   } /* */
}



