/*************************************************
*
*   edges.c
*  
*   routine to build datastructure for
*   adsorption edges energies for element 2 -94
*
*   data from
*   T.P. Thinh, J. Leroux
*   X-Ray Spectrom. 8 (1979) 85-91
*
*   M.Bos
*   november 1997
*
******************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"

extern double edges[TOTAL][10] ;

void bldedges(FILE *fp, double edges[][10])
{
    int i, j ;
    char regel[120];
    char **pntregel ;
    char *pntr;
  
    pntregel=&pntr;
   /* 2-He : 10-Ne */

   for (i=2; i < 11; i++)
   {
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /*  line  Z = .. */
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /* element name  */
      fgets(regel,120,fp);  /* data */
      for (j=0; j < 1 ; j++)
        sscanf(regel,"%lf",&edges[i][j]);
      for (j=1 ; j < 10; j++)
         edges[i][j]=0.0 ;       /* unused */
   }

   /* 11-Na : 27-Co */
   for (i=11; i < 28; i++)
   {
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /*  line  Z = .. */
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /* element name  */
      fgets(regel,120,fp);
      pntr=regel;
      for (j=0; j < 2 ; j++)
       {
         edges[i][j]=strtod(pntr, pntregel);
       }
      for (j=2 ; j < 10; j++)
         edges[i][j]=0.0 ;       /* unused */
   }


   /* 28-Ni : 29-Cu */
   for (i=28; i < 30; i++)
   {
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /*  line  Z = .. */
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /* element name  */
      fgets(regel,120,fp);  /* data */
      pntr=regel;
      for (j=0; j < 3 ; j++)
      {
        edges[i][j]=strtod(pntr,pntregel);
      }
      for (j=3 ; j < 10; j++)
         edges[i][j]=0.0 ;       /* unused */
   }

   /* 30-Zn : 51-Sb */
   for (i=30; i < 52; i++)
   {
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /*  line  Z = .. */
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /* element name  */
      fgets(regel,120,fp);
      pntr=regel;
      for (j=0; j < 5 ; j++)
      {
         edges[i][j]=strtod(pntr,pntregel);
      }
      for (j=5 ; j < 10; j++)
         edges[i][j]=0.0 ;       /* unused */
   }

   /* 52-Te : 54-Xe */
   for (i=52; i < 55 ; i++)
   {
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /*  line  Z = .. */
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /* element name  */
      fgets(regel,120,fp);
      pntr=regel;
      for (j=0; j < 6 ; j++)
      {
        edges[i][j]=strtod(pntr,pntregel);
      }
      for (j=6 ; j < 10; j++)
         edges[i][j]=0.0 ;       /* unused */
   }

   /* 55-Cs */
   for (i=55; i < 56; i++)
   {
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /*  line  Z = .. */
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /* element name  */
      fgets(regel,120,fp);
      pntr=regel;
      for (j=0; j < 7 ; j++)
      {
        edges[i][j]=strtod(pntr,pntregel);
      }
      for (j=7 ; j < 10; j++)
         edges[i][j]=0.0 ;       /* unused */
   }

   /* 56-Ba : 60-Nd */
   for (i=56; i < 61 ; i++)
   {
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /*  line  Z = .. */
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /* element name  */
      fgets(regel,120,fp);
      pntr=regel;
      for (j=0; j < 8 ; j++)
      {
        edges[i][j]=strtod(pntr,pntregel);
      }
      for (j=8 ; j < 10; j++)
         edges[i][j]=0.0 ;       /* unused */
   }

   /* 61-Pm : 94-Pu */
   for (i=61; i < 95; i++)
   {
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /*  line  Z = .. */
      fgets(regel,80,fp);   /* blank line */
      fgets(regel,80,fp);   /* element name  */
      fgets(regel,120,fp);
      pntr=regel;
      for (j=0; j < 10 ; j++)
      {
        edges[i][j]=strtod(pntr,pntregel);
      }
   }
  
 /* debugging info 

   for (i=2; i < 95; i++)
   {
    for (j=0; j <  10; j++)
     printf("edges[%d][%d] = %e ",i,j, edges[i][j]);
    printf("\n");
   } /* */
}
