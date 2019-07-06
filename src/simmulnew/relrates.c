/******************
*
*     relrates.c
*     vs 2.0
*
*     routines to build
*     relative emission rates data structure
*     plus its access routine
*
*     M.Bos
*     october 1997
*
*     K-line , L-1 L-2 L3 data from
*     S.I. Salem, S.L. Panossian, R.A. Krause,
*     Atomic Data and Nuclear Data Tables 14 (1974) 91-109
*
*     Data for  M-lines  from
*     E.P. Bertin, "Principles and Practice of X-Ray
*     Spectrometric Analysis", 2nd edn, Plenum New York, 1975
*     Appendix 1
*
*
***************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"

extern char *linenames[27];
extern FILE *own_fopen(char *);

void bldrelrates( double relrates[][TOT_RELRATES])
{
   char regel[LINELENGTH];
   int i,j, k_l_lines ;
   FILE *fp ;

   /* k-lines: ka, ka1, ka2, kb1 kb3, kb2 */

  fp=own_fopen("krelrates.dat");

  for (j=3; j < 101; j++)
  {
   fscanf(fp,"%s",regel);
/*   printf("relrates %s \n",regel); */
   for (i=0; i < NR_KLINES; i++)
   {
      fscanf(fp, "%lf",&relrates[j][i]);
   }
  }
  fclose(fp);

   /* l1-lines: lb3, lb4, lg2, lg3, lg4 */

   fp=own_fopen("l1relrts.dat");
   for (j=23; j < 96; j++)
   {
     fscanf(fp,"%s",regel);
     for (i=0; i < NR_L1LINES-1; i++)
     {
       fscanf(fp, "%lf", &relrates[j][NR_KLINES+i]);
     }
     relrates[j][NR_KLINES+NR_L1LINES-1]=0.0;
   }
   fclose(fp);

   /* l2-lines: lb1, lg1, lg6, leta */
   
   fp=own_fopen("l2relrts.dat");
   for (j=16; j < 97; j++)
   {
     fscanf(fp,"%s",regel);
     for (i=0; i < NR_L2LINES; i++)
     {
       fscanf(fp, "%lf", &relrates[j][NR_KLINES+NR_L1LINES+i]);
     }
   }
   fclose(fp);


   /* l3 -lines: la1, la2, lb2, lb5, ll */

   fp=own_fopen("l3relrts.dat");
   for (j=21; j < 95; j++)
   {
     fscanf(fp,"%s",regel);
     for (i=0; i < NR_L3LINES; i++)
     {
       fscanf(fp, "%lf", &relrates[j][NR_KLINES+NR_L1LINES+
                                      NR_L2LINES+i]);
     }
   }
   fclose(fp);

    fp=own_fopen("relrates.dat");
    k_l_lines = NR_KLINES + NR_L1LINES + NR_L2LINES + NR_L3LINES;

    /* M2-line */

    fscanf(fp, "%lf", &relrates[0][k_l_lines]);
    
    /* M3-line */
    fscanf(fp, "%lf", &relrates[0][k_l_lines+1]);

    /* M4 -lines : mb, mx2 */

    for (i=0; i < NR_M4LINES; i++)
    {
       fscanf(fp, "%lf",&relrates[0][k_l_lines+2+i]);
    }

    /* M5 -lines : ma1, ma2, mx1 */
    
    for (i=0; i < NR_M5LINES; i++)
    {
       fscanf(fp, "%lf", &relrates[0][k_l_lines+2+NR_M4LINES+i]);
    }
    fclose(fp);
    for (i=1; i <94 ; i++)
      for (j=0; j < TOT_RELRATES; j++)
      {
       /* debugging
       printf("relrates[%d][%d]  is %e\n", i,j, relrates[i][j]); /* */
      }
}


double getrelrate(int z, char linesymb[], double relrates[][TOT_RELRATES])
{
    int i;
    i=0;
    while ( (strcmp(linesymb,linenames[i]) != 0)&& (i <= TOT_RELRATES))
             i++;

    if ( i > TOT_RELRATES)
    {
        printf(" %s unknown linbesymbol \n",linesymb);
        exit(1);
    }
    if (i < 20)
       return(relrates[z][i]);   
    else
     {
      return(relrates[0][i]);  /* voorlopig alles zelfde Z=0 !!!!! */
     }
}
 
