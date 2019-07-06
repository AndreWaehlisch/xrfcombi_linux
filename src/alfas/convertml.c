/***********************************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*                    convertml.c
*                    Main program for quantitative 
*                    X-ray fluorescence analysis
*
*                    M.Bos
*                    University of Twente
*                    Faculty of Chem. Technology
*                    Department of Chemical Analysis
*                   
*
*                    converts ./data/simsampl.dat to elemnt data
*                    Nov 20th 2002
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"
#include "xrfdefs.h"

FILE *own_fopen(char *filenaam) ;
extern elem_list  lijst_els ;

int main(int argc,char *argv[])
{
   extern multicomplayer compsample ;   /* sample in terms of compounds */
   extern multilayer sample ;          /* sample in terms of elements */
   char el[3];
   char lin[6];
   FILE *fp, *fp1;
   char filename[80]="simsampl.dat"; 
   char dirname[80]="./data/" ;
   extern char title[132]; 
   int nr_compounds;

   fp1=own_fopen("periotbl.dat");
   bld_atlist(fp1,&lijst_els);

/*   printf("Give filename with sample data ");
   gets(filename) ; */

   strcpy(dirname+7,filename);
   fp=own_fopen(dirname);
   bldcompos(fp, &compsample);
   /* debug */
   printf("Na bldcompos \n"); /* */
   fclose(fp);
   convcomp(&compsample, &sample);
   checkdata(&sample);
   return(0);
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


