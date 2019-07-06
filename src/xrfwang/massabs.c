/*
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*       file massabs.c
*
*       calculates mass absorption coefficient
*       of sample in a given file
*       for a characteristic line
*       or for a given lambda
*
*       M.Bos
*       april 1998
*
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"
#include "xrfdefs.h"

int main(void)
{
   extern multicomplayer compsample ;   /* sample in terms of compounds */
   extern multilayer sample ;          /* sample in terms of elements */
   double lambda_in, mu , energy;
   int i, j ;
   char el[3];
   char lin[6];
   FILE *fp;
   layer *plaag;
   char filename[80]; 
   char dirname[80]="./data/";
   extern char title[132]; 

   /* debugging
   printf("voor init tables \n"); /* */
   init_tables();
   /* debugging
   printf("na init tables \n"); /* */


   printf("Give filename with sample data ");
   gets(filename) ;
   strcpy(dirname+7,filename);
   fp=own_fopen(dirname);
   bldcompos(fp, &compsample);
   fclose(fp);
   convcomp(&compsample, &sample);

   printf("(L)ambda of (C)haracteristic line?\n");
   gets(filename);
   if (filename[0] == 'L')
   {
       printf("Give lambda in aengstrom ");
       scanf("%lf", &lambda_in);
   }
   else
   {
       printf("Give element symbol ");
       gets(el);
       printf("Give line symbol ");
       gets(lin);
       lambda_in=characwl(lin, el, &lijst_els);
   }
   energy=12.396/lambda_in;
   plaag=(sample.lagen)[0];
   mu=calc_mumix(energy, plaag, &lijst_els, &abstabel);
   printf("Mass abs. coeff. is %g \n", mu);
}


