/*******************************************************
*
*   Copyright (C) 1999 M.Bos, for details dee file copying
*
*   File bscatter.c
*
*   calculates backscatter coefficient for use
*   in calculating tubespectra
*
*   cf H.Ebel, P.A. Pella
*   Advances X-Ray Spectrom. 35 (1992) 721
*
*   coded by M.Bos
*   University of Twente
*   Feb 1999
*
********************************************************/
#include <stdio.h>
#include <math.h>
#include "xrfluor.h"

double backscatter(int z, double u0)
{
   int i,j ;
   double result, t1 ;
   double  a[6][6];

   a[1][1]= 0.5580848699E-2;
   a[1][2]=0.2709177328E-3;
   a[1][3]=-0.5531081141E-5;
   a[1][4]=0.5955796251E-7;
   a[1][5]=-0.3210316856E-9;    
   a[2][1]=0.3401533559E-1;
   a[2][2]=-0.1601761397E-3;
   a[2][3]=0.2473523226E-5;
   a[2][4]=-0.3020861042E-7;
   a[3][1]=0.9916651666E-1;
   a[3][2]=-0.4615018255E-3;
   a[3][3]=-0.4332933627E-6;
   a[4][1]=0.1030099792;
   a[4][2]=-0.3113053618E-3;
   a[5][1]=0.3630169747E-1;

   t1=0.0 ;

   for (j=1;j <=5 ; j++)
     for (i=1; i <=j; i++)
      {
         t1+= a[i][j-i+1]*pow((1.0/u0 - 1.0),i)*pow(z,(j-i+1));
      }
   result=t1+1.0;
      
   return(result);
}
 
