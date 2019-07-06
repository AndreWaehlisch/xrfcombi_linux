/***********
/
*      characwl.c
*
*      routines to determine
*      characteristic wavelength for 
*      line and element
*      both given as symbols
*
*      M.Bos
*      sept 1997
*
*************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"

extern char *lines[27];
extern double klines[98][6];
extern double llines[98][14];
extern double mlines[98][7];
 

void bld_k_lines(FILE *fp, double klines[][6])
{
    int i;
    char regel[LINELENGTH];
    char *position ;
  

    i=2 ;  /* table starts at 3-Li */
   
    while (fgets(regel,LINELENGTH,fp) != NULL)
    {
      position=strtok(regel, " ");      /* elem symbol */
      /* debugging
      printf("%s \n",position);        /* debug info */
      position=strtok(NULL, " " );      /* K-alfa   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&klines[i][0]);
      position=strtok(NULL, " " );      /* K-alfa1   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&klines[i][1]);
      position=strtok(NULL, " " );      /* K-alfa2   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&klines[i][2]);
      position=strtok(NULL, " " );      /* K-beta1   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&klines[i][3]);
      position=strtok(NULL, " " );      /* K-beta3   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&klines[i][4]);
      position=strtok(NULL, " " );      /* K-beta2   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&klines[i][5]);
      i++ ;
    }
}
    
void bld_l_lines(FILE *fp, double llines[][14])
{
    int i;
    char regel[LINELENGTH];
    char *position ;

    i=15 ;  /* table starts at 16-S */
   
    while (fgets(regel,LINELENGTH,fp) != NULL)
    {
      position=strtok(regel, " ");      /* elem symbol */
      /* debugging
      printf("%s \n",position);        /* debug info */
      position=strtok(NULL, " " );      /* L-alfa1   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][0]);
      position=strtok(NULL, " " );      /* L-alfa2   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][1]);
      position=strtok(NULL, " " );      /* L-beta1   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][2]);
      position=strtok(NULL, " " );      /* L-beta2   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][3]);
      position=strtok(NULL, " " );      /* L-beta3   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][4]);
      position=strtok(NULL, " " );      /* L-beta4   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][5]);
      position=strtok(NULL, " " );      /* L-beta5   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][6]);
      position=strtok(NULL, " " );      /* L-gamma1   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][7]);
      position=strtok(NULL, " " );      /* L-gamma2   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][8]);
      position=strtok(NULL, " " );      /* L-gamma3   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][9]);
      position=strtok(NULL, " " );      /* L-gamma4   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][10]);
      position=strtok(NULL, " " );      /* L-gamma6   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][11]);
      position=strtok(NULL, " " );      /* L-l   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][12]);
      position=strtok(NULL, " " );      /* L-eta   */
      /* debugging
      printf("leta %s \n",position);         /* debug info */
      sscanf(position,"%lf",&llines[i][13]);
      i++ ;
    }
}


void bld_m_lines(FILE *fp, double mlines[][7])
{
    int i;
    char regel[LINELENGTH];
    char *position ;

    i=34 ;  /* table starts at 35-Br */
   
    while (fgets(regel,LINELENGTH,fp) != NULL)
    {
      position=strtok(regel, " ");      /* elem symbol */
      /* debugging
      printf("%s \n",position);        /* debug info */
      position=strtok(NULL, " " );      /* M-alfa1   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&mlines[i][0]);
      position=strtok(NULL, " " );      /* M-alfa2   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&mlines[i][1]);
      position=strtok(NULL, " " );      /* M-beta   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&mlines[i][2]);
      position=strtok(NULL, " " );      /* M-gamma   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&mlines[i][3]);
      position=strtok(NULL, " " );      /* M-Xi1   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&mlines[i][4]);
      position=strtok(NULL, " " );      /* M-Xi2   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&mlines[i][5]);
      position=strtok(NULL, " " );      /* M-noname   */
      /* debugging
      printf("%s \n",position);         /* debug info */
      sscanf(position,"%lf",&mlines[i][6]);
      i++ ;
    }
}


double characwl(char *linesymb, char *elemsymbol, elem_list *pelemlist)
{
   int nr_in_series, i, shell, z ;
   i=0;

   while ( strcmp(linesymb, lines[i]) != 0)
   {
     if ( i < 27) i++;
     else
       {
         printf("Non-existing line in characwl is %s \n", linesymb);
         exit(1);
       }
   }
   if ( i < 6) 
   {
      nr_in_series=i ;
      shell = 1;
   } else
   { 
      if ( i < 20)
      {
        nr_in_series= i-6;
        shell = 2;
      } else
        {
          nr_in_series = i-20;
          shell= 3;
         }
    }


      z=get_z(elemsymbol,pelemlist);  

  switch (shell)
   {
     case 1 :
      return( klines[z-1][nr_in_series]);
     break;
     case 2 : 
              return( llines[z-1][nr_in_series]);
     break;
     case 3 :
             return( mlines[z-1][nr_in_series]);
     break;
     default : printf("Inconsistency in detmn charac wavelength\n");
     exit(1);
   }
}

 
