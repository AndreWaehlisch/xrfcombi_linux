/******************
*
*     COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*     omega.c
*
*     routines to build
*     fluorescence yield datatable
*     plus its access routine
*
*     M.Bos
*     october 1997
*
*     omega K L1, L2 and L3 from Krause, M.O., J.Phys.Che.ref.Data,
*     vol 8 No 2, 1979, p 307
* 
*     omega M4 and M5 D.K.G. de Boer, SpectroChim. Acta Vol 44B, 1989,
*     pp 1171-1190 
*
*
*     omega M2 and M3 lines: W.Bambynek et.al.
*     Rev. Modern Phys. 44 (1972) 716
*
***************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"


void bldomega(FILE *fpkl, FILE *fpm45, FILE *fpm23, double omegas[][NR_OMEGAS])
{
   char regel[LINELENGTH];
   double omegak, omegal1, omegal2, omegal3, omegam4, omegam5 ;
   double omegam2, omegam3 ;
   int i,j ;
   char *positie ;
   /* zero array omega */
   for (i=0; i < TOTAL ; i ++)
      for (j=0; j < NR_OMEGAS; j++)
        omegas[i][j]=0.0;
   j=5;
   while ((fgets(regel,LINELENGTH,fpkl) != NULL)&& (j< TOTAL))
   {
     positie=strtok(regel," \t");  /* To Z of entry */
     /* debugging
     printf("%s \n",positie); /* debugging */
     positie=strtok(NULL," \t");  /* To element symbol */
     /* debugging
     printf("%s \n", positie); /* debugging */
     positie=strtok(NULL," \t");  /* to omega-K */
     sscanf(positie,"%lf",&omegak);
     /* debugging
     printf("omega-k is %f \n",omegak); /* */
     positie=strtok(NULL," \t"); /* to omega L1 */
     sscanf(positie,"%lf",&omegal1);
     positie=strtok(NULL," \t"); /* L2 */
     sscanf(positie,"%lg",&omegal2);
     positie=strtok(NULL," \t"); /* L3 */
     sscanf(positie,"%lg",&omegal3);
     /* debugging
     printf("omega-l1 is %g omega-l2 %g omega-l3 %g\n",
             omegal1, omegal2,omegal3); /* */
     /* now put these in the datastructure */
     omegas[j-1][0]=omegak ;
     omegas[j-1][1]=omegal1 ;
     omegas[j-1][2]=omegal2 ;
     omegas[j-1][3]=omegal3 ;
     j++ ;
   }
   j=57;
   while ((fgets(regel,LINELENGTH,fpm45) != NULL)&&( j < TOTAL))
   {
     positie=strtok(regel," \t");  /* To Z of entry */
     /* debugging
     printf("%s \n",positie); /* */
     positie=strtok(NULL," \t");  /* to omega-m4 */
     sscanf(positie,"%lf",&omegam4);
     /* debugging
     printf("omega-m4 is %f \n",omegam4); /* */
     positie=strtok(NULL," \t"); /* to omega m5 */
     sscanf(positie,"%lf",&omegam5);
     /* debugging
     printf("omega-m5 is %f \n", omegam5); /* */
     /* now put these in the datastructure */
     omegas[j-1][4]=0.0 ;  /* omegam2 */
     omegas[j-1][5]=0.0;  /*  omegam3 */
     omegas[j-1][6]=omegam4 ;
     omegas[j-1][7]=omegam5 ;
     j++ ;
   }
   j=20;
   while ((fgets(regel,LINELENGTH,fpm23) != NULL)&& ( j < 93))
   {
     positie=strtok(regel," \t");  /* to Z of entry */
     /* debugging
     printf("%s \n",positie); /* */
     positie=strtok(NULL," \t");  /* to omega-2 */
     sscanf(positie,"%lf",&omegam2);
     /* debugging
     printf("omega-m2 is %e\n",omegam2); /* */
     positie=strtok(NULL," \t");  /* to omega m3 */
     sscanf(positie,"%lf",&omegam3);
     /* debugging
     printf("omega-m3 is %e \n",omegam3); /* */
     /* now put these in the datastructure */
     omegas[j-1][4]=omegam2 ;  /* omegam2 */
     omegas[j-1][5]=omegam3;  /*  omegam3 */
     j++ ;
   }
}



double get_omega(int z, char *shell, double omegas[][NR_OMEGAS])
{
   int i;
   char *shellnames[]={ "k", "l1", "l2", "l3", "m2",
                         "m3", "m4", "m5"};

    i=0;
    while (strcmp(shell,shellnames[i])!=0)
    {
       if (i > 7)
       {
         printf("%s is a non existing shell \n",shell);
         exit(1);
        }
        else i++;
     }
    return( omegas[z-1][i]);
}
 
  
char *line_to_shell(char *line, char *linenames[27], char *shellnames[27])
{
   int i,j ;
   i=0;
   while (strcmp(linenames[i],line)!=0)
   {
           if (i < 27) i++;
           else 
              { 
                printf("Non-existing line in line to shell %s \n",
                        line);
               for(j=0; j < (int)strlen(line) ; j++)
               {
                 printf("%d ",line[j]); 
               }
                printf("\n");
                exit(1);

              }
    }
  return(shellnames[i]);
}

/*****************
*     
*   get_hole_nr()
*
*   returns 1 K shell
*           2 LI shell
*           3 LII shell
*           4 LIII shell
*           5 M2 shell
*           6 M3
*           7 M4
*           8 M5
****************************/

int get_hole_nr(char *shellname)
{
    int i ;
    char *shells[]={"k", "l1", "l2", "l3", "m2","m3","m4","m5"};
    i=0;
    while ((i < 9)&& strcmp(shellname, shells[i])) i++ ;   
    if (i==9)
     {
       printf("Unknown shellname %s \n",shellname);
       exit(1);
     }
    return (i+1);
}

