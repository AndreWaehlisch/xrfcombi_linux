/************************************************
*
*  COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*    enhanc.c
*
*    finds elements and their lines in a specific layer that
*    enhance primary radiation of a measured element line
*
*    M.Bos
*    oct 1997
*
*    corrected october 1998 with condition
*    labda enhancing line > lambda   
* 
*    corrected for use absorption edge instead of lambda
*    April 2001
*
*****************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include "xrfluor.h"

extern elem_list lijst_els ;
extern char *linenames[27];
extern char *lines[27];
extern char *shellreeks[27] ;
extern double klines [98][6];
extern double  llines [98][14];
extern double  mlines [98][7];
extern elem_line *enhances[MAX_ENH_LINES];
extern double edges[TOTAL][10];

int find_enhanc(char *el_i, char *line_i, layer *player, double lambda)
{
   int i, j, zi, zj ;
   int nr_elements_in_layer, enh_count;
   elem_line *pline ;
   double labda_i, edge_i, edge_j  ;
   char *shellname_i, *shellname_j ;
   int hole_nr ;
   
   labda_i=characwl(line_i, el_i, &lijst_els);
/*   printf("in find_enhanc el is %s \n",el_i); */
   zi= get_z(el_i,&lijst_els);
   nr_elements_in_layer=(*player).nr_elements ;
   shellname_i=line_to_shell(line_i, linenames, shellreeks);
   hole_nr = get_hole_nr(shellname_i);
   edge_i= 12.396/edges[zi][hole_nr-1];
    
/*
*    now we know the element and the characteristics
*    of the measured line
*
*    now start to find enhancing lines
*    based on higher Z same lines from series will enhance 
*/
  enh_count = 0;
  for (i=0; i < nr_elements_in_layer ; i++)
  {
     zj=get_z( ((*player).elementen[i]).name, &lijst_els);
     if (zj != zi)
     {
     /* now k - shell enhancing element */
     for (j=0; j < 6 ; j++)  
     {
      if (( zj < 13)|| (j != 0)) /* from Aluminium ka1 and ka2 lines listed */
       if ((klines[zj-1][j] < edge_i) &&
                    (12.396/edges[zj][0] > lambda) 
                     &&( klines[zj-1][j] != 0))
        {
          pline=(elem_line *) malloc(((sizeof(elem_line)+15)/16)*16);
          enhances[enh_count]=pline; 
          enh_count++;
          (*pline).z= zj;
          strcpy( (*pline).line, lines[j]);
           /* printf("enhanced by Z = %d, line %s, labda is %f \n",
                  zj,lines[j],klines[zj-1][j]); /*  */
        }
      }
         /* now l -shell */
        for (j=0; j < 14 ; j++)
        {
           shellname_j=line_to_shell(lines[6+j],linenames,shellreeks);
           hole_nr=get_hole_nr(shellname_j);
           edge_j=12.396/edges[zj][hole_nr-1];
           if ((llines[zj-1][j] < edge_i) &&
               (edge_j > lambda) && (llines[zj-1][j] != 0)) 
           {
              pline=(elem_line *)malloc(((sizeof(elem_line)+15)/16)*16);
              if (pline==NULL)
              {
                 printf("Cannot alloc mem for elem in enhanc \n");
                 exit(1);
              }
              enhances[enh_count]=pline; 
              enh_count++;
              strcpy( (*pline).line,lines[6+j]);
              (*pline).z=zj;
           /*  printf("enhanced by Z = %d, line %s, labda is %f \n",
                     (*pline).z,lines[6+j],llines[zj-1][j]);/*  */
            }
         }
         /* now m -shell */
        for (j=0; j < 7 ; j++)
        {
           shellname_j=line_to_shell(lines[20+j],linenames,shellreeks);
           hole_nr=get_hole_nr(shellname_j);
           edge_j=12.396/edges[zj][hole_nr-1];
           if ((mlines[zj-1][j] < edge_i) &&
               (edge_j > lambda) && (mlines[zj-1][j] !=0))
           {
              pline=(elem_line *)malloc(((sizeof(elem_line)+15)/16)*16);
              enhances[enh_count]=pline; 
              enh_count++;
              (*pline).z=zj;
              strcpy( (*pline).line,lines[20+j]);
         /* printf("enhanced by Z = %d, line %s labda is %f \n",
              zj,lines[20+j],mlines[zj-1][j]); /*  */
            }
         } 
      }
  } 
/*  show_enhances(enh_count,enhances);  */
  return(enh_count); 
}


void clean_up_enh(int enh_count)
{
   int i;
   if (enh_count==0) return;
   for (i=0; i < enh_count ; i++)
    free(enhances[i]);
}

/**************************
*
*  debugging help
*
*********************************************/
void show_enhances(int count, elem_line **pnt)
{
   int i;
   for (i=0; i < count; i++)
   printf("z %d line %s \n", pnt[i]->z, pnt[i]->line);
}

