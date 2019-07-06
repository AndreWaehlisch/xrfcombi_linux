/*****************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*     file phase2.c
*     calculates elemental composition of a sample
*     based on iterating first equation in Table 1
*     from paper by R.M. Rousseau and M. Bouchard
*     X-Ray Spectrom. 15 (1986) 207
*
*     M.Bos
*     April 1998
*
***************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern double *rxi_meas;
extern double **alpha_ij, **rho_ij ;
extern  int nr_of_meas_lines ;
extern meas_line lines_to_meas[MAX_LINES];
extern elem_list lijst_els;
double attn_factor(int, int, elem *, int, int *);
double enh_factor(int, int, elem *, int, int *);


void phase2(multilayer *sample)
{
  int i, j ;
  int nr_elem;
  double delta;
  layer *ptop;
  elem *pelements;
  char meas_elem[3];
  int z_meas_elem;
  int *seq_elms_in_layer;
  double newconc, attnfactor, enhfactor ;
  int no_oxygen;
  no_oxygen=1 ; /* default suppose there is no oxygen */

  /* debugging 
  printf("Nr of measured lines is %d\n", nr_of_meas_lines); /* */
  ptop=(*sample).lagen[0];
  nr_elem= (*ptop).nr_elements ;
  pelements= (*ptop).elementen ;
  seq_elms_in_layer= (int *) malloc(nr_of_meas_lines*sizeof(int));

  /* set up correspondence of elements to measured line */

  for (i=0; i < nr_of_meas_lines ; i++)
  {
     z_meas_elem= (lines_to_meas[i]).z ;
     strcpy(meas_elem, get_asymb(z_meas_elem, &lijst_els));
     /* debugging 
     printf("z_meas_elem %d symbol %s\n", z_meas_elem, meas_elem); /* */
     j=0;
     while (( strcmp(meas_elem, (pelements[j]).name) !=0 )&& (j <= nr_elem)) 
     {
        /* debugging 
        printf("%s \n", (pelements[j]).name);  /* */
        j++;
     }
     if (j> nr_elem)
     { 
      printf("error in finding element %s  belonging to line %d\n",
              meas_elem, i);
      exit(1);
     }
     seq_elms_in_layer[i]=j;     
   }


  do
  {
    delta =0.0;

    for (i=0; i < nr_of_meas_lines; i++)
    {
      /* debugging 
      printf("phase2 line nr is %d\n",i); /* */
      attnfactor=attn_factor(i, nr_of_meas_lines,pelements,
                                nr_elem, seq_elms_in_layer);
      newconc =rxi_meas[i]*attnfactor;
      /* debugging 
      printf("attn is %g \n",attnfactor); /* */
      enhfactor=enh_factor(i, nr_of_meas_lines, pelements,
                                nr_elem, seq_elms_in_layer);
      newconc /=enhfactor;
      /* debugging 
      printf("enh_fact is %g \n",enhfactor); /* */
      delta+= fabs(newconc- (pelements[seq_elms_in_layer[i]]).wfract);
      (pelements[seq_elms_in_layer[i]]).wfract=newconc;
      /* debugging 
      printf("newconc %d %g \n", i, newconc);  /* */
    }
  } while (delta > 1.0E-10) ;
  printf("----------------------------------------------\n");
  printf("Estimated conc. after 2nd step rousseau algorithm\n");
  for (i=0; i < nr_of_meas_lines;i++)
  {
    printf("conc. %s is %g\n", (pelements[seq_elms_in_layer[i]]).name,
                               (pelements[seq_elms_in_layer[i]]).wfract);
  }

  printf("-----------------------------------------------\n");
}


double attn_factor(int line_nr, int tot_lines, elem *pelem,
                  int nr_elem,  int *sequence)
{
   int i,  l ;
   double result, cj ;

   result=1.0;
   l=sequence[line_nr];
   for (i=0; i < nr_elem ; i++) 
   {
      cj= (pelem[i]).wfract;
      /* debugging 
      printf("in attn_factor line_nr is %d i is %d l is %d  cj is %f \n",line_nr, i,l,  cj); /* */

      if (i != l )
      {
        /* debugging
        printf("Sum term alpha_ij %d %d %g  cj %g\n",
                 line_nr, i, alpha_ij[line_nr][i], cj); /* */
        result += alpha_ij[line_nr][i]*cj;
      }
   }
   return(result);
}
  
     
double enh_factor(int line_nr, int tot_lines, elem *pelem,
                  int nr_elem,  int *sequence)
{
   int i, l ;
   double result, cj ;

   result=1.0;
   l=sequence[line_nr];
   for (i=0; i < nr_elem ; i++)  
   {
      cj= (pelem[i]).wfract;
      if (i != l )
      {
        /* debugging
        printf("Sum term rho_ij %d %d %g  cj %g\n",
                 line_nr, i, rho_ij[line_nr][i], cj);  /* */
        result += rho_ij[line_nr][i]*cj;
      }
   }
   return(result);
}


