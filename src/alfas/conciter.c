/**********************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*     conciter.c
*
*     calculates sample composition 
*     from relative intensities
*     using influence coefficients  aij, aijj, aijk
*     and eqn (11) from
*     R.M.Rouseau , X-ray Spectrom. 13 (1984) 121
*
*     Input parameter is pointer to a multilayer sample
*     On output this structure contains the result 
*     For now only top layer is inspected and used
*
*     M.Bos
*     jan 1998
*
**********************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern int nr_of_meas_lines ;
extern multicomplayer compsample ;
extern struct spectrum_pair **pntrspectra;
extern meas_line lines_to_meas[MAX_LINES];
extern elem_list lijst_els;
extern double *rxi_meas;
extern double *rxi_calc;
extern double *intprep_calc;
extern double *intbulk_calc;
extern int nr_free_vars;
extern char anode_el[4];
extern double take_off, berryliumd ;
extern double **aij;
extern double **aijj;
extern double ***aijk;



void conciter(multilayer *sample)
{
  int i,j ;
  layer *ptop ;
  int nr_elem ;
  double delta;
  elem *pelements ;
  char meas_elem[3] ;
  int z_meas_elem ;
  int *seq_elms_in_layer ;
  double newconc;

  ptop=(*sample).lagen[0];
  nr_elem= (*ptop).nr_elements ;
  pelements= (*ptop).elementen ;
  seq_elms_in_layer= (int *) malloc(nr_of_meas_lines*sizeof(int));

  /* set up correspondence of elements to measured line */

  for (i=0; i < nr_of_meas_lines ; i++)
  {
     z_meas_elem= (lines_to_meas[i]).z ;
     strcpy(meas_elem, get_asymb(z_meas_elem, &lijst_els));
     j=0;
     while (( strcmp(meas_elem, (pelements[j]).name) !=0 )&& (j <= nr_elem)) 
     {
        /* debugging
        printf("%s \n", (pelements[j]).name); /* */
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
      newconc=rxi_meas[i]*inflfactor(i, nr_of_meas_lines,pelements,
                                      nr_elem, seq_elms_in_layer);
      delta+= fabs(newconc- (pelements[seq_elms_in_layer[i]]).wfract);
      (pelements[seq_elms_in_layer[i]]).wfract=newconc;
      /* debugging
      printf("newconc %d %g \n", i, newconc); /* */
    }
  } while (delta > 1.0E-10) ;
  printf("----------------------------------------------\n");
  printf("Estimated conc. after first step rousseau algorithm\n");
  for (i=0; i < nr_of_meas_lines;i++)
  {
    printf("conc. %s is %g\n", (pelements[seq_elms_in_layer[i]]).name,
                               (pelements[seq_elms_in_layer[i]]).wfract);
  }
  printf("-----------------------------------------------\n");
}


double inflfactor(int line_nr, int tot_lines, elem *pelem,
                  int nr_elem,  int *sequence)
{
   int i, j, l;
   double result, cj, ck, ci ;

   result=1.0;
   l=sequence[line_nr];
   ci= (pelem[l]).wfract;
   for (i=0; i < nr_elem ; i++)  /* calc Cj(aij+ aijjCM) */
   {
      cj= (pelem[i]).wfract;
      if (i != sequence[line_nr])
      {
        /* debugging
        printf("Sum term ij %d %d %g ci %g cj %g\n",
                 line_nr, i, aij[line_nr][i], ci, cj); /* */
        result += (aij[line_nr][i] +aijj[line_nr][i]*(1.0 - ci))*cj;
      }
   }
 
   for (i=0; i < nr_elem; i++)
   {
      cj= (pelem[i]).wfract;
      if ( i != sequence[line_nr])
      {
        for (j=0; j < nr_elem ; j++)
        if (( j != sequence[line_nr]) && ( j!= i))
        {
          ck=(pelem[j]).wfract;
          /* debugging
          printf("Sum term ijk %d %d %d %g cj %g ck %g\n",
                  line_nr, i, j, aijk[line_nr][i][j], cj, ck); /* */
          result += aijk[line_nr][i][j]*cj*ck;
        }
      }
    }

    return(result);
}
  
     
