/***********************************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*                    rouss.c
*                    routine to calculate all binary
*                    and ternary interaction coeffs
*                    for LTQC procedure 
*                    calculations based on fundamental parameter
*                    approach
*                    M.Bos
*                    jan 1998
*
*************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "xrfluor.h"

extern elem_list lijst_els ;
extern meas_line  lines_to_meas[MAX_LINES] ;
extern double take_off, berryliumd; 
extern double muzlambda[MAX_LINES][TOTAL][PNTS_SPECTRUM];
extern char anode_el[3];
extern double **aij, **aijj , ***aijk;
extern struct spectrum_pair **pntrspectra;   
extern double *intbulk_calc, *intprep_calc ;

void calc_alphas(multilayer *sample, int nr_lines)
{
   int  nr_layers, tot_elements;
   layer **plagen, *plaag, *pbulk, **plagenbin, *plaagbin;
   layer **plagentert, *plaagtert;
   elem *pelem_in_laag;
   int z_elem_in_laag, z_meas, z_filt, z_derde_elem;
   int i,j, k,l, m ;
   multilayer *binsmpl1, *binsmpl2, *tertsmpl1, *tertsmpl2, *tertsmpl3 ;
   char eli[3], elj[3], elk[3], lin[6] ;
   multilayer bulksample;
   double totint_bulk, lambda_in, primint_bulk, mth_filt, kev;
   double primint, secint, totint;
   double rxilow, rxihigh , rximedium ;
   double  fi_2_8_0, fi_8_2_0, fi_3_35_35, fi_3_7_0, fi_3_0_7;

   tot_elements=0; 
   nr_layers=(*sample).nr_layers;
   plagen= (*sample).lagen;
   for (i=0; i < nr_layers; i++)
   {
    plaag=plagen[i];
    tot_elements += (*plaag).nr_elements;
   }
   /* debugging
   printf("Totaal aantal elementen in preparaat %d \n",tot_elements); /* */
   
   /* alloc space for influence coeffs */

   aij= (double **) malloc(nr_lines*sizeof(double * ));
   for (i=0; i < nr_lines ; i++)
      aij[i]= (double *) malloc(tot_elements*sizeof(double));
   aijj= (double **) malloc(nr_lines*sizeof(double * ));
   for (i=0; i < nr_lines ; i++)
      aijj[i]= (double *) malloc(tot_elements*sizeof(double));
    
   
   aijk= (double ***)malloc(nr_lines*sizeof(double **));
   for (i=0; i < nr_lines; i++)
   {
      aijk[i]= (double **)malloc(tot_elements*sizeof(double *));
      for (j=0; j < tot_elements ; j++)
       {
         aijk[i][j]= (double *)malloc(tot_elements*sizeof(double));
       }
   }
   
for (k=0; k < nr_lines; k++)
   {
    z_meas=(lines_to_meas[k]).z;
    strcpy(eli, get_asymb(z_meas, &lijst_els));
    strcpy(lin, (lines_to_meas[k]).line);
    z_filt=(lines_to_meas[k]).filt_z;
    mth_filt=(lines_to_meas[k]).mth_filt;
    kev=(lines_to_meas[k]).kev;
    pntrspectra[k]=gen_spec(anode_el, take_off, berryliumd, kev,
              PNTS_SPECTRUM, z_filt, mth_filt);
    pbulk=mkbulklay(eli);
    bulksample.nr_layers=1;
    bulksample.lagen=&pbulk;
    totint_bulk =0.0;
    for (i=0; i < PNTS_SPECTRUM; i++)
    {
      primint_bulk= pilabda(k,i, pbulk)*(pntrspectra[k][i]).intens ;
/*      printf("intens by %d is %g\n",i, (pntrspectra[k][i]).intens); */
      totint_bulk += primint_bulk;
    }
    intbulk_calc[k]=totint_bulk;

    printf("Bulk int %s %s %e \n", eli, lin, totint_bulk);

    
    for (i=0; i < nr_layers; i++)
    {
     plaag=plagen[i];
     pelem_in_laag= (*plaag).elementen;

     /* start calc of interactions  */
     for (j=0; j< (*plaag).nr_elements ; j ++)
     {
       z_elem_in_laag= get_z(((*plaag).elementen[j]).name, &lijst_els);
       if ( z_elem_in_laag != z_meas)   /* only for oher elems */
       {
         /* start with binary coefs */
         strcpy(elj, get_asymb(z_elem_in_laag, &lijst_els)); 
         binsmpl1=mkbinbulk(eli,elj, 0.2);
         plagenbin= (*binsmpl1).lagen;
         plaagbin=plagenbin[0];
         totint=0.0;
         for (l=0; l < PNTS_SPECTRUM ; l++)
         {  
           lambda_in =(pntrspectra[k][l]).lambda;
           primint= pilabda(k,l, plaagbin)*(pntrspectra[k][l]).intens; 
           secint=intrasecflu(eli,lin, lambda_in, binsmpl1, 0, primint,
                              k,l);
           totint += primint + secint;
         }
         intprep_calc[k]=totint;
         rxilow= intprep_calc[k]/intbulk_calc[k];
 /*        printf("rxi low voor lijn %s van %s met %s is %g %g\n",
                 lin, eli, elj, rxilow, intprep_calc[k]);  /* */
         binsmpl2=mkbinbulk(eli,elj, 0.8);
         plagenbin= (*binsmpl2).lagen;
         plaagbin=plagenbin[0];
         totint=0.0;
         for (l=0; l < PNTS_SPECTRUM ; l++)
         {  
           lambda_in =(pntrspectra[k][l]).lambda;
           primint= pilabda(k,l, plaagbin)*(pntrspectra[k][l]).intens; 
           secint=intrasecflu(eli,lin, lambda_in, binsmpl2, 0, primint,
                              k,l);
           totint += primint + secint;
         }
         intprep_calc[k]=totint;
         rxihigh= intprep_calc[k]/intbulk_calc[k];
/*         printf("rxi high voor lijn %s van %s met %s is %g %g\n",
                 lin, eli, elj,rxihigh, intprep_calc[k]); */
         fi_2_8_0= (0.2/rxilow - 1.0)/0.8;
         fi_8_2_0= (0.8/rxihigh - 1.0)/0.2;

         aij[k][j]=(-fi_2_8_0 + 4.0*fi_8_2_0)/3.0 ;
         aijj[k][j]= 5.0*(fi_2_8_0 - fi_8_2_0)/3.0;
         /* debugging
         printf("Binary interaction coeffs aii %s %s %g, aijj %g\n",
                 eli, elj, aij[k][j], aijj[k][j]); /* */

        }
         /* Now calculate ternary interactions */
      
         for (l=j; l < (*plaag).nr_elements ; l++)
         {
           z_derde_elem= get_z(((*plaag).elementen[l]).name, &lijst_els);
         if ((z_elem_in_laag != z_meas) && (z_elem_in_laag != z_derde_elem) && (z_meas != z_derde_elem))  
         {
           strcpy(elk, get_asymb(z_derde_elem, &lijst_els)); 
           tertsmpl1=mktertbulk(eli, elj, elk, 0.3, 0.35);
           plagentert= (*tertsmpl1).lagen ;
           plaagtert= plagentert[0];
           totint=0.0;
           for (m=0; m < PNTS_SPECTRUM;  m++)
           {
             lambda_in= (pntrspectra[k][m]).lambda;
             primint= pilabda(k, m, plaagtert)*(pntrspectra[k][m]).intens;
             secint=intrasecflu(eli,lin, lambda_in, tertsmpl1, 0, primint,
                                k, m);
             totint += primint+secint;
           }
           intprep_calc[k]=totint;
           rximedium=intprep_calc[k]/intbulk_calc[k];
           tertsmpl2=mktertbulk(eli, elj, elk, 0.3, 0.7);
           plagentert= (*tertsmpl2).lagen ;
           plaagtert= plagentert[0];
           totint=0.0;
           for (m=0; m < PNTS_SPECTRUM;  m++)
           {
             lambda_in= (pntrspectra[k][m]).lambda;
             primint= pilabda(k, m, plaagtert)*(pntrspectra[k][m]).intens;
             secint=intrasecflu(eli,lin, lambda_in, tertsmpl2, 0, primint,
                                k, m);
             totint += primint+secint;
           }
           intprep_calc[k]=totint;
           rxihigh=intprep_calc[k]/intbulk_calc[k];
           tertsmpl3=mktertbulk(eli, elj, elk, 0.3, 0.0);
           plagentert= (*tertsmpl3).lagen ;
           plaagtert= plagentert[0];
           totint=0.0;
           for (m=0; m < PNTS_SPECTRUM;  m++)
           {
             lambda_in= (pntrspectra[k][m]).lambda;
             primint= pilabda(k, m, plaagtert)*(pntrspectra[k][m]).intens;
             secint=intrasecflu(eli,lin, lambda_in, tertsmpl3, 0, primint,
                                k, m);
             totint += primint+secint;
           }
           intprep_calc[k]=totint;
           rxilow=intprep_calc[k]/intbulk_calc[k];
           fi_3_35_35= (0.3/rximedium -1.0)/0.35 ;
           fi_3_7_0= (0.3/rxihigh -1.0)/0.7 ;
           fi_3_0_7 = (0.3/rxilow -1.0)/0.7 ;   /* ???????? */
           
           aijk[k][j][l]=(20.0/7.0)*(fi_3_35_35 - fi_3_7_0 -fi_3_0_7);
           /* debugging 
           printf("Ternary influence coef analyte %s with %s and %s is %g \n",
                    eli, elj, elk, aijk[k][j][l]);  /* */
        }    
       }
     }  
    }
  }
  
}

