/*********************************************************
*
*   File conprof.c
*
*   Calculates composition of 50 -layer sample
*   with exponentially decreasing Ag2O molfraction
*
*   input data from stdin : 
*   nr of compounds   value of exponent for decay of conc.
*   per compound  Compound Formula, Molwt, density and mol fraction
*   OF COMPOUNDS in  blank bulk glass
*   last item  should be the compound exchanged with Ag2O 
*   (in this case Na2O)
*  
*   M.Bos
*   May 2000
*
***********************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NRLAYERS 10
#define THICKNESS 2.50e-3    /* ANALYZING DEPTH in cm */

int main(void)
{
    int i, j, k, nr_compounds, retvalscanf ;
    double decay_exp;
    char cmpdnames[100] [100] ;  /* far too big, who cares? */
    double molwt[100], wtfraction[100], density[100], molfrac[100] ;
    double weight[100], weightf_layers[NRLAYERS][100] , mean_density[NRLAYERS] ;
    double moles_tot, molfag2o, massthickn[NRLAYERS];
   
    retvalscanf= scanf("%d", &nr_compounds);
    if (retvalscanf==EOF)
    {
      printf("error in input, no nr of compounds found\n");
      exit(1);
    }

    retvalscanf= scanf("%lf", &decay_exp);
    if (retvalscanf==EOF)
    {
      printf("error in input, no decay exponent\n");
      exit(1);
    }
    /* now get data per compound */

    for (i=0; i < nr_compounds; i++)
    {
       
       retvalscanf=scanf("%s", cmpdnames[i]);
       if (retvalscanf==EOF)
       {
          printf("error in input file no compoundname found\n");
          exit(2);
       }
       retvalscanf=scanf("%lf", &(molwt[i]));
       if (retvalscanf==EOF)
       {
          printf("error in input file no molwt found\n");
          exit(2);
       }
       retvalscanf=scanf("%lf", &(density[i]));
       if (retvalscanf==EOF)
       {
          printf("error in input file no density found\n");
          exit(2);
       }
       retvalscanf=scanf("%lf", &wtfraction[i]);
       if (retvalscanf==EOF)
       {
          printf("error in input file no wtfraction found\n");
          exit(2);
       }
       /* debugging 
       printf("wtfracion[%d] %e ", i, wtfraction[i]); /* */
    }
     molwt[nr_compounds] = 231.736  ;  /* molwt Ag2O */
     strcpy(cmpdnames[nr_compounds],"Ag2O");
     density[nr_compounds]= 7.14 ;   /* density Ag2O */

   /* now start to setup outputfile */

   printf("%d\n",NRLAYERS);              /* 5 layers */

     /* now calculate moles of all compounds in  sample */

     moles_tot = 0.0;
     for (j=0; j < nr_compounds; j++)
     {
       moles_tot += (wtfraction[j])/molwt[j] ;

       /* debugging 
       printf("moles tot %e wtfraction[%d] %e molwt[%d]  %e\n",
                 moles_tot, j, wtfraction[j],j,  molwt[j]); /* */ 
     }

     /* molfractions now */

     for (j=0; j < nr_compounds; j++)
     {
       molfrac[j]= (wtfraction[j])/(molwt[j]*moles_tot);

       /* debugging 
       printf("molfrac[%d] is %e\n", j, molfrac[j]); /* */
     }


 
     molfag2o= molfrac[nr_compounds-1];
   /* produce data per layer */

   for (i=0; i < NRLAYERS ; i++)
   {
     printf("%d ", nr_compounds+1); /* now including Ag2O */

     /* first calculate Ag2O via exponential decay */
     molfrac[nr_compounds]= molfag2o*exp(-decay_exp*
                         (THICKNESS/(2.0*NRLAYERS) + 
                i*(THICKNESS/(NRLAYERS))));
     molfrac[nr_compounds-1]=molfag2o-molfrac[nr_compounds];

     /* now calculate total weight of layer */
     weight[i] = 0.0;

     for (j=0; j <= nr_compounds  ; j++)
     {
        weight[i] += molfrac[j]* molwt[j];
     }
     /* debugging 
     printf("weight[%d] is %e\n",i,weight[i]); /* */
     mean_density[i] = 0.0 ;
     /* now fill data weightfractions all layers */
     for (j=0; j <= nr_compounds; j++)
     {
       weightf_layers[i][j]=molfrac[j]*molwt[j]/weight[i];
       mean_density[i] += weightf_layers[i][j]*density[j];
     }
     massthickn[i]=(THICKNESS/NRLAYERS)*mean_density[i] ; /* per slice xx mm */

     printf(" %e ", massthickn[i]);
     /* now print data */
     for (j=0; j <= nr_compounds; j++)
     {
       printf("%s %le ", cmpdnames[j], molfrac[j]*molwt[j]/weight[i]);
     }
     printf("\n");
  }
  exit(0);
}

