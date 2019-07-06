/**********************
*
*        COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*        jmpedg.c
*
*        routine to determine
*        absorption jump factors
*
*        for the moment based 
*        on tables from
*      
*        E.P. Bertin, "principles of X-ray Spectrometric Analysis"
*        2nd edn, Plenum New York 1975
*
*        later to be adapted to
*        D.K.G. de Boer
*        Spectrochim. Acta, 44B (1989) 1171
*
*        M.Bos
*        August 1997
*
*********************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"

extern double edges[TOTAL][10];
extern double nkl123[TOTAL][3];
extern double fij[TOTAL][3];
extern double omegas[TOTAL][NR_OMEGAS];

void bldjmpedg(FILE *fpjmpdata,double jmpdata[][2])
{
  char regel[LINELENGTH];
  double kedge_r, L3edge_r ;
  int  j ;
  char *positie;
  j=3;
  while (fgets(regel,LINELENGTH,fpjmpdata) != NULL)
  {
    positie=strtok(regel, " ");  /* elem symbol */
    /* debugging
    printf("%s \n",positie); /* */
    positie=strtok(NULL, " "); /* now r of kedge */
    /* debugging
    printf("%s \n",positie); /* */
    sscanf(positie,"%lf",&kedge_r);
    positie=strtok(NULL," ");
    /* debugging   
    printf("%s \n",positie); /* */
    positie=strtok(NULL," ");  
    /* debugging
    printf("%s \n",positie); /* */ 
    if ( *positie=='-')
       L3edge_r=1.0;
    else
       sscanf(positie,"%lf",&L3edge_r);
    jmpdata[j][0]=kedge_r;
    jmpdata[j][1]=L3edge_r;
    j++ ;
  }
}

/*******************************************
*
*       routine to determine jumpfactors
*
*       for creation of specific hole in atom
*
*       holes are numbered:
*       1   K 
*       2   LI
*       3   LII
*       4   LIII
*       5   M2 
*       6   M3
*       7   M4
*       8   M5
*
************************************************/

double getjmpfact(int z, int hole_nr, double energy, double jmpdata[][2])
{
   double jmpratio,jmpfactor;
   double rl1, rl2, rl3, rk;
   double rm1, rm2, rm3, rm4, rm5 ;
   double omegam2, omegam3, omegam4;

   rl1=1.2 ;
   rl2=1.4 ;
   if (( z < 4 ) || (z >94))
   {
     printf("element number out of range in call to getjmpfact\n");
     exit (1);
   }
  
   if (hole_nr > 8)
   {
    printf("Creation of holes other than K or L123 M2345 not implemented\n");  
     return(0.0);
   }

   if (hole_nr==1)
   {
       jmpratio=jmpdata[z-1][0];
       jmpfactor=(jmpratio-1.0)/jmpratio; 
       return(jmpfactor);
   }
   if (hole_nr==2)     /* L1 */
   {
      if (energy < edges[z][1])
        return(0.0);
      if (energy < edges[z][0])
        return( (rl1- 1.0)/rl1);
      else
       {
         rk= jmpdata[z-1][0];
         return( 0.16667/rk + (rk-1.0)*nkl123[z][0]/rk); 
       }
    }

   if (hole_nr==3)    /* L2 */
   {
     if (energy < edges[z][2])
       return(0.0);
     if (energy < edges[z][1])
       {
         return( (rl2-1.0)/rl2);
       }
     if (energy < edges[z][0])
       {
           jmpfactor= 0.2857142/rl1;
           jmpfactor+=0.1666666666667*fij[z][0];
           return(jmpfactor);
       }
     else                              /* E > Ek */
       {
         jmpfactor= 0.2857142/rl1;
         jmpfactor+=0.1666666666667*fij[z][0];
         rk=jmpdata[z-1][0];
         jmpfactor /= rk; 
         jmpfactor += ((rk-1.0)/rk)*(nkl123[z][0]*fij[z][0]+ nkl123[z][1]);
         return(jmpfactor); 
       }
      
    }

   if (hole_nr==4)  /* L3 hole */
   {
      if (energy < edges[z][3])
        return(0.0);
      if (energy < edges[z][2])
       { 
         rl3=jmpdata[z-1][1];
         /* printf("rl3 is %g \n",rl3); */
         jmpfactor=(rl3-1.0)/rl3;
     /*     printf("rl3 sec jmpfactor %g \n",jmpfactor);  */
         return(jmpfactor); 
       }
      if (energy < edges[z][1])
       {
         rl3=jmpdata[z-1][1];
       /*  printf("t\rl3 is %g \n",rl3); */
         jmpfactor=(rl3-1.0)/rl3;
         jmpfactor /= rl2;
         jmpfactor += ((rl2-1.0)/rl2)*fij[z][2];
        /*  printf("rl3 ook rl2 jmpfactor is %g \n",jmpfactor);  */
         return(jmpfactor);
       }
      if (energy < edges[z][0])
       {
         /* printf("hole 4 (L3) edges %g energy %g\n",edges[z][0],energy); */
         rl3=jmpdata[z-1][1];
    /*      printf("rl3 is %g \n", rl3); */
         jmpfactor=(rl3-1.0)/rl3;
         jmpfactor /= rl2;
         jmpfactor += ((rl2-1.0)/rl2)*fij[z][2];
         jmpfactor /= rl1;
         /* printf("Z %d rl1 %g rl2 %g rl3 %g f12 %g f13 %g f23 %g\n",
                  z, rl1, rl2, rl3, fij[z][0],fij[z][1],fij[z][2]); */
         jmpfactor += ((rl1-1.0)/rl1)*(fij[z][1]+fij[z][0]*fij[z][2]);
      /*     printf("jmpfactor %g met rl3 rl2 en rl1\n",jmpfactor);  */
         return(jmpfactor);
       }
      else
       {
         rl3=jmpdata[z-1][1];
       /*   printf("rl3 is %g \n", rl3); */
         jmpfactor=(rl3-1.0)/rl3;
         jmpfactor /= rl2;
         jmpfactor += ((rl2-1.0)/rl2)*fij[z][2];
         jmpfactor /= rl1;
         jmpfactor += ((rl1-1.0)/rl1)*(fij[z][1]+fij[z][0]*fij[z][2]);
         rk=jmpdata[z-1][0];
         jmpfactor /=  rk;
         jmpfactor += ((rk-1.0)/rk)*(nkl123[z][2]+nkl123[z][1]*fij[z][2]+
                                      nkl123[z][0]*fij[z][1]);
     /*     printf("Meest complexe jmpfact %g \n",jmpfactor);  */
         return(jmpfactor);
       }
   }

   rm1=1.1;
   rm2=1.1 ;
   rm3=1.2;
   rm4=1.5;
   rm5= 225.0/z - 0.35 ;
   omegam2=get_omega(z, "m2", omegas);
   omegam3=get_omega(z, "m3", omegas);
   omegam4=get_omega(z, "m4", omegas);
  
   if (hole_nr==8)  /* M5 */
   {
      if (energy < edges[z][8])  
        return(0.0);           /* E < EM5 */
      if (energy < edges[z][7])    
        {
          jmpfactor=(rm5-1.0)/rm5 ;   /*   EM4 < E < EM5 */
          return(jmpfactor);
        }
      if (energy < edges[z][6])     /*  EM3 < E < EM4 */
       {
         jmpfactor=(rm5-1.0)/rm5;
         jmpfactor /= rm4;
         jmpfactor+= ((rm4-1.0)/rm4)*fij[z][12];  /* JM4*f45 */
         return(jmpfactor);
       }
      if (energy < edges[z][5])   /* EM3 < E < EM2 */
       {
         jmpfactor= (rm5-1.0)/(rm3*rm4*rm5);
         jmpfactor += ((rm3-1.0)/rm3)*(fij[z][11]+fij[z][10]*fij[z][12]); 
         return(jmpfactor);
       }
      if (energy < edges[z][4])    /* EM2 < E < EM1 */
       {
         jmpfactor= (rm5-1.0)/(rm2*rm3*rm4*rm5);
         jmpfactor += ((rm2-1.0)/rm2)*(fij[z][9]+fij[z][8]*fij[z][12]+
                        fij[z][7]*fij[z][10]*fij[z][12]);
         return(jmpfactor);
       }
      if (energy < edges[z][3])  /* EM1 < E < EL3 */
       {
          jmpfactor= (rm5-1.0)/(rm1*rm2*rm3*rm4*rm5);
          jmpfactor += ((rm1-1.0)/rm1)*(fij[z][6]+ fij[z][5]* fij[z][12]+
                       fij[z][4]*(fij[z][11]+fij[z][10]*fij[z][12])+
             fij[z][3]*(fij[z][9]+fij[z][8]*fij[z][12]+
              fij[z][7]*(fij[z][11]+fij[z][10]*fij[z][12])));
          return(jmpfactor);
       }
       if (energy < edges[z][2]) /* EL3 < E < EL2 */
       {
         if ( z < 57) return(0.0);
         if (z  < 79) return(1.01);
         else return(0.92);
       }
       if (energy < edges[z][1])  /* EL2 < E < EL1 */
       {
           if (z < 57) return(0.0);
           if ( z < 79 ) return(0.86);
            else return(0.68);
       }
       if (energy < edges[z][0]) /* EL1 < E < EK */
       {    
          if (z < 57) return(0.0);
          if (z < 79) return(0.92);
          else return(0.83);
        }
        
       if (energy >= edges[z][0]) 
         return(0.24);
  }     
   
   if (hole_nr==7)   /* M4 */
   {
     if (energy < edges[z][7])    /* E < EM4 */
         return(0.0);
     if (energy < edges[z][6] )   /* EM4 < E < EM3 */
       {
          jmpfactor=(rm4-1.0)/rm4 ;
          return(jmpfactor);
       }
     if (energy < edges[z][5])    /* EM3 < E < EM2 */
       {
         jmpfactor=(rm4 -1.0)/rm4;
         jmpfactor /= rm3 ;
         jmpfactor += ((rm3-1.0)/rm3)*fij[z][10] ;
         return(jmpfactor);
       }
     if (energy < edges[z][4])  /* EM2 < E < EM1 */
       {
           jmpfactor = (rm4-1.0)/(rm2*rm3*rm4);
           jmpfactor += ((rm2-1.0)/rm2)*(fij[z][8]+fij[z][7]*fij[z][10]);
           return(jmpfactor);
       }
     if (energy < edges[z][3])   /* EM1 < E < EL3 */
       {
           jmpfactor= (rm4-1.0)/(rm4*rm3*rm2*rm1);
           jmpfactor+= ((rm1-1.0)/rm1)*(fij[z][5]+ fij[z][4]*fij[z][10]+
                       fij[z][3]*(fij[z][8]+fij[z][7]*fij[z][10]));
           return(jmpfactor);
       }

     if (energy < edges[z][2])  /* EL3 < E < EL2 */
     {
        jmpfactor=0.20/(omegam2+omegam3+omegam4);
        jmpfactor *= omegam4;
        return(jmpfactor);
     }
     if (energy < edges[z][1]) /* EL2 < E < EL1 */
     {  
       jmpfactor=0.54/(omegam2+omegam3+omegam4);
       jmpfactor *=omegam4;
       return(jmpfactor);
     }
     if (energy < edges[z][0]) /* EL1 < E < EK */
     {
       jmpfactor=0.63/(omegam2+omegam3+omegam4);
       jmpfactor *=omegam4;
       return(jmpfactor);
     }
     if (energy >= edges[z][0])
     {
       jmpfactor=0.20/(omegam2+omegam3+omegam4);
       jmpfactor *=omegam4;
       return(jmpfactor);
     }
   }  

   if (hole_nr==5) /* M2 */
   {
     if ( energy < edges[z][5])   /*   E < EM2 */
      return(0.0);
     if (energy < edges[z][4])    /* EM2 < E < EM1 */
     {  
        jmpfactor= (rm2-1.0)/rm2 ;
        return (jmpfactor);
     }
     if (energy < edges[z][3])    /* EM1 < E < EL3 */
     {
       jmpfactor= (rm2-1.0)/(rm1*rm2);
       jmpfactor+= ((rm1-1.0)/rm1)*fij[z][3];
       return(jmpfactor);
     }
     if (energy < edges[z][2]) /*  EL3 < E < EL2 */
     {
        jmpfactor= 0.20 /(omegam2+omegam3+omegam4);
        jmpfactor *= omegam2;
        return(jmpfactor);
     }
     if (energy < edges[z][1]) /* EL2 < E < EL1 */
     {  
       jmpfactor=0.54/(omegam2+omegam3+omegam4);
       jmpfactor *=omegam2;
       return(jmpfactor);
     }
     if (energy < edges[z][0])  /* EL1 < E < EK */
      {
        jmpfactor= 0.63 /(omegam2+omegam3+omegam4);
        jmpfactor *= omegam2;
        return(jmpfactor);
      }
    if (energy >= edges[z][0] )  /* E > EK */
     {
        jmpfactor= 0.20 /(omegam2+omegam3+omegam4);
        jmpfactor *= omegam2;
        return(jmpfactor);
     }
   printf("Problems !!\n");
   }
   if (hole_nr==6)  /* M3  */
   {
      if (energy < edges[z][6])   /* E < EM3 */
       return(0.0);

      if ( energy < edges[z][5])  /* EM3 < E < EM2 */
      {
        jmpfactor = (rm3-1.0)/rm3;
        return(jmpfactor);
      }
      if ( energy < edges[z][4]) /* EM2 < E < EM1 */
      {
        jmpfactor = (rm3-1.0)/(rm3*rm2);
        jmpfactor += ((rm2-1.0)/rm2)*fij[z][7];
        return(jmpfactor);
      }
      if (energy < edges[z][3])  /* EM1 < E < EL3 */
      {
          jmpfactor= (rm3-1.0)/(rm1*rm2*rm3);
          jmpfactor += ((rm1-1.0)/rm1)*(fij[z][4] * fij[z][3]*fij[z][7]);
          return(jmpfactor);
      }
      if (energy < edges[z][2])  /* EL3 < E < EL2 */
      {
         jmpfactor= 0.20/(omegam2+omegam3+omegam4);
         jmpfactor *= omegam3;
         return(jmpfactor);
      }
      if (energy < edges[z][1]) /* EL2 < E < EL1 */
      {
         jmpfactor= 0.54/(omegam2+omegam3+omegam4);
         jmpfactor *= omegam3 ;
         return(jmpfactor);
      }
      if (energy < edges[z][0]) /* EL1 < E < EK */
      {
         jmpfactor= 0.63/(omegam2+omegam3+omegam4);
         jmpfactor *= omegam3 ;
         return(jmpfactor);
      }
      if (energy >= edges[z][0]) /*  E >= EK */
      {
         jmpfactor= 0.20/(omegam2+omegam3+omegam4);
         jmpfactor *= omegam3 ;
         return(jmpfactor);
      }
   }
}
