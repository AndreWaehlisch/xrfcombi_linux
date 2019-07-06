/*****************************************************
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*	compdata.c
*	datstructure and access routines
*	for multilayer thin film samples
*	consisting of chemically defined 
* 	compounds
*
*	M.Bos
*	october 1997
*
*******************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include "xrfluor.h"

extern FILE *own_fopen(char *);

/********************************************************
*  
* routine to build structure for sample data
* in terms of compounds 
*
*  inputs:
*  (1)  file pointer to file with
*  ascii representation of structure:
*  number of layers
*  per layer:
*  number of compounds, massthickness, compdata
*  compdata: compound symbol, weight fraction
*  (2) pointer to output variable of type multicomplayer
*
*  datstructures are mallocated from memory
*  and elements in the multicomplayer variable
*  are filled with adresses or data
*
********************************************************/

void bldcompos(FILE *fp, multicomplayer *alles)
{
	int nr_layers ;
   	int nr_comps ;
        int i,j, k, nr_atoms, cnt_at ;
        int fix_thickn, fix_wfract ;
        char *symbol_1;
        elem_cnt *el_pairs ;
        char regel[LINELENGTH];
        char *position;
        FILE *fpcomplist;
  	char sym_compound[MAX_CHAR_COMP] ;
	double weight_fraction, thickness;
        complayer **setlayers ;
        compound *reeks ;
        extern int nr_free_vars ;

        nr_free_vars =0 ; 
        fscanf(fp, "%d", &nr_layers);
        /* debugging
        printf("nr of layers with compounds is %d\n",nr_layers); /* */
        setlayers = (complayer * * )malloc(nr_layers*sizeof(complayer));
        for (i=0; i < nr_layers; i++)
        {
  	  fscanf(fp,"%d",&nr_comps);
          fscanf(fp,"%lf",&thickness);
          fscanf(fp,"%d",&fix_thickn);
          if (fix_thickn ==0 ) nr_free_vars++ ;
          /* debugging 
          printf("nr %d thickn %g fix %d \n",nr_comps, thickness,fix_thickn);/* */
	  reeks=(compound *)malloc(nr_comps*sizeof(compound));
          for (j=0; j < nr_comps; j++)
          {
             fpcomplist=own_fopen("./data/complist.dat");
             fscanf(fp,"%s %lf %d ",sym_compound, &weight_fraction,
                    &fix_wfract);
             if (fix_wfract==0) nr_free_vars++;
             /* debugging 
             printf("gelezen is %s %g %d\n",sym_compound, weight_fraction,
                    fix_wfract);  /* */
             reeks[j].wfract=weight_fraction;
             reeks[j].wfixed=fix_wfract;
             strcpy(reeks[j].name,sym_compound);

             /* now read from file complist to match compounds name */
 
             while ( (position=fgets(regel,LINELENGTH,fpcomplist)) != NULL)
             { 
               symbol_1=strtok(regel," ");   /* to symbol */
               if ( strcmp(symbol_1,sym_compound) ==0)
                break;
             }

             /* we now have the matching line or are past the end */

             if (position==NULL)
             {
               printf("unknown compound %s \n",sym_compound);
               exit(1);
             }

             position=strtok(NULL, " "); /* to density , throw away*/
             position=strtok(NULL," ");  /* AT or WT ? */
             if (strcmp(position, "AT") != 0)
             {
               printf("In ./data/complist.dat only AT entries for now\n"); 
               exit(1);
             }
             position=strtok(NULL," ");  /* to number of atoms */
             nr_atoms=atoi(position);
             el_pairs=(elem_cnt *) malloc(nr_atoms*sizeof(elem_cnt));
             reeks[j].nr_elems=nr_atoms;
             reeks[j].elementen=el_pairs;

             for (k=0; k < nr_atoms; k++)
             {
               position=strtok(NULL," ");   /* symbol of element */
               /* convert to lower case to match rest of system */
/*               if ((*(position+1) < 0x64) && ((*(position+1) != 0x20)) 
               *(position+1)= *(position+1)+0x20; */
               strcpy(el_pairs[k].name,position); /* name of element in place */
               position=strtok(NULL," \n");  /* atom count */
               cnt_at=atoi(position);
               el_pairs[k].count=cnt_at;
             }
           }
           setlayers[i]= (complayer *) malloc(sizeof(complayer));
           setlayers[i]->nr_compounds = nr_comps;
           setlayers[i]->comps=reeks;
           setlayers[i]->massthickness= thickness;
           setlayers[i]->mthfixed=fix_thickn;
        }
        (*alles).nr_layers=nr_layers;
        (*alles).lagen=setlayers;
}
/*******
*  everything is now filled 
* here follow the access routines
*/

int get_nr_clay(multicomplayer *mclayer)
{
	return( (*mclayer).nr_layers);
}

int get_nr_comps(complayer *clayer)
{
    return(*clayer).nr_compounds;
}

complayer *get_clayer(multicomplayer *mclayer, int layer_nr)
{
      complayer *layer;
      layer = (*mclayer).lagen[layer_nr];
      return(layer);
}

char *get_comp_name(complayer *layer, int comp_nr)
{
     return ( (*layer).comps[comp_nr].name);
}

double get_cweightf(complayer *layer, int comp_nr)
{
    return( (*layer).comps[comp_nr].wfract);
}

int get_fixwght(complayer *layer, int comp_nr)
{
    return( (*layer).comps[comp_nr].wfixed);
}

double get_cmth(complayer *layer)
{
	return( (*layer).massthickness);
}

int get_fixmth(complayer *layer)
{
	return( (*layer).mthfixed);
}


int get_params(multicomplayer *mclayer, double *free_params)
{
   int i, j, nr_layers, cnt_free_par, nr_comps ;
   complayer *laag;

   cnt_free_par=0;
   nr_layers=get_nr_clay(mclayer);
   /* debugging 
   printf("nr_layers in get params is %d \n",nr_layers); /* */
   for (i=0; i < nr_layers ; i++)
   {
     laag=get_clayer(mclayer, i);
     if (get_fixmth(laag)==0) 
      {
        free_params[cnt_free_par]=get_cmth(laag);
        /* debugging 
        printf("fixmth is 0 \n"); /* */
        cnt_free_par++;
      }
     nr_comps=get_nr_comps(laag);
     /* debugging 
     printf("nr_comps in get_params is %d\n",nr_comps); */
     for (j=0; j < nr_comps ; j++)
     {
       if (get_fixwght(laag, j)==0)
        {
          free_params[cnt_free_par]=get_cweightf(laag,j);
          cnt_free_par++;
        }
     } 
  }  
  return(cnt_free_par);
}
    

int put_params(multicomplayer *mclayer, double *free_params)
{
   int i, j, nr_layers, cnt_free_par, nr_comps ;
   complayer *laag;

   cnt_free_par=0;
   nr_layers=get_nr_clay(mclayer);
   for (i=0; i < nr_layers ; i++)
   {
     laag=get_clayer(mclayer, i);
     if (get_fixmth(laag)==0) 
      {
        (*laag).massthickness= free_params[cnt_free_par];
        cnt_free_par++;
      }
     nr_comps=get_nr_comps(laag);
     for (j=0; j < nr_comps ; j++)
     {
       if (get_fixwght(laag, j)==0)
        {
           (*laag).comps[j].wfract=free_params[cnt_free_par];
          cnt_free_par++;
        }
     }
  }
  return(cnt_free_par);
}

void show_params(double *free_pars)
{
  int i;
  extern int nr_free_vars ;
  for (i=0; i < nr_free_vars ; i++)
  {
    printf(" Par %d %g \n", i, free_pars[i]);
  }
   printf("Respons is %g \n", free_pars[nr_free_vars]);
}
  


/***************************
*
*  checkcomps()
*
*  for debugging purposes
* checks is access routines work
*  and that structures contain the required data
*
*************************************************/

void checkcomps(multicomplayer * mclayer)
{
	complayer *player ;
        int i, j ;
	int aant_comps ;
        int cnt_lay;
        cnt_lay=get_nr_clay(mclayer);
        printf("Number of layers with compounds is %d \n",cnt_lay);
        for (i=0; i < cnt_lay; i++) 
        {
	  player=get_clayer(mclayer, i);
          aant_comps=get_nr_comps(player);
          printf("Thickness layer %d (fix %d) is %g with %d compounds\n",
                  i, get_fixmth(player), get_cmth(player),aant_comps);
          for (j=0; j< aant_comps; j++)
          {
            printf("Compounds %s - fraction %g (fix %d) \n",
                    get_comp_name(player,j), get_cweightf(player,j),
                    get_fixwght(player,j));
           }
         }
}

 
