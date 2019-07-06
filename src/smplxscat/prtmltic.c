/*
*
*    COPYRIGHT (C) M.BOS 1998 , for details see copying
*
*
*      file prtmltic.c
*      prints composition of multilayer defined as
*      compounds to given file
*
*      M.Bos
*      june 1998
*
************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xrfluor.h"

extern elem_list lijst_els;

int prtmcomp(FILE *fpout, multicomplayer *multicomplay)
{
    int i, j;
    int nr_layers, nr_compounds;
    double wfract, mth ;
    complayer *compslay;
    char *comp_name ;
    nr_layers=get_nr_clay(multicomplay);
    for (i=0; i < nr_layers ; i++)
    {
        compslay = ((*multicomplay).lagen[i]);
        nr_compounds=(*compslay).nr_compounds;
        mth=(*compslay).massthickness ;
        fprintf(fpout,
        "Layer %d has %d compounds and a massthickness of %g  g/cm2\n",
         i, nr_compounds, mth);
        for (j=0; j < nr_compounds; j++)
        {
          comp_name=(*compslay).comps[j].name;
          wfract=(*compslay).comps[j].wfract;
          fprintf(fpout,"Compound %s wfract %g  \n", 
                 comp_name, wfract);
        }
    
     }
    return(nr_compounds);
}
