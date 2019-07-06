/*********************************************
*
*   COPYRIGHT 1999  (C) M.BOS , for details see file copying
*	
*   File bldpis.c
*
*   Program to build datastructure
*   for characteristic fluorescence lines
*   containing all measurement condition data
*   plus factors Pi cf Wang, X-Ray Spectrom. 25 (1996) 245
*
*   M.Bos
*   jan 1999
*
*  layout file pifcts.dat
*
*  Element, line name, tube_el, kV, mA, take_off, tube_beryllium,
   filt_elem, filt_mth, XTAL , collim, counter
*  pulseHigh, pulseLow,  factor Pi
*
*  collim = C(oarse) or F(ine)
*  Counter = F(low), S(cintillation), M(ixed)
*
******************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xrfluor.h"

extern struct line_details *el_details[TOTAL] ;

void bld_setpis(FILE *fp)
{

    int i, j, k;
    char *position ;
    char regel[LINELENGTH];
    int z;
    char el_name[3];
    char tube_el[3];
    char line_symb[6];
    char filt_el[5];
    char xtal_name[7];
    double filt_mth;
    double tube_volts, tube_millis , tube_takeoff, tube_beryllium ;
    char collim ;
    char counter;
    double pulslo, pulshi, pi;
    struct line_details *plindtls;
    struct tube_details *ptubedtls ;
    struct filter_details *pfiltdtls;
    struct xtal_details *pxtaldtls;
    struct collim_details *pcollimdtls ;
    struct counter_details *pcounterdtls;
    

    /* initialize pointer array for elements NULL */
    for (i=0; i < TOTAL; i++) el_details[i]=NULL;
   
   while ( fgets(regel,LINELENGTH,fp) != NULL)
   {
       /* debugging 
       fputs(regel, stdout);  /* */
       position=strtok(regel," ");
       strcpy(el_name, position);
       /* debugging 
       printf("Element %s \n",position); /* */
       z=get_z(position, &lijst_els); 
      
       /* now get pointer to node of element with z as atomic number */
      plindtls=  el_details[z];
      if (plindtls == NULL)
      {
        /* element not yet in dataset */
        /* debugging 
        printf("element %s not yet in dataset \n", position);
        /* now build datastructure for it */
        /* alles hagelnieuw ! */
        plindtls= ( struct line_details *) malloc(sizeof(struct line_details));
        if (plindtls== NULL)
        {
           printf("memory exhausted \n");
           exit(1);
        }
        el_details[z]=plindtls;
        /* get line_name  */
        position=strtok(NULL," ");
        /* debugging 
        printf("line symbol %s \n",position); /* */
        strcpy(line_symb, position);

        /* now put it in place */
        strcpy(plindtls->line_name, position);
        plindtls->next=NULL ;

        /* create tube details structure */
        ptubedtls= (struct tube_details *) malloc(sizeof(struct tube_details));
        plindtls->tube= ptubedtls;

        /* read the tube details now */        
        position=strtok(NULL," ");
        strcpy(tube_el,position);
        strcpy(ptubedtls->tube_el, position);
        position=strtok(NULL," ");
        sscanf(position,"%lf",&tube_volts);
        ptubedtls->voltage=tube_volts;
        position=strtok(NULL," ");
        sscanf(position,"%lf",&tube_millis);
        ptubedtls->current=tube_millis;
        position=strtok(NULL," ");
        sscanf(position,"%lf",&tube_takeoff);
        ptubedtls->take_off=tube_takeoff;
        position=strtok(NULL," ");
        sscanf(position,"%lf",&tube_beryllium);
        ptubedtls->d_beryllium=tube_beryllium;
        ptubedtls->next=NULL;
 
        /* debugging 
        printf("%s %s %s %f %f %f %f\n", el_name, line_symb, tube_el,
                tube_volts, tube_millis, tube_takeoff, tube_beryllium); /* */


        /* create filter data block and fill it */
        pfiltdtls= (struct filter_details *)malloc(sizeof(struct
                                        filter_details));
        if (pfiltdtls==NULL)
        {
           printf("memory exhausted in allocation filter details\n");
           exit(1);
        }     
        ptubedtls->filter=pfiltdtls;
        position=strtok(NULL," ");
        strcpy(filt_el,position); 
        position=strtok(NULL," ");
        sscanf(position,"%lf",&filt_mth);
        strcpy(pfiltdtls->filt_elem,filt_el);
        pfiltdtls->filt_mth=filt_mth;
        pfiltdtls->next=NULL;

        /* create xtal data block and fill it */
        pxtaldtls= (struct xtal_details *)malloc(sizeof(
               struct xtal_details));
        if (pxtaldtls==NULL)
        {
          printf("Memory exhausted in allocation xtal details\n");
          exit(1);
        }
        pfiltdtls->xtal=pxtaldtls;

        position=strtok(NULL," ");
        strcpy(xtal_name,position);
        strcpy(pxtaldtls->xtal_name,position);
        pxtaldtls->next = NULL;

        /* create collimator data block and fill it */
        pcollimdtls=(struct collim_details *) malloc(
                          sizeof(struct collim_details));

        if (pcollimdtls==NULL)
        {
          printf("Memory exhausted in alloc for collim_details\n");
          exit(1);
        }
        pxtaldtls->collimator=pcollimdtls;

        position=strtok(NULL," ");
        collim=position[0];
        pcollimdtls->type=position[0];
        pcollimdtls->next=NULL;

        /* create counter data block and fill it */
        pcounterdtls=(struct counter_details *)malloc(sizeof(
                      struct counter_details));
        if (pcounterdtls==NULL)
        {
          printf("Memory exhausted in alloc for counter details\n");
          exit(1);
        }
        position=strtok(NULL," ");
        counter=position[0];
        pcollimdtls->counter=pcounterdtls;
        position=strtok(NULL," ");
        sscanf(position,"%lf",&pulshi);
        position=strtok(NULL," ");
        sscanf(position,"%lf",&pulslo);
        position=strtok(NULL," \n");
        sscanf(position,"%lf",&pi);
        pcounterdtls->type=counter;
        pcounterdtls->pulslo=pulslo;
        pcounterdtls->pulshi=pulshi;
        pcounterdtls->pi=pi;
        pcounterdtls->next=NULL;
        
        /* debugging 
        printf("%s %s %s %6.2f %6.2f %6.2f %6.2f %s %6.2f"
               " %s %c %c %6.2f %6.2f %6.2e\n",
               el_name, line_symb, tube_el, tube_volts, tube_millis,
               tube_takeoff, tube_beryllium,
               filt_el, filt_mth, xtal_name, collim, counter,
               pulshi, pulslo, pi); /* */
     }

     else /* element was already present */
     {
        /* now look for the line_name */
        position=strtok(NULL," ");
        /* debugging 
        printf("line symbol %s \n",position); /* */
        strcpy(line_symb, position);

        /* now look for line_name entry in linked list line_details */

        while (strcmp(plindtls->line_name, line_symb))
        {
           if (plindtls->next != NULL)
            plindtls=plindtls->next;


        /* now put it in place */
        strcpy(plindtls->line_name, position);





        plindtls->next=NULL ;

        /* create tube details structure */
        ptubedtls= (struct tube_details *) malloc(sizeof(struct tube_details));
        plindtls->tube= ptubedtls;

        /* read the tube details now */        
        position=strtok(NULL," ");
        strcpy(tube_el,position);
        strcpy(ptubedtls->tube_el, position);
        position=strtok(NULL," ");
        sscanf(position,"%lf",&tube_volts);
        ptubedtls->voltage=tube_volts;
        position=strtok(NULL," ");
        sscanf(position,"%lf",&tube_millis);
        ptubedtls->current=tube_millis;
        position=strtok(NULL," ");
        sscanf(position,"%lf",&tube_takeoff);
        ptubedtls->take_off=tube_takeoff;
        position=strtok(NULL," ");
        sscanf(position,"%lf",&tube_beryllium);
        ptubedtls->d_beryllium=tube_beryllium;
        ptubedtls->next=NULL;
 
        /* debugging 
        printf("%s %s %s %f %f %f %f\n", el_name, line_symb, tube_el,
                tube_volts, tube_millis, tube_takeoff, tube_beryllium); /* */


        /* create filter data block and fill it */
        pfiltdtls= (struct filter_details *)malloc(sizeof(struct
                                        filter_details));
        if (pfiltdtls==NULL)
        {
           printf("memory exhausted in allocation filter details\n");
           exit(1);
        }     
        ptubedtls->filter=pfiltdtls;
        position=strtok(NULL," ");
        strcpy(filt_el,position); 
        position=strtok(NULL," ");
        sscanf(position,"%lf",&filt_mth);
        strcpy(pfiltdtls->filt_elem,filt_el);
        pfiltdtls->filt_mth=filt_mth;
        pfiltdtls->next=NULL;

        /* create xtal data block and fill it */
        pxtaldtls= (struct xtal_details *)malloc(sizeof(
               struct xtal_details));
        if (pxtaldtls==NULL)
        {
          printf("Memory exhausted in allocation xtal details\n");
          exit(1);
        }
        pfiltdtls->xtal=pxtaldtls;

        position=strtok(NULL," ");
        strcpy(xtal_name,position);
        strcpy(pxtaldtls->xtal_name,position);
        pxtaldtls->next = NULL;

        /* create collimator data block and fill it */
        pcollimdtls=(struct collim_details *) malloc(
                          sizeof(struct collim_details));

        if (pcollimdtls==NULL)
        {
          printf("Memory exhausted in alloc for collim_details\n");
          exit(1);
        }
        pxtaldtls->collimator=pcollimdtls;

        position=strtok(NULL," ");
        collim=position[0];
        pcollimdtls->type=position[0];
        pcollimdtls->next=NULL;

        /* create counter data block and fill it */
        pcounterdtls=(struct counter_details *)malloc(sizeof(
                      struct counter_details));
        if (pcounterdtls==NULL)
        {
          printf("Memory exhausted in alloc for counter details\n");
          exit(1);
        }
        position=strtok(NULL," ");
        counter=position[0];
        pcollimdtls->counter=pcounterdtls;
        position=strtok(NULL," ");
        sscanf(position,"%lf",&pulshi);
        position=strtok(NULL," ");
        sscanf(position,"%lf",&pulslo);
        position=strtok(NULL," \n");
        sscanf(position,"%lf",&pi);
        pcounterdtls->type=counter;
        pcounterdtls->pulslo=pulslo;
        pcounterdtls->pulshi=pulshi;
        pcounterdtls->pi=pi;
        pcounterdtls->next=NULL;
        
        /* debugging 
        printf("%s %s %s %6.2f %6.2f %6.2f %6.2f %s %6.2f"
               " %s %c %c %6.2f %6.2f %6.2e\n",
               el_name, line_symb, tube_el, tube_volts, tube_millis,
               tube_takeoff, tube_beryllium,
               filt_el, filt_mth, xtal_name, collim, counter,
               pulshi, pulslo, pi); /* */
     }
      
   }
}
}


