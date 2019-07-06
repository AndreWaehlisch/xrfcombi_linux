/*************************************************************
*
*  COPYRIGHT 1999 (C) M.Bos, for details see the file copying
*
*  File bldsetpis.c
*
*  Program to build datastructure
*  for characteristic fluorescence lines
*  containing all measurement condition data
*  plus factors Pi cf Wang, X-Ray spectrom. 25 (1996) 245
*
*  M.Bos
*  febr 1999
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
struct counter_details  *new_counter(char, double, double, double);
struct collim_details *new_collim(char, char, double, double, double);
struct xtal_details *new_xtal(char *, char, char,
                              double, double, double);
struct filter_details *new_filter(char *, double,
                                  char *, char, char,
                             double, double, double);
struct tube_details *new_tube(char *, double, double, double,
                        double, char *, double,
                        char *, char, char,
                        double, double, double);
struct line_details *new_line(char *, char *, double, double, double,
                        double, char *, double,
                        char *, char, char,
                        double, double, double);

void bld_setpis(FILE *fp)
{

    int i;
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
       
       /* get line_name  */
       position=strtok(NULL," ");
       /* debugging 
       printf("line symbol %s \n",position); /* */
       strcpy(line_symb, position);
        
       /* read the tube details now */        
       position=strtok(NULL," ");
       strcpy(tube_el,position);
       position=strtok(NULL," ");
       sscanf(position,"%lf",&tube_volts);
       position=strtok(NULL," ");
       sscanf(position,"%lf",&tube_millis);
       position=strtok(NULL," ");
       sscanf(position,"%lf",&tube_takeoff);
       position=strtok(NULL," ");
       sscanf(position,"%lf",&tube_beryllium);
 
        /* debugging 
        printf("%s %s %s %f %f %f %f\n", el_name, line_symb, tube_el,
                tube_volts, tube_millis, tube_takeoff, tube_beryllium); /* */

       /* get filter data */ 
        position=strtok(NULL," ");
        strcpy(filt_el,position); 
        position=strtok(NULL," ");
        sscanf(position,"%lf",&filt_mth);

        /* get xtal data  */

        position=strtok(NULL," ");
        strcpy(xtal_name,position);

        /* get collimator data  */

        position=strtok(NULL," ");
        collim=position[0];

        /* get counter data  */
        position=strtok(NULL," ");
        counter=position[0];
        position=strtok(NULL," ");
        sscanf(position,"%lf",&pulshi);
        position=strtok(NULL," ");
        sscanf(position,"%lf",&pulslo);
        position=strtok(NULL," \n");
        sscanf(position,"%lf",&pi);
        
        /* debugging 
        printf("%s %s %s %6.2f %6.2f %6.2f %6.2f %s %6.2f"
               " %s %c %c %6.2f %6.2f %6.2e\n",
               el_name, line_symb, tube_el, tube_volts, tube_millis,
               tube_takeoff, tube_beryllium,
               filt_el, filt_mth, xtal_name, collim, counter,
               pulshi, pulslo, pi); /* */

     /* now get pointer to node of element with z as atomic number */
     plindtls=el_details[z];
     if (plindtls==NULL)
     {
       /* element not yet is dataset */
       plindtls=new_line(line_symb, tube_el, tube_volts, tube_millis,
                         tube_takeoff, tube_beryllium, filt_el, filt_mth,
                         xtal_name, collim, counter, pulshi, pulslo, pi);
       el_details[z]=plindtls;
     }
     else
     {
       /* now look for line_name entry in linked list line_details */
       while (strcmp(plindtls->line_name, line_symb))
       {
         if (plindtls->next != NULL)
           plindtls=plindtls->next;
         else
          {
            plindtls->next=new_line(line_symb, tube_el, tube_volts, tube_millis,
                         tube_takeoff, tube_beryllium, filt_el, filt_mth,
                         xtal_name, collim, counter, pulshi, pulslo, pi);
            plindtls=plindtls->next;
          }
       }
       ptubedtls=plindtls->tube;
       while (strcmp(tube_el, ptubedtls->tube_el) ||
               (tube_volts != ptubedtls->voltage) ||
               (tube_millis != ptubedtls->current) ||
               (tube_takeoff != ptubedtls->take_off) ||
               (tube_beryllium != ptubedtls->d_beryllium))
       {
         if (ptubedtls->next != NULL)
           ptubedtls=ptubedtls->next;
         else
          {
            ptubedtls->next=new_tube(tube_el, tube_volts, tube_millis, 
                                    tube_takeoff, tube_beryllium, filt_el,
                                    filt_mth, xtal_name, collim, counter,
                                    pulshi, pulslo, pi);
            ptubedtls=ptubedtls->next;
          }
       }
       pfiltdtls=ptubedtls->filter;

       while (strcmp(filt_el, pfiltdtls->filt_elem) ||
              (filt_mth != pfiltdtls->filt_mth))
       {
         if (pfiltdtls->next != NULL)
            pfiltdtls=pfiltdtls->next;
         else
          {
            pfiltdtls->next=new_filter(filt_el, filt_mth, xtal_name, collim,
                                      counter, pulshi, pulslo, pi);
            pfiltdtls=pfiltdtls->next;
          }
       }
       pxtaldtls=pfiltdtls->xtal;
  
       while (strcmp(xtal_name, pxtaldtls->xtal_name))
       {
         if (pxtaldtls->next != NULL)
           pxtaldtls=pxtaldtls->next;
         else
          {
            pxtaldtls->next=new_xtal(xtal_name, collim, counter,
                                     pulshi, pulslo, pi);
            pxtaldtls=pxtaldtls->next;
          }
       }
       pcollimdtls=pxtaldtls->collimator;

       while (collim != pcollimdtls->type)
       {
         if (pcollimdtls->next != NULL)
           pcollimdtls=pcollimdtls->next;
         else
         {
          pcollimdtls->next=new_collim(collim, counter, pulshi, pulslo, pi);
          pcollimdtls=pcollimdtls->next;
         }
       }
       pcounterdtls=pcollimdtls->counter;

       while (( counter != pcounterdtls->type) ||
              (pulshi != pcounterdtls->pulshi) ||
              (pulslo != pcounterdtls->pulslo))
       {
         if (pcounterdtls->next != NULL)
            pcounterdtls=pcounterdtls->next;
         else
          {
            pcounterdtls->next=new_counter(counter, pulshi, pulslo, pi);
            pcounterdtls=pcounterdtls->next;
          }
       } 
    }
   }
}



struct counter_details  *new_counter(char counter, double pulshi,
                double pulslo, double pifactor)
{
    struct counter_details *pcounterdtls ;

    /* create counter data block and fill it */

    pcounterdtls=
      (struct counter_details *)malloc(sizeof( struct counter_details));
    if (pcounterdtls==NULL)
    {
      printf("Memory exhausted in alloc for counter details\n");
      exit(1);
    }
   
    /* now fill data structure */
    pcounterdtls->type=counter;
    pcounterdtls->pulslo=pulslo;
    pcounterdtls->pulshi=pulshi;
    pcounterdtls->pi=pifactor;
    pcounterdtls->next=NULL;
    return(pcounterdtls);
}

struct collim_details *new_collim(char collim, char counter,
        double pulshi, double pulslo, double pifactor)
{
 
   struct collim_details *pcollimdtls;
   struct counter_details *pcounterdtls;
   /* create collimator data block and fill it */

   pcollimdtls=(struct collim_details *)malloc(
                 sizeof(struct collim_details));
   if (pcollimdtls==NULL)
   {
     printf("Memory exhausted in alloc for collim_details\n");
     exit(1);
   }
   pcollimdtls->type=collim;
   pcollimdtls->next=NULL;
   pcounterdtls=new_counter(counter, pulshi, pulslo, pifactor);
   pcollimdtls->counter=pcounterdtls;
   return(pcollimdtls);
}


   
struct xtal_details *new_xtal(char *xtal_name, char collim, char counter,
        double pulshi, double pulslo, double pifactor)
{
  
   struct xtal_details *pxtaldtls;
   struct collim_details *pcollimdtls;

   /* create xtal datablock and fill it */ 
   pxtaldtls= (struct xtal_details *)malloc(sizeof( struct xtal_details));
   if ( pxtaldtls==NULL)
   {
      printf("Memory exhausted in alloc xtal details \n");
      exit(1);
   }

   strcpy(pxtaldtls->xtal_name,xtal_name);
   pxtaldtls->next=NULL;
   pcollimdtls=new_collim(collim, counter,pulshi, pulslo, pifactor);
   pxtaldtls->collimator=pcollimdtls;
   return(pxtaldtls);
}

struct filter_details *new_filter(char *filt_el, double filt_mth,
                                  char *xtal_name, char collim, char counter,
                             double pulshi, double pulslo, double pifactor)
{
   struct filter_details *pfiltdtls;
   struct xtal_details *pxtaldtls;

   /* create datablock for filter and fill it */

   pfiltdtls= (struct filter_details *)malloc(sizeof(struct filter_details));
   if (pfiltdtls==NULL)
   {
     printf("Memory exhausted in alloc filter details\n");
     exit(1);
   }
   strcpy(pfiltdtls->filt_elem, filt_el);
   pfiltdtls->filt_mth=filt_mth;
   pfiltdtls->next=NULL;
   pxtaldtls=new_xtal(xtal_name, collim,  counter,
                      pulshi,  pulslo, pifactor);
   pfiltdtls->xtal=pxtaldtls;
   return(pfiltdtls);
}


struct tube_details *new_tube(char *tube_el, double tube_volts,
                        double tube_millis, double tube_takeoff,
                        double tube_beryllium,
                        char *filt_el, double filt_mth,
                        char *xtal_name, char collim, char counter,
                        double pulshi, double pulslo, double pifactor)
{
   struct tube_details *ptubedtls;
   struct filter_details *pfiltdtls;

   /* create data block for tube details and fill it */
   ptubedtls=(struct tube_details *)malloc(sizeof(struct tube_details));
   if (ptubedtls==NULL)
   {
      printf("memory exhausted in alloc for tube details \n");
      exit(1);
   }

   strcpy(ptubedtls->tube_el,tube_el);
   ptubedtls->voltage=tube_volts;
   ptubedtls->current=tube_millis;
   ptubedtls->take_off=tube_takeoff;
   ptubedtls->d_beryllium=tube_beryllium;
   ptubedtls->next=NULL;
   pfiltdtls=new_filter(filt_el, filt_mth, xtal_name, collim,  counter,
                        pulshi,  pulslo,  pifactor);
   ptubedtls->filter=pfiltdtls;
   return(ptubedtls);
}


struct line_details *new_line(char *line_symb, char *tube_el,
                        double tube_volts,
                        double tube_millis, double tube_takeoff,
                        double tube_beryllium,
                        char *filt_el, double filt_mth,
                        char *xtal_name, char collim, char counter,
                        double pulshi, double pulslo, double pifactor)
{
    struct line_details *plindtls;
    struct tube_details *ptubedtls;

    /* create data block for line details and fill it */

    plindtls=(struct line_details *)malloc(sizeof(struct line_details));
    if (plindtls==NULL)
    {
      printf("Memory exhausted in alloc for line details\n");
      exit(1);
    }
    strcpy(plindtls->line_name, line_symb);
    plindtls->next=NULL;
    ptubedtls=new_tube(tube_el, tube_volts, tube_millis,  tube_takeoff,
                       tube_beryllium, filt_el,  filt_mth,
                       xtal_name, collim, counter,
                       pulshi, pulslo,  pifactor);
    plindtls->tube=ptubedtls;
    return(plindtls);
}

    


