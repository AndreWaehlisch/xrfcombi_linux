/*****************************************************
*
*    xrfdefs.h
*    global vars with data
*    for XRFLUOR
*
*    M.Bos
*    Sept 1997
*
******************************************************/
    

elem_list lijst_els ;

char *linenames[27] ={ "ka", "ka1", "ka2", "kb1", "kb3", "kb2",
                            "lb3","lb4","lg2","lg3","lg4",
                            "lb1","lg1","lg6","leta",
                            "la1", "la2", "lb2","lb5","ll",
                            "mx3","mg",
                            "mb","mx2",
                            "ma1", "ma2", "mx1"};

char *lines[27] ={ "ka", "ka1", "ka2", "kb1", "kb3", "kb2",
                  "la1", "la2", "lb1", "lb2", "lb3", "lb4",
                  "lb5", "lg1", "lg2", "lg3", "lg4", "lg6",
                  "ll", "leta",
                  "ma1", "ma2", "mb","mg", "mx1", "mx2", "mx3"} ;  

char *shellreeks[27]={ "k","k", "k", "k", "k","k",
                         "l1","l1","l1","l1", "l1",
                         "l2","l2", "l2", "l2",
                         "l3","l3", "l3", "l3", "l3",
                          "m2","m3","m4","m4", "m5", "m5","m5"};

/* vars to hold spectrometer settings */
double Gfactor, psi1, psi2  ;
 
double klines[98][6], llines[98][14], mlines[98][7];
double relrates[TOTAL][TOT_RELRATES] ;
double jmpdata[ALL_ELMS][2];
double omegas[TOTAL][NR_OMEGAS] ;
double jmpdata[ALL_ELMS][2];
emudata abstabel;
elem_line *enhances [MAX_ENH_LINES] ;
char anode_el[3];
double berryliumd, take_off;
double pellas_consts[4][3] = { {3.22E6 , 9.76E4, -0.39},
                              {5.13E5 , 2.05E5, -0.014},
                              {2.02E7 , 2.65E6 , 0.21},
                              {1.76E7, 6.05E6, -0.09}};

meas_line lines_to_meas[MAX_LINES];
int nr_of_meas_lines ;
double edges[TOTAL][10] ;
double nkl123[TOTAL][3];
double fij[TOTAL][13];
double sumfracts[MAX_LAYERS];
