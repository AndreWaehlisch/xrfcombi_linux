/*****************************************************
*
*
* Copyright (C) 1998, M.Bos, for details see copying
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
double Gfactor, psi1, psi2, kev ;
 
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
int nr_free_vars ;
double *rxi_meas, *rxi_calc, *rxi_calc2, *intbulk_calc, *intprep_calc ;
multicomplayer compsample ;
multilayer sample ;
multilayer sample2 ;
struct spectrum_pair **pntrspectra ;
struct pstruct *p, pcent, **p_p, pmin ;
double exit_test, mean_func, rms_func, test, rms_data;
int iter, maxiter, nparm, nvert, ndata, nfree ;
int maxquad_skip, prt_cycle , nquad_skip, quad_cycle;
double quad_test ;
struct qstruct q ;
double yzero, ymin, ypmin, mse ;
double qmat [NPARM][NPARM];
char  title[132];
double muzlambda[MAX_LINES][TOTAL][PNTS_SPECTRUM];
double *kvlines ;
double **aij, **aijj, ***aijk ;
double **alpha_ij, **rho_ij;

	/* extra globals for double integral calcn */
double sinpsi2, cospsi2,  cospsi1,  EEi, muci, muclabda, sinpsi1 ;
layer *ptrlaag ;
int f2cnt ;
double store_labda, store_labdai ;
int store_z ;
double store_tauilambda ;
double store_energy ;
double store_thickness;
double fluor_grens;

double cross [100][60][100];
double energies [100][60];
double thetas [100][100];

