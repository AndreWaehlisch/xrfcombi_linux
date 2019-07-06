/*************************
*
*
* Copyright (C) 1998, M.Bos, for details see copying
*
*      xrfluor.h
*      prototypes of functions
*      for XRFLUOR
*
*      M.Bos
*      Oct 1997
*
****************************/
#ifndef _XRFLUOR
#define _XRFLUOR

#define LINELENGTH 300
#define MAXCOEFS 200
#define HIGHEST_Z 94
#define ALL_ELMS 94
#define TOTAL 103
#define TOT_RELRATES 28
#define NR_KLINES 6
#define NR_L1LINES 5
#define NR_L3LINES 5
#define NR_L2LINES 4
#define NR_M5LINES 3
#define NR_M4LINES 2
#define NR_M3LINES 1
#define NR_M2LINES 1
#define NR_OMEGAS 8
#define MAX_ENH_LINES 1000
#define MAX_OTHER_ELMS 40
#define LAMBDA_HIGHEST 25.0 /* highest wavelength to be used in excit */
#define PNTS_SPECTRUM  800
#define MAX_LINES 30
#define MAX_CHAR_COMP 60
#define NPARM  25    /* max number of free parameters to be determined */
#define PI 4*atan(1.0) /* */


typedef struct 
{
	char name[3];
	double wfract ;
} elem 
;

typedef struct
{
   char name[3];
   int count ;
} elem_cnt
;

typedef struct
{
   char name[MAX_CHAR_COMP];
   double wfract;
   int wfixed;
   int nr_elems;
   elem_cnt *elementen;
} compound
;



typedef struct 
{
	int nr_elements ;
	double massthickness ;
        elem *elementen   ; 
} layer
;

typedef struct
{
	int nr_layers ;
	layer **lagen;
} multilayer
;

typedef struct
{
   int nr_compounds;
   double massthickness ;
   int mthfixed;
   compound *comps;
} complayer
;

typedef struct
{
  int nr_layers ;
  complayer **lagen;
} multicomplayer
;

typedef struct 
{
	double energy;
	double abscoef ;
	double slope,A;
} emupair 
;

typedef struct 
{
        int z;
	int nr_emupairs ;
        emupair  *paren   ; 
} atom
;

typedef struct
{
	int at_nrs ;
	atom **atomcoefs;
} emudata
;


typedef struct
{
	int z;
	char sym[3];
	double at_weight;
	double density;
} atom_dat ;

typedef struct
{
        int nr_elements ;
        atom_dat *patom_dat ;
		char inverse_lookup[1024];
}  elem_list ;

typedef struct
{
         int z;  /* atomic number */
         char line[6];
} elem_line ;

typedef struct
{
     int z;
     char line[6];
     double kev;
     int filt_z;
     double mth_filt;
} meas_line ;

typedef struct
{
   char *key;
   int z_elem ;
 } ENTRY;



struct spectrum_pair
  {
     double lambda;
     double intens;
   };

struct pstruct {
   double val;
   double parm[NPARM]; 
} ;

struct qstruct {
   int parmndx[NPARM];
   double q[NPARM];
   double yplus[NPARM];
   double yminus[NPARM];
   double a[NPARM];
   double bdiag[NPARM];
   double inv_bdiag[NPARM];
   double std_dev[NPARM];
} ;


struct counter_details
{
   char type ;
   double pulslo ;
   double pulshi ;
   double pi ;   
   struct counter_details *next;
} ;

struct collim_details
{
   char type ;
   struct counter_details *counter; 
   struct collim_details *next ;
} ;

struct xtal_details
{
   char xtal_name[7];
   struct collim_details *collimator ;
   struct xtal_details *next;
} ;

struct filter_details
{
   char filt_elem[5] ;
   double filt_mth ;
   struct xtal_details *xtal;
   struct filter_details *next;
} ;

struct tube_details
{
   char tube_el[3];
   double voltage ;
   double current ;
   double take_off;
   double d_beryllium ;
   struct filter_details *filter;
   struct tube_details *next;
} ;

struct line_details
{
   char line_name[6] ;
   struct tube_details *tube;
   struct line_details *next;
} ;


void buildlay(FILE *, multilayer *);
void checkdata(multilayer *);
int get_nr_layers(multilayer *);
int get_nr_elements(layer *);
layer *get_layer(multilayer *, int);
char *get_elem(layer *, int);
int get_nr_elem(layer *, char *);
double get_weight_fraction(layer *,int);
double get_mthickness(layer *);
void check(multilayer *);
void bldemutable(FILE *, emudata *);
void chkdata(emudata *);
double calc_muilambda(int, double, emudata *);
double calc_mumix(double, layer *, elem_list *, emudata *);
double get_weightfract(layer *, int);
void bld_atlist(FILE *, elem_list *);
void chk_atlist(elem_list *);
char *get_asymb(int, elem_list * ); 
extern char reverse_lookup[32767];
#define get_z(a,b) reverse_lookup[(((int)a[0])<<8)+(int)a[1]]
/*int get_z(char *, elem_list *);
*/
double get_rho(int, elem_list *);
double get_atw(int, elem_list *);
double attn(double mualamba, double psi1, double muai, double psi2, double da);
double primintmono(double G, double psi1, double psi2,
                          double Ci, double Jilambda ,double omegai,
                          double gi, double tauilambda,
                          double mu1lambda, double mu1i);
double selfattn(double mu1lamba,  double mu1i,  double d);
void bldrelrates(double relrates[][TOT_RELRATES]);
double getrelrate(int , char *, double relrates[][TOT_RELRATES]);
void bldjmpedg(FILE *, double jmpdata[][2]);
double getjmpfact( int z, int hole_nr, double energy, double jmpdata[][2]);
void bldomega(FILE *, FILE *, FILE *, double omegas[][NR_OMEGAS]);
double get_omega(int, char *, double omegas[][NR_OMEGAS]);
char *line_to_shell(char *line, char *linenames[27], char *shellnames[27]);
int get_hole_nr(char *);
void bld_k_lines( FILE *fpapdx1, double klines[][6]);
void bld_l_lines( FILE *fpapdx2, double llines[][14]);
void bld_m_lines( FILE *fpapdx3, double mlines[][7]);
double characwl(char *, char *, elem_list *);
void init_tables(void);
void set_spectrom(void);
double E1z(double, int );
double Eiz(double, int );
double dfunc(double);
double lzero(double, double, double, double);
double ltot(double, double, double, double);
double pilabda(int , int , layer *);
double sijlambda(char *, char *, char *, char *,
                 layer *, double, double, int, int );
double linfty(double,double,double);
int find_enhanc(char *, char *, layer *, double);
double intrasecflu(char *, char *, double, multilayer *, int, double, int, int);
struct spectrum_pair *gen_spec(char *, double, double, double,int,int, double);
int get_mslines(FILE *);
layer *mkbulklay(char *);
void show_enhances(int, elem_line **);
void clean_up_enh(int);
void bldcompos(FILE *, multicomplayer *);
int get_nr_clay(multicomplayer *);
int get_nr_comps(complayer *);
complayer *get_clayer(multicomplayer *, int);
char *get_comp_name(complayer *, int);
double get_cweightf(complayer *, int);
int get_fixwght(complayer *, int);
double get_cmth(complayer *);
int get_fixmth(complayer *);
double convcomp(multicomplayer *,multilayer *);
void checkcomps(multicomplayer *);
void bldedges(FILE *, double edges[][10]);
void bldnkls(FILE *, double nkl123[][3]);
void bldfij(FILE *, FILE *,double fij[][13]);
void get_imeas(FILE *, int, double *);
FILE *own_fopen(char *);
int put_params(multicomplayer *, double *);
int get_params(multicomplayer *, double *);
void show_params(double *);
int func(struct pstruct *);
int pvalcomp(int , int);
int p_swap(int, int);
int centroid(struct pstruct *, struct pstruct *psrc[], int);
int ptrial(struct pstruct *, struct pstruct *, struct pstruct *,
              double , double );
void fpprint(FILE *, struct pstruct *);
void fsprint(FILE *, char *);
void pcopy(struct pstruct *, struct pstruct *);
void fmatprint(FILE *, double nep[][NPARM]);
void fqprint(FILE *, double *);
void bsort(int);
void show_pars(struct pstruct *);
int simpfit(FILE *);
void ffitprint(FILE *);
void fquadprint(FILE *);
int simpdev(FILE *);
void free_multil(multilayer *);
void fillmulabda(multilayer *, int);
void check_mu_el(int, int);
double calc_mumix_snel(int, int, layer *, elem_list *, emudata *);
multilayer *mkbinbulk(char *, char *, double);
multilayer *mktertbulk( char *, char *, char *, double, double);
void conciter(multilayer *);
double inflfactor(int, int, elem *, int, int *);
double newton(multilayer *, int);
void copysample(multilayer *, multilayer *);
double func_el(multilayer *);
double sum_el(multilayer *, int);
double wilabdak(multilayer *, int, int, int);
double k_j(int, int, int, char *, int);
double delta_ij(int, int, int, int, multilayer *);
double calc_alfa_ij(multilayer *, int, int, int);
double calc_rho_ij(int, int, int, multilayer *);
double p_ij(int, int, int, char *, int, multilayer *);
void calc_alphas(multilayer *, int);
void criss(multilayer *);
int prtmcomp(FILE *, multicomplayer *);
void alphrho(multilayer *, int);
void phase2(multilayer *);
void get_pifacts(FILE *, int, double *);
int get_gifactors(FILE *);
int get_nr_lines_meas(FILE *);
void calc_rikmeas(double *, int, int);
#endif

