/*************************
*
*      COPYRIGHT (C) M.BOS 1998, FOR DETAILS SEE COPYING
*
*
*      xrfluor.h
*      prototypes of functions
*      for XRFLUOR
*
*      M.Bos
*      Oct 1997
*
****************************/

#define LINELENGTH 1000
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
#define PNTS_SPECTRUM 800
#define MAX_LINES 100
#define MAX_CHAR_COMP 60
#define PI (4*atan(1.0))    /*  */
#define NMAX1  24	
#define NMAX2 24
#define NR_ENERGIES 58
#define NR_THETAS 98
#define E_ELECTRON 511.0  /* mC2 electron at rest */
#define R0 2.8179380e-13 
#define Na 6.02486e23   /* Avogadro */
#define  Compton       0.0242631  /* compton wavelength of electron in A */
#define MAX_SMPLS 40

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
}  elem_list ;

typedef struct
{
         int z;  /* atomic number */
         char line[6];
} elem_line ;

typedef struct
{
         int z;  /* atomic number */
         char line[6];
         double kev;
         int filt_z ;
         double mth_filt;
} meas_line ;

struct spectrum_pair
  {
     double lambda;
     double intens;
   };


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
int get_z(char *, elem_list *);
double get_rho(int, elem_list *);
double get_atw(int, elem_list *);
double attn(double mualamba, double psi1, double muai, double psi2, double da);
double primintmono(double G, double psi1, double psi2,
                          double Ci, double Jilambda ,double omegai,
                          double gi, double tauilambda,
                          double mu1lambda, double mu1i);
double primintfull(double, double, double,  double, double, double, double);
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
double pilabda(char *, char *, double, layer *);
double sijlambda(char *, char *, char *, char *,
                 layer *, double, double );
double linfty(double,double,double);
int find_enhanc(char *, char *, layer *, double);
double intrasecflu(char *, char *, double, multilayer *, int, double);
struct spectrum_pair *gen_spec(char *, double, double, double, int, int, double);
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
void checkcomps(multicomplayer *);
void bldedges(FILE *, double edges[][10]);
void bldnkls(FILE *, double nkl123[][3]);
void bldfij(FILE *, FILE *,double fij[][13]);
void def_sample(FILE *, multilayer *, multicomplayer *);
void convcomp(multicomplayer *, multilayer *);
double backscatter(int, double);
double comptfot_func1(double, double);
double comptfot_func2(double, double);
double raylfot_func1(double, double);
double raylfot_func2(double, double);
double fotrayl_func1(double, double);
double fotrayl_func2(double, double);
double comptcrossmix(double, layer *, elem_list *, double);
double raycrossmix(double, layer *, elem_list *, double);
double comptisesr(char *, char *, double, layer *, int, int);
double raylisesr(char *, char *, double, layer *, int, int);
double ispfr(char *, char *, double, layer *, int, int);
double simpsdbl(double (*f)(double, double), double, double, double, double,
			    int, int);
int storeffs(void);
double calcsxz(int, double, double);
double calcff2(int, double, double);
