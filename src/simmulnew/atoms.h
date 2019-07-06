/*
*	atoms.h
*	datastructures for access routines
*	for data of the elements
*	especially abs. coeffs
*
*	M.Bos
*	august 1997
*/
#define HIGHEST_Z 18

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


