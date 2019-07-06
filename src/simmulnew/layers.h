/*
*	layers.h
*	datastructure for access routines
*	for multi-layer thin film samples
*
*	M.Bos
*	august 1997
*/

typedef struct 
{
	char name[3];
	double wfract ;
} elem 
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

