#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include <string.h>
#include <string>
#include <omp.h>
#include <cstdlib> 
#include <ctime>
#include <vector>
#include <time.h>
#include <complex>
using namespace std;

#define caseNum 	3

// Variables
#define XLENGTH		110
#define YLENGTH		87
#define ZLENGTH		103
#define i_main3     XLENGTH*YLENGTH*ZLENGTH // node °³¼ö


// Fluid properties
#define tau0		1.0 // relaxation time

// BC, IC
#define uMax		0.01 // BC
#define RHO_0		1.0 // IC

// LBM settings
#define Q			19 // D3Q19
#define w0			1.0 / 3.0
#define w1			1.0 / 18.0
#define w2			1.0 / 36.0
#define T			30001 // simulation time (duration)
#define T_save		500 // save period

// Carreau-Yasuda model (blood viscosity)
#define CY_n		0.3568
#define CY_mu0		0.056 // 2.1867
#define CY_muinf	0.0035 // 0.1367
#define CY_lambda	3.313 // 8484.512
#define CY_a		2

// Domain
#define i_main2		(XLENGTH*YLENGTH*ZLENGTH*Q+1) // distribution
#define i_main3		(XLENGTH*YLENGTH*ZLENGTH+1) // node
#define i_main4		(XLENGTH*YLENGTH*ZLENGTH*3*3+1) // StrainRateTensor

// Macro
#define SQ(x)		((x)*(x)) // square


// Pointer
double	*f		= new double[i_main2];
int 	*BCnode = new int[i_main3]; // 0: wall, 1: fluid, 2: inlet, 3: outlet


// Modules
void 	plt_load(void);
void 	Log_timestep(void);

int main(void)
{
	plt_load();
	Log_timestep();
}

void 	plt_load(void)
{
    FILE *f0;
    char buffer[100];
    int x, y, z;
    int x1, y1, z1;
    int i;
    int i_nd;

    if((f0 = fopen("coronary.plt", "r")) == NULL) // vox_dim 1.0, cropped ÈÄ
    {
        system("PAUSE");
        return;
    }

    fgets(buffer, sizeof(buffer), f0);
	fgets(buffer, sizeof(buffer), f0);

    // initialize
    for (z=0; z<ZLENGTH; z++) for (y=0; y<YLENGTH; y++) for (x=0; x<XLENGTH; x++)
    {
        i_nd = z + ZLENGTH * (y + YLENGTH * x);
        *(BCnode + i_nd) = 0;
        fscanf(f0, "%d\t%d\t%d\t%d\n", &x1, &y1, &z1, &*(BCnode + i_nd));
    }

    fflush(f0);
	fclose(f0);


}

void 	Log_timestep(void)
{
    char FileName_dummy[100];
    FILE *fP3;

    sprintf(FileName_dummy, "results/case%d/case_%d_timestep.txt", caseNum, caseNum);
    fP3 = fopen(FileName_dummy, "a");
    fprintf(fP3,"%d\n", T_save);
    fflush(fP3);
    fclose(fP3);

	return;
}