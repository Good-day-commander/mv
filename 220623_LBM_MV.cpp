#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cmath>
	/*Care for streaming in tables and update x in the extrapolation */
#include <string.h>
#include <string>
#include <omp.h>
#include <cstdlib> 
#include <ctime>
#include <vector>
#include <time.h>
#include <complex>


// 0.4um�� 1ĭ
// 8um¥��


#define		XLENGTH2     	110
#define		YLENGTH2     	87
#define		ZLENGTH2     	103


#define		XLENGTH     	ZLENGTH2
#define		YLENGTH     	YLENGTH2
#define		ZLENGTH     	XLENGTH2


#define		tau0           1.0 // relax time

#define		uMax           0.01//0.001//(0.01333) OO 0.03
#define     forcepowerAC    (0) // IBM ?

// RBC (IBM) ///////////////////////////
/*


#define		rbcn   1
#define		p_num       500
#define     f_main      996
#define     l_main       1494

#define		tp_num       (rbcn*p_num)
#define     tf_main       (rbcn*f_main)
#define     tl_main       (rbcn*l_main)
*/
//////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define		Q           	19
#define		w0              1.0 / 3.0
#define		w1              1.0 / 18.0
#define		w2              1.0 / 36.0
#define		RHO_0	        1.0 // IC


#define		T				   30001
#define		Twrite			   1
#define		ppptime             500//200//   data
#define		afterprint          1//
#define		reintime        0


#define		seng	        0.3568
#define		mu0             2.1867//0.056
#define		muinf           0.1367//0.0035
#define		ramda           8484.512//3.313

#define	    i_main2         (XLENGTH*YLENGTH*ZLENGTH*Q+1)
#define     i_main3         (XLENGTH*YLENGTH*ZLENGTH+1)
#define     i_main4         (XLENGTH*YLENGTH*ZLENGTH*3*3+1)
#define     i_main3_2       (XLENGTH2*YLENGTH2*ZLENGTH2+1)



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void        tables                              ( void );
void        propagate                           ( void );
void        initialise                          ( void );
void        calc_obs                            ( void );
void        project                             ( void );
void        density_data                        ( void );
void		timestep_data						( void );
void	    pulsatile                           ( void );
void		constructimage						( void );


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double		*b = new double[i_main2];
double		*back = new double[i_main2];
double		*beq = new double[i_main2];


int	       next_x         [XLENGTH+1][Q];
int	       next_y	      [YLENGTH+1][Q];
int	       next_z	      [ZLENGTH+1][Q];
int         cx                          [Q] = { 0,-1,-1,-1,-1,-1, 0, 0, 0, 0,1,1,1, 1, 1,0,0, 0,0};
int         cy                          [Q] = { 0,-1, 0, 0, 0, 1,-1,-1,-1, 0,1,0,0, 0,-1,1,1, 1,0};
int         cz                          [Q] = { 0, 0,-1, 0, 1, 0,-1, 0, 1,-1,0,1,0,-1, 0,1,0,-1,1};

int	         v	        [Q][3];
int	        *bnode = new int[i_main3];
int	        *b_y = new int[i_main3];
int	        *b_z = new int[i_main3];
int	        *b_b = new int[i_main3];
int	        *bnode2 = new int[i_main3_2];
int	        *bnode3 = new int[i_main3_2];

double		t                           [Q] = {w0,w2,w2,w1,w2,w2,w2,w1,w2,w1,w2,w2,w1,w2,w2,w2,w1,w2,w1};

double      *rhoN = new double[i_main3];    
double	    *ux	= new double[i_main3];         
double	    *uy = new double[i_main3];          
double	    *uz = new double[i_main3];  
double	    *uu = new double[i_main3]; 
double	    *rho = new double[i_main3];        


double		*RGradient	  = new double[XLENGTH];

double      *tau = new double[i_main3];  
double      *Shear_stress = new double[i_main3];  
double      *mu = new double[i_main3];  
double      *summ = new double[i_main3];  
double      *strate = new double[i_main4];
double      *real_shear_x = new double[i_main3]; 
double      *real_shear_y = new double[i_main3];
double      *real_shear_z = new double[i_main3];
double      *real_shear_A = new double[i_main3];


double      *OSI_x1 = new double[i_main3]; 
double      *OSI_y1 = new double[i_main3];
double      *OSI_z1 = new double[i_main3];
double      *OSI_x2 = new double[i_main3]; 
double      *OSI_y2 = new double[i_main3];
double      *OSI_z2 = new double[i_main3];
double      *OSI = new double[i_main3];

double	    *press = new double[i_main3];
double		*press_diff = new double[i_main3];

double      pulse           [20];



//Real IBM
double	    *FforceX = new double[i_main3];          
double	    *FforceY = new double[i_main3];          
double	    *FforceZ = new double[i_main3];  
	



////////////////////////////////////////////////////////////////////////////////////
int         Timek;
int			pulse_step;
int			c_num; // case number

double	outlet_velocity;
double	outlet_node;
double	Resistance;

double pressure1;


using namespace std;



   
int	main	(void)
{
	c_num = 3; // case number
	outlet_node = 0;
	outlet_velocity = 0;


    //int count;
	int jj;
	pulse_step=0;	
	Timek = 0;
    tables		();


    initialise	();
	constructimage();
	pulsatile();


	for (Timek = reintime; Timek < T; Timek++)
	{

		if(Timek%1000==0)
		{
			pulse_step ++;
			if(pulse_step>=20)
			{
				pulse_step = 0;
			}
		}



		propagate();
		calc_obs();
		project();


		if (Timek == 1)
		{
			density_data();
		}

		if ( Timek > afterprint)
		{
			if (Timek >= Twrite) if (Timek % ppptime == 0)
		{
			density_data();
			timestep_data();

		}
		}

		//if(Timek%10==0) timestep_data();


		
    }


    scanf("%d",&Timek);



	delete[]b;
	delete[]back;
	delete[]press;
	delete[]rhoN;

	delete[]rho;
	delete[]bnode;
	delete[]b_y;
	delete[]b_z;
	delete[]ux;
	delete[]uy;
	delete[]uz;

    delete[]beq;


	delete[]FforceX;
	delete[]FforceY;
	delete[]FforceZ;


	delete[]Shear_stress;

		return  0;
}

void    density_data (void)
{   int     x, y, z, i_nd, a0;
	double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    FILE    *f2;
    char FileNameVelocity[100];
    char string[100] = {'0'};

    sprintf(FileNameVelocity, "results/case%d/case_%d_data2_%d.plt", c_num, c_num, Timek);

   f2 = fopen(FileNameVelocity, "w");
   fprintf(f2, "Variables = X, Y, Z, Boundary, Ux, Uy, Uz, Pressure, Shear Stress, mu, Force, RealShear, OSI\n");
   fprintf(f2,"Zone I=   %d,J=   %d, k=   %d,F=POINT\n",XLENGTH, YLENGTH, ZLENGTH);
   for (z = 0; z < ZLENGTH; z++)
  {
   for ( y = 0; y < YLENGTH; y++ )
      {
        for ( x = 0; x < XLENGTH; x++ )
         {
			                  
                  {
         i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;          
		 a0 = *(bnode3 + i_nd);
         a1 = *(ux + i_nd);
		 a2 = *(uy + i_nd);
		 a3 = *(uz + i_nd);
		 a4 = *(press + i_nd);


		 if(*(bnode + i_nd)!=0)  
		 {
			a1  = 0.0;
			a2  = 0.0;
			a3  = 0.0;
			a4  = 1.0 / 3.0;
		 }
		 
		 //if(*(cnode + i_nd)==1)
		 if(*(bnode + i_nd)==0) 
		 {
			 a5 = *(Shear_stress + i_nd); 
			 a6 = *(mu + i_nd);
			 a7 = *(press + i_nd) * *(RGradient + x + 1);
 			 a8 = *(real_shear_A + i_nd);
			 a9 = *(OSI + i_nd);

		 }
		 else
		 {
			 a5  = 0.0;
			 a6  = 0.0;
			 a7  = 0.0;
			 a8  = 0.0;
			 a9  = 0.0;
		 }
		 
		 fprintf(f2,"%d\t%d\t%d\t%d\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\n",x,y,z,a0, a1,a2,a3,a4,a5, a6, a7, a8, a9);
		 //fprintf(f2,"%d\t%d\t%d\t%d\t%d\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\n",x,y,z, *(b_b + i_nd),a0, a1,a2,a3,a4);
		 //fprintf(f2,"%d\t%d\t%d\t%d\t%d\t%d\n",x,y,z,a0, *(b_b + i_nd), *(b_z + i_nd));
         fprintf(f2, "\n");
			 }
       }
   }
   }
   
 
        fflush(f2);
    fclose(f2);


    return;
}

void 	timestep_data()
{
	FILE *fP3;
	FILE *fP4;
	FILE *fP5;
	FILE *fP6;

	int x1, x2,y,z, i_in1, i_in2, i;
	double ffr1, ffr2, ffr3;
	int y1, y2, y3;
	int i_nd11, i_nd12, i_nd21, i_nd22, i_nd31, i_nd32;

	x1 = 10;	
	x2 = XLENGTH-10;
	y1 = 35;
	y2 = 15;
	y3 = 55;
	z = 35;

	i_nd11 = z+ZLENGTH*(y1+YLENGTH*x1)+1; 
	i_nd12 = z+ZLENGTH*(y1+YLENGTH*x2)+1;

	i_nd21 = z+ZLENGTH*(y2+YLENGTH*x1)+1; 
	i_nd22 = z+ZLENGTH*(y2+YLENGTH*x2)+1;

	i_nd31 = z+ZLENGTH*(y3+YLENGTH*x1)+1; 
	i_nd32 = z+ZLENGTH*(y3+YLENGTH*x2)+1;

	ffr1 = *(mu + i_nd11) / *(mu + i_nd12);
	ffr2 = *(mu + i_nd21) / *(mu + i_nd22);
	ffr3 = *(mu + i_nd31) / *(mu + i_nd32);

	fP3 = fopen("timestep.txt", "a");
	fprintf(fP3,"%d\n", Timek);
	fflush(fP3);
	fclose(fP3);

	fP4 = fopen("ffr1.txt", "a");
	fprintf(fP4,"%lf\n", ffr1);
	fflush(fP4);
	fclose(fP4);

	fP5 = fopen("ffr2.txt", "a");
	fprintf(fP5,"%lf\n", ffr2);
	fflush(fP5);
	fclose(fP5);

	fP6 = fopen("ffr3.txt", "a");
	fprintf(fP6,"%lf\n", ffr3);
	fflush(fP6);
	fclose(fP6);

	return;

}

void	tables	(void)
{
    int x;
    int y;
    int z;
    int i;





    for (x = 1; x <= XLENGTH; x++)
    {
		 for (y = 1; y <= YLENGTH; y++)
      {  for (z = 1; z <= ZLENGTH; z++)
        for (i = 0; i < Q; i++)
	     {   next_x [x][i] = x + cx[i];
		 		if (next_x[x][i] > XLENGTH)
				next_x [x][i] = 1;
	 			if (next_x[x][i] < 1)
				next_x [x][i] = XLENGTH;

             next_y [y][i] = y + cy[i];
		  	  	if (next_y[y][i] > YLENGTH)
				next_y [y][i] = 1;
	 		  	if (next_y[y][i] < 1)
				next_y [y][i] = YLENGTH;
			
             next_z [z][i] = z + cz[i];
		 		if (next_z[z][i] > ZLENGTH)
				next_z [z][i] = 1;
	 			if (next_z[z][i] < 1)
				next_z [z][i] = ZLENGTH;	
 				}
     }
    }
    return;
}

void	initialise (void)
{   
	int        x, y, z, i, x0, y0, z0, a, i_nd, i_in, a1, a2, a3, a4, a5;
   


                  for (i = 0; i < Q; i++)
               {
                   v[0][i] = cx[i];
                   v[1][i] = cy[i];                   
                   v[2][i] = cz[i];                   
               }
   

	for (x = 0; x < XLENGTH; x++)
      for (y = 0; y < YLENGTH; y++)
        for (z = 0; z < ZLENGTH; z++)
      {
           i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;
           *(tau + i_nd) = 1.0;
		   *(real_shear_x + i_nd)=0;
			*(real_shear_y + i_nd)=0;
			*(real_shear_z + i_nd)=0;
			*(real_shear_A + i_nd)=0;
			*(OSI + i_nd)=0;
			*(OSI_x1 + i_nd)=0;
			*(OSI_y1 + i_nd)=0;
			*(OSI_z1 + i_nd)=0;
			*(OSI_x2 + i_nd)=0;
			*(OSI_y2 + i_nd)=0;
			*(OSI_z2 + i_nd)=0;
      }

		for (x = 0; x < XLENGTH; x++) *(RGradient + x + 1) = 1.0;



				  
				  
#pragma omp parallel
		{

			/*
				for (x = 1; x <= XLENGTH; x++)
				for (y = 1; y <= YLENGTH; y++)
				for (z = 1; z <= ZLENGTH; z++)
				{
				i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;
				*(tau + i_nd) = 1.0;
				}
				*/

#pragma omp	for private(i_nd, y, z)
			for (x = 1; x <= XLENGTH; x++)
				for (y = 1; y <= YLENGTH; y++)
					for (z = 1; z <= ZLENGTH; z++)
					{
				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
				*(bnode + i_nd) = 0;
					}

#pragma omp	for private(i_nd, y, z,i,i_in)
			for (x = 1; x <= XLENGTH; x++)
				for (y = 1; y <= YLENGTH; y++)
					for (z = 1; z <= ZLENGTH; z++)
					{
				i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
				//*(bnode + i_nd) = 0.0;

				for (i = 0; i < Q; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					*(b + i_in) = RHO_0 * t[i];
				}

					}
		}

   for (z = 0; z < ZLENGTH; z++) for ( y = 0; y < YLENGTH; y++ ) for ( x = 0; x < XLENGTH; x++ )
   {
			                  
         i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;          
		 *(uu + i_nd) = 0;
   }


    return;
}

void	calc_obs ( void )
{   int         x, y, z, i, i_nd, i_in;
    double      sum1, sum2, sum3, sum4;

#pragma omp parallel
		{
#pragma omp	for private(i_nd, y, z, i, i_in, sum1, sum2, sum3, sum4)
    for (x = 1; x <= XLENGTH; x++)
      for (y = 1; y <= YLENGTH; y++)
        for (z = 1; z <= ZLENGTH; z++)
      {

		     i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;
*(rho + i_nd) =  *(ux + i_nd) = *(uy + i_nd) = *(uz + i_nd) = 0.0;
    
if (*(bnode + i_nd) == 0)
         {
			     sum1 = 0.0;
				 sum2 = 0.0;
				 sum3 = 0.0;
				 sum4 = 0.0;
				 

 	  for (i = 0; i < Q; i++)
          {
			i_in = Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+i+1;
	        sum1	 += *(b + i_in);
            sum2	 += *(b + i_in)*double(cx[i]);
            sum3	 += *(b + i_in)*double(cy[i]);
            sum4	 += *(b + i_in)*double(cz[i]);
	      }
	  *(rho + i_nd) = sum1;
      *(ux + i_nd) = sum2+0.5*(*(FforceX + i_nd) + forcepowerAC);
	  *(uy + i_nd) = sum3+0.5**(FforceY + i_nd);
	  *(uz + i_nd) = sum4+0.5**(FforceZ + i_nd);
         }

         if ( *(rho + i_nd) !=0 )
          {
            *(rhoN + i_nd) = - *(rho + i_nd) / *(rho + i_nd);
	        *(ux + i_nd) /= *(rho + i_nd);
 	        *(uy + i_nd) /= *(rho + i_nd);
 	        *(uz + i_nd) /= *(rho + i_nd);
			*(uu + i_nd) = *(ux + i_nd) * *(ux + i_nd) +  *(uy + i_nd) * *(uy + i_nd) +  *(uz + i_nd) * *(uz + i_nd);
          }




	}
	


		}


		 // outlet velocity
		
		
			x = XLENGTH-1;
			outlet_node = 0;
			outlet_velocity = 0;
			for (y = 1; y <= YLENGTH; y++) for (z = 1; z <= ZLENGTH; z++)
			{
				i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;
				if (*(bnode + i_nd) == 0)
				{
					outlet_node += 1.0;
					outlet_velocity += pow(*(ux + i_nd)* *(ux + i_nd) + *(uy + i_nd)**(uy + i_nd) + *(uz + i_nd)**(uz + i_nd),0.5);

				}
			}
			
		outlet_velocity /= outlet_node;
		
		

    return;
}

void	project ( void )
{   int         x, y, z, i, i_nd, i_in2, i_in, as, bs, x1, x2, i_nd1, i_nd2;
	int			y1, y2, z1, z2;
    double	    udotc, f0, u2, u_mag1, u_mag2, du,LF;

#pragma omp parallel
  {
#pragma omp	for private(y, z, i_in, i_nd, u2, udotc, f0, u_mag1, u_mag2, du,i)
	for (x = 1; x <= XLENGTH; x++)
	  for (y = 1; y <= YLENGTH; y++)
        for (z = 1; z <= ZLENGTH; z++)
           {   
             i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;
     if(*(bnode + i_nd ) != 15)
     	 {
             *(press + i_nd) = *(rho + i_nd)/3.0;
			  *(press_diff + i_nd) = 0;

			  if(x!=1 && x!=XLENGTH)
			 {
				 x1 = x-1;
				 x2 = x+1;
				 i_nd1 = z+ZLENGTH*(y+YLENGTH*x1)+1;
				 i_nd2 = z+ZLENGTH*(y+YLENGTH*x2)+1;

				 if(*(bnode + i_nd ) ==0) if(*(bnode + i_nd1 ) ==0) if(*(bnode + i_nd2 ) ==0)
				 {
					 *(press_diff + i_nd) = (*(press + i_nd2) - *(press + i_nd1))/2;
				 }
			 }



             u2    = (*(ux + i_nd)) * (*(ux + i_nd)) + (*(uy + i_nd)) * (*(uy + i_nd)) + (*(uz + i_nd)) * (*(uz + i_nd));

                if (*(bnode + i_nd) == 0)
                   {	
                  	  for (i = 0; i < Q; i++)
                       {
                         udotc  = (*(ux + i_nd)) * double(cx[i]) + (*(uy + i_nd)) * double(cy[i]) + (*(uz + i_nd)) * double(cz[i]);
			       	     i_in = Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+i+1;
                         *(beq + i_in) =  t[i] * *(rho + i_nd) * (1.0 + 3.0 * udotc + 4.5 * udotc * udotc - 1.5 * u2);
                        }
                    }
            }        
		}

				
		//tau ���

	#pragma omp	for private(y, z, as, bs, i_nd, i_in)
	for (x = 0; x < XLENGTH; x++)
	  for (y = 0; y < YLENGTH; y++)
        for (z = 0; z < ZLENGTH; z++)
           {   
             i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(bnode + i_nd ) != 15)
		  {
               *(summ + i_nd) =0.0;
          for (as = 0; as < 3; as++)
            for (bs = 0; bs < 3; bs++)
              {  
     			i_in = 3*(as+3*(z+ZLENGTH*(y+YLENGTH*x)))+bs+1;          			 
                *(strate + i_in) = 0.0;
              }
          }
		}

	#pragma omp	for private(y, z, i, as, bs, i_in, i_in2, i_nd)
	for (x = 0; x < XLENGTH; x++) for (y = 0; y < YLENGTH; y++) for (z = 0; z < ZLENGTH; z++)
           {   
             i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;
			      if(*(bnode + i_nd ) != 15)
				 {
          for (as = 0; as < 3; as++) for (bs = 0; bs < 3; bs++)
              {  
     			i_in2 = 3*(as+3*(z+ZLENGTH*(y+YLENGTH*x)))+bs+1;          			
					for (i = 0; i < Q; i++)
                   { 
                   i_in = Q*(z+ZLENGTH*(y+YLENGTH*x))+i+1;
                   *(strate + i_in2) = *(strate + i_in2) - 3.0/2.0/ *(tau + i_nd)*(*(b + i_in) - *(beq + i_in))*v[i][as]*v[i][bs] ;   
                    }
                   *(summ + i_nd) = *(summ + i_nd) +(*(strate + i_in2) * *(strate + i_in2)); 
              }
              *(mu + i_nd) = muinf + (mu0 - muinf) * pow((1.0+  (ramda * ramda * 4.0 * *(summ + i_nd))),(seng-1.0)/2.0);
              *(tau + i_nd) = tau0;//3.0 * *(mu + i_nd) / *(rho + i_nd) + 0.5; //tau �� ����
			  *(Shear_stress + i_nd) = pow(*(summ + i_nd),0.5) * 2 * *(mu + i_nd);
			  //*(real_shear_x + i_nd)=0;
			  //*(real_shear_y + i_nd)=0;
			  //*(real_shear_z + i_nd)=0;
			  //*(real_shear_A + i_nd)=0;
			  //*(OSI + i_nd)=0;
            }
		}

	double shear1, shear2;
	double osi_a1, osi_b1, osi_c1, osi_d1;
	double osi_a2, osi_b2, osi_c2, osi_d2;

	#pragma omp	for private(y, z, as, bs, i_nd, i_in)
	for (x = 1; x < XLENGTH-1; x++) for (y = 1; y < YLENGTH-1; y++) for (z = 1; z < ZLENGTH-1; z++)
    {   
		i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;

		if(*(bnode + i_nd ) == 0)
		{
			i_nd1 = z+ZLENGTH*(y+YLENGTH*(x+1))+1;
			i_nd2 = z+ZLENGTH*(y+YLENGTH*(x-1))+1;
			if(*(bnode + i_nd1 ) != 0) shear1 = *(ux + i_nd1) - *(ux + i_nd);
			else shear1 = 0 - *(ux + i_nd);
			if(*(bnode + i_nd2 ) != 0) shear2 = *(ux + i_nd) - *(ux + i_nd2);
			else shear1 =  *(ux + i_nd);
			if((shear1 - shear2)!=0) *(real_shear_x + i_nd) = (shear1 + shear2)/2;
			else *(real_shear_x + i_nd)=0;

			i_nd1 = z+ZLENGTH*((y+1)+YLENGTH*x)+1;
			i_nd2 = z+ZLENGTH*((y-1)+YLENGTH*x)+1;
			if(*(bnode + i_nd1 ) != 0) shear1 = *(uy + i_nd1) - *(uy + i_nd);
			else shear1 = 0 - *(uy + i_nd);
			if(*(bnode + i_nd2 ) != 0) shear2 = *(uy + i_nd) - *(uy + i_nd2);
			else shear1 =  *(uy + i_nd);
			if((shear1 - shear2)!=0) *(real_shear_y + i_nd) = (shear1 + shear2)/2;
			else *(real_shear_y + i_nd)=0;

			i_nd1 = z+1+ZLENGTH*(y+YLENGTH*x)+1;
			i_nd2 = z-1+ZLENGTH*(y+YLENGTH*x)+1;
			if(*(bnode + i_nd1 ) != 0) shear1 = *(uz + i_nd1) - *(uz + i_nd);
			else shear1 = 0 - *(uz + i_nd);
			if(*(bnode + i_nd2 ) != 0) shear2 = *(uz + i_nd) - *(uz + i_nd2);
			else shear1 =  *(uy + i_nd);
			if((shear1 - shear2)!=0) *(real_shear_z + i_nd) = (shear1 + shear2)/2;
			else *(real_shear_z + i_nd)=0;

			*(real_shear_A + i_nd) += *(real_shear_x + i_nd)* *(real_shear_x + i_nd);
			*(real_shear_A + i_nd) += *(real_shear_y + i_nd)* *(real_shear_y + i_nd);
			*(real_shear_A + i_nd) += *(real_shear_z + i_nd)* *(real_shear_z + i_nd);
			//*(real_shear_A + i_nd) *= 10000;
			//*(real_shear_A + i_nd) = sqrt(*(real_shear_A + i_nd));
			//*(real_shear_A + i_nd) /= 100;
			//if(Timek>50000)
			{
				*(OSI_x1 + i_nd) += *(real_shear_x + i_nd);
				*(OSI_y1 + i_nd) += *(real_shear_y + i_nd);
				*(OSI_z1 + i_nd) += *(real_shear_z + i_nd);
				*(OSI_x2 + i_nd) += abs(*(real_shear_x + i_nd));
				*(OSI_y2 + i_nd) += abs(*(real_shear_y + i_nd));
				*(OSI_z2 + i_nd) += abs(*(real_shear_z + i_nd));
				if(*(OSI_x1 + i_nd)!=0) osi_a1 = *(OSI_x1 + i_nd)/(Timek-0); else osi_a1=0;
				if(*(OSI_y1 + i_nd)!=0) osi_b1 = *(OSI_y1 + i_nd)/(Timek-0); else osi_b1=0;
				if(*(OSI_z1 + i_nd)!=0) osi_c1 = *(OSI_z1 + i_nd)/(Timek-0); else osi_c1=0;
				if(*(OSI_x2 + i_nd)!=0) osi_a2 = *(OSI_x2 + i_nd)/(Timek-0); else osi_a2=0;
				if(*(OSI_y2 + i_nd)!=0) osi_b2 = *(OSI_y2 + i_nd)/(Timek-0); else osi_b2=0;
				if(*(OSI_z2 + i_nd)!=0) osi_c2 = *(OSI_z2 + i_nd)/(Timek-0); else osi_c2=0;
				osi_d1 = pow((osi_a1*osi_a1)+(osi_b1*osi_b1)+(osi_c1*osi_c1),0.5);
				osi_d2 = pow((osi_a2*osi_a2)+(osi_b2*osi_b2)+(osi_c2*osi_c2),0.5);
				if(osi_d2!=0) *(OSI + i_nd) = osi_d1 / osi_d2; else *(OSI + i_nd)=0;
				//*(OSI + i_nd) = osi_d1;
				 //*(OSI + i_nd) = (1-*(OSI + i_nd))*0.5;
			}



		}





	}

	




    #pragma omp	for private(y, z, i_in, i_nd, i,LF)
	for (x = 1; x <= XLENGTH; x++)
	  for (y = 1; y <= YLENGTH; y++)
        for (z = 1; z <= ZLENGTH; z++)
        {   i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;                         
          for (i = 0; i < Q; i++)
               {
                   if (*(bnode + i_nd) == 0)
                      { 
							*(tau + i_nd) = tau0; //Newtonian
                      		i_in = Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+i+1;
							LF = (1.0 - 0.5* *(tau + i_nd))*t[i] * ((3.0 * (double(cx[i]) - *(ux + i_nd)) + 9.0 * (double(cx[i]) * *(ux + i_nd) + double(cy[i]) * *(uy + i_nd) + double(cz[i]) * *(uz + i_nd)) * double(cx[i])) * ((*(FforceX + i_nd) + forcepowerAC)) + (3.0 * (double(cy[i]) - *(uy + i_nd)) + 9.0 * (double(cx[i]) * *(ux + i_nd) + double(cy[i]) * *(uy + i_nd) + double(cz[i]) * *(uz + i_nd)) * double(cy[i])) * *(FforceY + i_nd) + (3.0 * (double(cz[i]) - *(uz + i_nd)) + 9.0 * (double(cx[i]) * *(ux + i_nd) + double(cy[i]) * *(uy + i_nd) + double(cz[i]) * *(uz + i_nd)) * double(cz[i])) * *(FforceZ + i_nd));
					
							*(b + i_in) = (1.0-1.0 / *(tau + i_nd)) * *(b + i_in) + 1.0 / *(tau + i_nd) * *(beq + i_in) + LF;
							
                      }
                }
        } 



   }

  
     return;

}

void	propagate(void)
{
	int	x, y, z, i, j, new_x, new_y, new_z, l, h, m, i_in, i_in2, i_in3, i_in4, i_nd;
	double aux, aux1, aux2, rho0, ru, u_Poi, out_pressure, yDiff, zDiff, rr, rrr, xxx, yyy,zzz, u_p;
	double out_pressure2;
	int newstep, newstep2;

	
	

	out_pressure = 1.1;

#pragma omp parallel
	{
#pragma omp	for private(y, z, i, new_x, new_y, new_z,i_in,i_in2,i_in3)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		for (i = 0; i < Q; i++)
		{
			i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;

			i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
			if (*(bnode + i_in2) != 15)
			{
				new_x = next_x[x][i];
				new_y = next_y[y][i];
				new_z = next_z[z][i];
				i_in3 = Q*(new_z - 1 + ZLENGTH*(new_y - 1 + YLENGTH*(new_x - 1))) + i + 1;
				*(back + i_in3) = *(b + i_in);

			}
		}

#pragma omp	for private(y, z, i_in,i_nd,i)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		{
			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
			if (*(bnode + i_nd) != 15)
			{

				for (i = 0; i < Q; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					*(b + i_in) = *(back + i_in);
				}
			}
		}

#pragma omp	for private(y, z, aux, i, i_nd, i_in)   
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z <= ZLENGTH; z++)
		{
			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
			if (*(bnode + i_nd) == 1)
			{
				for (i = 1; i <= 1; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					aux = *(b + i_in);
					*(b + i_in) = *(b + i_in + 9);
					*(b + i_in + 9) = aux;
				}

				for (i = 2; i <= 2; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					aux = *(b + i_in);
					*(b + i_in) = *(b + i_in + 9);
					*(b + i_in + 9) = aux;
				}

				for (i = 3; i <= 3; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					aux = *(b + i_in);
					*(b + i_in) = *(b + i_in + 9);
					*(b + i_in + 9) = aux;
				}

				for (i = 4; i <= 4; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					aux = *(b + i_in);
					*(b + i_in) = *(b + i_in + 9);
					*(b + i_in + 9) = aux;
				}

				for (i = 5; i <= 5; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					aux = *(b + i_in);
					*(b + i_in) = *(b + i_in + 9);
					*(b + i_in + 9) = aux;
				}

				for (i = 6; i <= 6; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					aux = *(b + i_in);
					*(b + i_in) = *(b + i_in + 9);
					*(b + i_in + 9) = aux;
				}

				for (i = 7; i <= 7; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					aux = *(b + i_in);
					*(b + i_in) = *(b + i_in + 9);
					*(b + i_in + 9) = aux;
				}

				for (i = 8; i <= 8; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					aux = *(b + i_in);
					*(b + i_in) = *(b + i_in + 9);
					*(b + i_in + 9) = aux;
				}

				for (i = 9; i <= 9; i++)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					aux = *(b + i_in);
					*(b + i_in) = *(b + i_in + 9);
					*(b + i_in + 9) = aux;
				}

			}
		}
		/*
#pragma omp	for private(y, i_in, i_nd, i)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		{
			i_nd = 1 + ZLENGTH*(y + YLENGTH*x) + 1;
			for (i = 0; i < Q; i++)
			{
				if (*(bnode + i_nd) == 0)
				{

					i_in = Q*(1 - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					*(b + i_in) = t[i] * (1.0 + 3.0*cx[i] * uMax);

				}
			}
		}
#pragma omp	for private(y, z, i_in, i_nd, i)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = ZLENGTH; z < ZLENGTH + 1; z++)
		{
			i_nd = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
			for (i = 0; i < Q; i++)
			{
				if (*(bnode + i_nd) == 0)
				{
					i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
					*(b + i_in) = t[i] * (1.0 - 3.0*cx[i] * uMax);

				}
			}
		}
		*/
		
		/*
		#pragma omp	for private(y, z, i_in, i_nd, i)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = 1; z < 2; z++)
		{   i_nd = z+ZLENGTH*(y+YLENGTH*x)+1;
		for (i = 0; i < Q; i++)
		{
		if (*(bnode + i_nd) == 0)
		{

		i_in = Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+i+1;
		*(b + i_in) = t[i]*(1.0+3.0*cx[i]*uMax);

		}
		}
		}
		#pragma omp	for private(y, z, i_in, i_nd, i)
		for (x = 1; x <= XLENGTH; x++)
		for (y = 1; y <= YLENGTH; y++)
		for (z = ZLENGTH; z < ZLENGTH+1; z++)
		{   i_nd = z-1+ZLENGTH*(y-1+YLENGTH*(x-1))+1;
		for (i = 0; i < Q; i++)
		{
		if (*(bnode + i_nd) == 0)
		{
		i_in = Q*(z-1+ZLENGTH*(y-1+YLENGTH*(x-1)))+i+1;
		*(b + i_in) = t[i]*(1.0-3.0*cx[i]*uMax);

		}
		}
		}

		*/


// inlet -> outlet 바꿔야 함
#pragma omp	for private(z, x, i, i_in, i_in2, l, h, m, u_Poi, rho0, ru,yDiff,zDiff, rr, rrr)          
		for (y = 1; y <= YLENGTH; y++)
		{
			for (z = 1; z <= ZLENGTH; z++)
			{
				for (x = 1; x <= 1; x++)
				{
					i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

					if (*(bnode + i_in2) == 0)
					{


						for (i = 0; i < 1; i++)
						{
							i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;


							rrr = (double(z) - (double(ZLENGTH) + 1.0) / 2.0)*(double(z) - (double(ZLENGTH) + 1.0) / 2.0) + (double(y) - (double(YLENGTH) + 1.0) / 2.0)*(double(y) - (double(YLENGTH) + 1.0) / 2.0);
							rr = pow(rrr, 0.5);

							 //*(ux + i_in2) = uMax*(1.0 + (rr / pipe_r)* (rr / pipe_r));
							//*(ux + i_in2) = uMax;
							*(ux + i_in2) = uMax*pulse[pulse_step];

							rho0 = (*(b + i_in) + *(b + i_in + 16) + *(b + i_in + 7) + *(b + i_in + 9) + *(b + i_in + 18) + *(b + i_in + 17) + *(b + i_in + 8) + *(b + i_in + 15) + *(b + i_in + 6) + 2.0 *(*(b + i_in + 3) + *(b + i_in + 1) + *(b + i_in + 5) + *(b + i_in + 4) + *(b + i_in + 2))) / (1.0 + *(ux + i_in2));

							yDiff = *(b + i_in + 15) + *(b + i_in + 17) + *(b + i_in + 16) - (*(b + i_in + 6) + *(b + i_in + 7) + *(b + i_in + 8));
							zDiff = *(b + i_in + 8) + *(b + i_in + 15) + *(b + i_in + 18) - (*(b + i_in + 6) + *(b + i_in + 9) + *(b + i_in + 17));

							*(b + i_in + 12) = *(b + i_in + 3) + (1.0 / 3.0)**(ux + i_in2)*rho0;
							*(b + i_in + 10) = *(b + i_in + 1) + (1.0 / 6.0)*rho0**(ux + i_in2) - (1.0 / 2.0) * (yDiff);
							*(b + i_in + 14) = *(b + i_in + 5) + (1.0 / 6.0)*rho0**(ux + i_in2) + (1.0 / 2.0) * (yDiff);
							*(b + i_in + 13) = *(b + i_in + 4) + (1.0 / 6.0)*rho0**(ux + i_in2) + (1.0 / 2.0) * (zDiff);
							*(b + i_in + 11) = *(b + i_in + 2) + (1.0 / 6.0)*rho0**(ux + i_in2) - (1.0 / 2.0) * (zDiff);
						}
					}
				}
			}
		}









		/*
#pragma omp	for private(z, x, i, i_in, i_in2, i_in3, i_in4, l, h, m, u_Poi, rho0, ru,yDiff,zDiff, rr, rrr)          
		for (y = 1; y <= YLENGTH; y++)
		{
			for (z = 1; z <= ZLENGTH; z++)
			{
				for (x = XLENGTH; x <= XLENGTH; x++)
				{
					i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

					if (*(bnode + i_in2) == 0)
					{


						for (i = 0; i < 1; i++)
						{


							i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
							i_in3 = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 2))) + i + 1;
							i_in4 = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 3))) + i + 1;

							*(b + i_in + 1) = 2.0 * *(b + i_in3 + 1) - *(b + i_in4 + 1);
							*(b + i_in + 2) = 2.0 * *(b + i_in3 + 2) - *(b + i_in4 + 2);
							*(b + i_in + 3) = 2.0 * *(b + i_in3 + 3) - *(b + i_in4 + 3);
							*(b + i_in + 4) = 2.0 * *(b + i_in3 + 4) - *(b + i_in4 + 4);
							*(b + i_in + 5) = 2.0 * *(b + i_in3 + 5) - *(b + i_in4 + 5);
						}
					}
				}
			}
		}
		*/



// outlet?
#pragma omp	for private(z, x, i, i_in, i_in2, i_in3, i_in4, l, h, m, u_Poi, rho0, ru,yDiff,zDiff, rr, rrr, pressure1)          
		for (y = 1; y <= YLENGTH; y++)
		{
			for (z = 1; z <= ZLENGTH; z++)
			{
				//for (x = XLENGTH; x <= XLENGTH; x++)
				for (x = 1; x <= XLENGTH; x++)
				{
					i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;
					Resistance = *(uu + i_in3)*0.3;//0.3;
					
					if (*(b_b + i_in2) == 7)
					{
						i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + 1;
						i_in3= Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 2))) + 1;
						*(ux + i_in3) = *(b + i_in3+10)+*(b + i_in3+11)+*(b + i_in3+12)+*(b + i_in3+13)+*(b + i_in3+14)-(*(b + i_in3+1)+*(b + i_in3+2)+*(b + i_in3+3)+*(b + i_in3+4)+*(b + i_in3+5));
						
						//pressure1 = 1.0 + Resistance;
						pressure1 = 1.2 + Resistance;
						u_p = 2*(*(b + i_in+10)+*(b + i_in+11)+*(b + i_in+12)+*(b + i_in+13)+*(b + i_in+14))+(*(b + i_in+0)+*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8)+*(b + i_in+9)+*(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)+*(b + i_in+18))-pressure1;
						yDiff = *(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)-(*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8));
			            zDiff = *(b + i_in+8)+*(b + i_in+15)+*(b + i_in+18)-(*(b + i_in+6)+*(b + i_in+9)+*(b + i_in+17)); 
						*(b + i_in+3)=*(b + i_in+12)-(1.0/3.0)*u_p;
						*(b + i_in+1)=*(b + i_in+10)-(1.0/6.0)*u_p  + (1.0/2.0) * (yDiff); 
						*(b + i_in+5)=*(b + i_in+14)-(1.0/6.0)*u_p  - (1.0/2.0) * (yDiff);  
						*(b + i_in+4)=*(b + i_in+13)-(1.0/6.0)*u_p  - (1.0/2.0) * (zDiff);   
						*(b + i_in+2)=*(b + i_in+11)-(1.0/6.0)*u_p  + (1.0/2.0) * (zDiff); 
						
						//new model // *(ux + i_in2) = -1.0 + (*(b + i_in+ 6)+ *(b + i_in+7)+*(b + i_in+8)+*(b + i_in+16)+*(b + i_in+0)+*(b + i_in+9)+*(b + i_in+17)+*(b + i_in+18)+*(b + i_in+15) + 2.0*(*(b + i_in+10)+*(b + i_in+11)+*(b + i_in+12)+*(b + i_in+13)+*(b + i_in+14)))/pressure1;
					}

					else if (*(b_b + i_in2) == 1)
					{
						i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + 1;
						i_in3= Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 2))) + 1;
						*(ux + i_in3) = *(b + i_in3+10)+*(b + i_in3+11)+*(b + i_in3+12)+*(b + i_in3+13)+*(b + i_in3+14)-(*(b + i_in3+1)+*(b + i_in3+2)+*(b + i_in3+3)+*(b + i_in3+4)+*(b + i_in3+5));
						
						//pressure1 = 1.0 + Resistance;
						pressure1 = 1.0 + Resistance;
						u_p = 2*(*(b + i_in+10)+*(b + i_in+11)+*(b + i_in+12)+*(b + i_in+13)+*(b + i_in+14))+(*(b + i_in+0)+*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8)+*(b + i_in+9)+*(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)+*(b + i_in+18))-pressure1;
						yDiff = *(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)-(*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8));
			            zDiff = *(b + i_in+8)+*(b + i_in+15)+*(b + i_in+18)-(*(b + i_in+6)+*(b + i_in+9)+*(b + i_in+17)); 
						*(b + i_in+3)=*(b + i_in+12)-(1.0/3.0)*u_p;
						*(b + i_in+1)=*(b + i_in+10)-(1.0/6.0)*u_p  + (1.0/2.0) * (yDiff); 
						*(b + i_in+5)=*(b + i_in+14)-(1.0/6.0)*u_p  - (1.0/2.0) * (yDiff);  
						*(b + i_in+4)=*(b + i_in+13)-(1.0/6.0)*u_p  - (1.0/2.0) * (zDiff);   
						*(b + i_in+2)=*(b + i_in+11)-(1.0/6.0)*u_p  + (1.0/2.0) * (zDiff); 
						
						//new model // *(ux + i_in2) = -1.0 + (*(b + i_in+ 6)+ *(b + i_in+7)+*(b + i_in+8)+*(b + i_in+16)+*(b + i_in+0)+*(b + i_in+9)+*(b + i_in+17)+*(b + i_in+18)+*(b + i_in+15) + 2.0*(*(b + i_in+10)+*(b + i_in+11)+*(b + i_in+12)+*(b + i_in+13)+*(b + i_in+14)))/pressure1;
					}

					else if (*(b_b + i_in2) == 2) //x-
					{
						i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + 1;
						i_in3= Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 0))) + 1;
						*(ux + i_in3) = *(b + i_in3+10)+*(b + i_in3+11)+*(b + i_in3+12)+*(b + i_in3+13)+*(b + i_in3+14)-(*(b + i_in3+1)+*(b + i_in3+2)+*(b + i_in3+3)+*(b + i_in3+4)+*(b + i_in3+5));
						*(ux + i_in3) = *(ux + i_in3)*(-1);
						pressure1 = 1.0 + Resistance;
						u_p = 2*(*(b + i_in+1)+*(b + i_in+2)+*(b + i_in+3)+*(b + i_in+4)+*(b + i_in+5))+(*(b + i_in+0)+*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8)+*(b + i_in+9)+*(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)+*(b + i_in+18))-pressure1;
						yDiff = *(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)-(*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8));
			            zDiff = *(b + i_in+8)+*(b + i_in+15)+*(b + i_in+18)-(*(b + i_in+6)+*(b + i_in+9)+*(b + i_in+17)); 
						*(b + i_in+12)=*(b + i_in+3)-(1.0/3.0)*u_p;
						*(b + i_in+10)=*(b + i_in+1)-(1.0/6.0)*u_p  - (1.0/2.0) * (yDiff); 
						*(b + i_in+14)=*(b + i_in+5)-(1.0/6.0)*u_p  + (1.0/2.0) * (yDiff);  
						*(b + i_in+13)=*(b + i_in+4)-(1.0/6.0)*u_p  + (1.0/2.0) * (zDiff);   
						*(b + i_in+11)=*(b + i_in+2)-(1.0/6.0)*u_p  - (1.0/2.0) * (zDiff);
					}

					else if (*(b_b + i_in2) == 3) //y-
					{
						i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x))) + 1;
						i_in3= Q*(z - 1 + ZLENGTH*(y - 0 + YLENGTH*(x))) + 1;
						*(ux + i_in3) = *(b + i_in3+5)+*(b + i_in3+10)+*(b + i_in3+15)+*(b + i_in3+16)+*(b + i_in3+17)-(*(b + i_in3+1)+*(b + i_in3+6)+*(b + i_in3+7)+*(b + i_in3+8)+*(b + i_in3+14));
						*(ux + i_in3) = *(ux + i_in3)*(-1);
						pressure1 = 1.0 + Resistance;
						u_p = 2*(*(b + i_in+1)+*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8)+*(b + i_in+14))+(*(b + i_in+0)+*(b + i_in+2)+*(b + i_in+9)+*(b + i_in+13)+*(b + i_in+3)+*(b + i_in+12)+*(b + i_in+4)+*(b + i_in+18)+*(b + i_in+11))-pressure1;
						yDiff = *(b + i_in+11)+*(b + i_in+12)+*(b + i_in+13)-(*(b + i_in+2)+*(b + i_in+3)+*(b + i_in+4));
			            zDiff = *(b + i_in+8)+*(b + i_in+15)+*(b + i_in+18)-(*(b + i_in+6)+*(b + i_in+9)+*(b + i_in+17)); 
						*(b + i_in+16)=*(b + i_in+7)-(1.0/3.0)*u_p;
						*(b + i_in+10)=*(b + i_in+1)-(1.0/6.0)*u_p  - (1.0/2.0) * (yDiff); 
						*(b + i_in+5)=*(b + i_in+14)-(1.0/6.0)*u_p  + (1.0/2.0) * (yDiff);  
						*(b + i_in+17)=*(b + i_in+8)-(1.0/6.0)*u_p  + (1.0/2.0) * (zDiff);   
						*(b + i_in+15)=*(b + i_in+6)-(1.0/6.0)*u_p  - (1.0/2.0) * (zDiff);
					}


					else if (*(b_b + i_in2) == 4) //z+
					{
						i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x))) + 1;
						i_in3= Q*(z - 2 + ZLENGTH*(y - 1 + YLENGTH*(x))) + 1;
						*(ux + i_in3) = *(b + i_in3+4)+*(b + i_in3+8)+*(b + i_in3+11)+*(b + i_in3+15)+*(b + i_in3+18)-(*(b + i_in3+2)+*(b + i_in3+6)+*(b + i_in3+9)+*(b + i_in3+13)+*(b + i_in3+17));
						pressure1 = 1.0 + Resistance;
						u_p = 2*(*(b + i_in+10)+*(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)+*(b + i_in+5))+(*(b + i_in+0)+*(b + i_in+2)+*(b + i_in+9)+*(b + i_in+13)+*(b + i_in+3)+*(b + i_in+12)+*(b + i_in+4)+*(b + i_in+18)+*(b + i_in+11))-pressure1;
						yDiff = *(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)-(*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8));
			            zDiff = *(b + i_in+10)+*(b + i_in+12)+*(b + i_in+14)-(*(b + i_in+1)+*(b + i_in+3)+*(b + i_in+5)); 

						*(b + i_in+9)=*(b + i_in+18)-(1.0/3.0)*u_p;
						*(b + i_in+6)=*(b + i_in+15)-(1.0/6.0)*u_p  + (1.0/2.0) * (yDiff); 
						*(b + i_in+17)=*(b + i_in+8)-(1.0/6.0)*u_p  - (1.0/2.0) * (yDiff);  
						*(b + i_in+13)=*(b + i_in+4)-(1.0/6.0)*u_p  - (1.0/2.0) * (zDiff);   
						*(b + i_in+2)=*(b + i_in+11)-(1.0/6.0)*u_p  + (1.0/2.0) * (zDiff); 

					}


					else if (*(b_b + i_in2) == 5) //y+
					{
						i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x))) + 1;
						i_in3= Q*(z - 1 + ZLENGTH*(y - 2 + YLENGTH*(x))) + 1;
						*(ux + i_in3) = *(b + i_in3+5)+*(b + i_in3+10)+*(b + i_in3+15)+*(b + i_in3+16)+*(b + i_in3+17)-(*(b + i_in3+1)+*(b + i_in3+6)+*(b + i_in3+7)+*(b + i_in3+8)+*(b + i_in3+14));
						pressure1 = 1.0 + Resistance;
						u_p = 2*(*(b + i_in+4)+*(b + i_in+8)+*(b + i_in+11)+*(b + i_in+15)+*(b + i_in+18))+(*(b + i_in+0)+*(b + i_in+1)+*(b + i_in+3)+*(b + i_in+5)+*(b + i_in+7)+*(b + i_in+10)+*(b + i_in+12)+*(b + i_in+14)+*(b + i_in+16))-pressure1;
						yDiff = *(b + i_in+11)+*(b + i_in+12)+*(b + i_in+13)-(*(b + i_in+2)+*(b + i_in+3)+*(b + i_in+4));
			            zDiff = *(b + i_in+8)+*(b + i_in+15)+*(b + i_in+18)-(*(b + i_in+6)+*(b + i_in+9)+*(b + i_in+17)); 
						*(b + i_in+7)=*(b + i_in+16)-(1.0/3.0)*u_p;
						*(b + i_in+1)=*(b + i_in+10)-(1.0/6.0)*u_p  + (1.0/2.0) * (yDiff); 
						*(b + i_in+14)=*(b + i_in+5)-(1.0/6.0)*u_p  - (1.0/2.0) * (yDiff);  
						*(b + i_in+8)=*(b + i_in+17)-(1.0/6.0)*u_p  - (1.0/2.0) * (zDiff);   
						*(b + i_in+6)=*(b + i_in+15)-(1.0/6.0)*u_p  + (1.0/2.0) * (zDiff); 

					}

//pressure Outlet - z0 direction
						/*
						i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + 1;
						i_in3= Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 2))) + 1;
						*(ux + i_in3) = *(b + i_in3+4)+*(b + i_in3+15)+*(b + i_in3+18)+*(b + i_in3+8)+*(b + i_in3+11)-(*(b + i_in3+2)+*(b + i_in3+6)+*(b + i_in3+9)+*(b + i_in3+17)+*(b + i_in3+13));
						*(ux + i_in3) = *(ux + i_in3)*(-1);
						pressure1 = 1.0 + (outlet_velocity) * (Resistance);
						u_p = 2*(*(b + i_in+4)+*(b + i_in+15)+*(b + i_in+18)+*(b + i_in+8)+*(b + i_in+11))+(*(b + i_in+0)+*(b + i_in+1)+*(b + i_in+3)+*(b + i_in+5)+*(b + i_in+7)+*(b + i_in+16)+*(b + i_in+10)+*(b + i_in+12)+*(b + i_in+14))-pressure1;
						yDiff = *(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)-(*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8));
			            zDiff = *(b + i_in+10)+*(b + i_in+12)+*(b + i_in+14)-(*(b + i_in+1)+*(b + i_in+3)+*(b + i_in+5)); 
						*(b + i_in+18)=*(b + i_in+9)-(1.0/3.0)*u_p;
						*(b + i_in+15)=*(b + i_in+6)-(1.0/6.0)*u_p  - (1.0/2.0) * (yDiff); 
						*(b + i_in+8)=*(b + i_in+17)-(1.0/6.0)*u_p  + (1.0/2.0) * (yDiff);  
						*(b + i_in+4)=*(b + i_in+13)-(1.0/6.0)*u_p  + (1.0/2.0) * (zDiff);   
						*(b + i_in+11)=*(b + i_in+2)-(1.0/6.0)*u_p  - (1.0/2.0) * (zDiff); 
						*/
						//pressure Outlet - zmax -> z�� ���̳ʺ� ������ ������ ��
						/*
						i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + 1;
						i_in3= Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 2))) + 1;
						*(ux + i_in3) = *(b + i_in3+4)+*(b + i_in3+15)+*(b + i_in3+18)+*(b + i_in3+8)+*(b + i_in3+11)-(*(b + i_in3+2)+*(b + i_in3+6)+*(b + i_in3+9)+*(b + i_in3+17)+*(b + i_in3+13));
						pressure1 = 1.0 + (outlet_velocity) * (Resistance);
						u_p = 2*(*(b + i_in+2)+*(b + i_in+6)+*(b + i_in+9)+*(b + i_in+17)+*(b + i_in+13))+(*(b + i_in+0)+*(b + i_in+1)+*(b + i_in+3)+*(b + i_in+5)+*(b + i_in+7)+*(b + i_in+16)+*(b + i_in+10)+*(b + i_in+12)+*(b + i_in+14))-pressure1;
						yDiff = *(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)-(*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8));
			            zDiff = *(b + i_in+10)+*(b + i_in+12)+*(b + i_in+14)-(*(b + i_in+1)+*(b + i_in+3)+*(b + i_in+5)); 
						*(b + i_in+9)=*(b + i_in+18)-(1.0/3.0)*u_p;
						*(b + i_in+6)=*(b + i_in+15)-(1.0/6.0)*u_p  + (1.0/2.0) * (yDiff); 
						*(b + i_in+17)=*(b + i_in+8)-(1.0/6.0)*u_p  - (1.0/2.0) * (yDiff);  
						*(b + i_in+13)=*(b + i_in+4)-(1.0/6.0)*u_p  - (1.0/2.0) * (zDiff);   
						*(b + i_in+2)=*(b + i_in+11)-(1.0/6.0)*u_p  + (1.0/2.0) * (zDiff); 
						*/

						//pressure Outlet - x+ direction
						/*
						i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + 1;
						i_in3= Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 2))) + 1;
						*(ux + i_in3) = *(b + i_in3+10)+*(b + i_in3+11)+*(b + i_in3+12)+*(b + i_in3+13)+*(b + i_in3+14)-(*(b + i_in3+1)+*(b + i_in3+2)+*(b + i_in3+3)+*(b + i_in3+4)+*(b + i_in3+5));
						pressure1 = 1.0 + (outlet_velocity) * (Resistance);
						u_p = 2*(*(b + i_in+10)+*(b + i_in+11)+*(b + i_in+12)+*(b + i_in+13)+*(b + i_in+14))+(*(b + i_in+0)+*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8)+*(b + i_in+9)+*(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)+*(b + i_in+18))-pressure1;
						yDiff = *(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)-(*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8));
			            zDiff = *(b + i_in+8)+*(b + i_in+15)+*(b + i_in+18)-(*(b + i_in+6)+*(b + i_in+9)+*(b + i_in+17)); 
						*(b + i_in+3)=*(b + i_in+12)-(1.0/3.0)*u_p;
						*(b + i_in+1)=*(b + i_in+10)-(1.0/6.0)*u_p  + (1.0/2.0) * (yDiff); 
						*(b + i_in+5)=*(b + i_in+14)-(1.0/6.0)*u_p  - (1.0/2.0) * (yDiff);  
						*(b + i_in+4)=*(b + i_in+13)-(1.0/6.0)*u_p  - (1.0/2.0) * (zDiff);   
						*(b + i_in+2)=*(b + i_in+11)-(1.0/6.0)*u_p  + (1.0/2.0) * (zDiff); 
						*/

						//pressure Outlet - x- direction
						/*
						i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + 1;
						i_in3= Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 2))) + 1;
						*(ux + i_in3) = *(b + i_in3+10)+*(b + i_in3+11)+*(b + i_in3+12)+*(b + i_in3+13)+*(b + i_in3+14)-(*(b + i_in3+1)+*(b + i_in3+2)+*(b + i_in3+3)+*(b + i_in3+4)+*(b + i_in3+5));
						*(ux + i_in3) = *(ux + i_in3)*(-1);
						pressure1 = 1.0 + (outlet_velocity) * (Resistance);
						u_p = 2*(*(b + i_in+1)+*(b + i_in+2)+*(b + i_in+3)+*(b + i_in+4)+*(b + i_in+5))+(*(b + i_in+0)+*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8)+*(b + i_in+9)+*(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)+*(b + i_in+18))-pressure1;
						yDiff = *(b + i_in+15)+*(b + i_in+16)+*(b + i_in+17)-(*(b + i_in+6)+*(b + i_in+7)+*(b + i_in+8));
			            zDiff = *(b + i_in+8)+*(b + i_in+15)+*(b + i_in+18)-(*(b + i_in+6)+*(b + i_in+9)+*(b + i_in+17)); 
						*(b + i_in+12)=*(b + i_in+3)-(1.0/3.0)*u_p;
						*(b + i_in+10)=*(b + i_in+1)-(1.0/6.0)*u_p  - (1.0/2.0) * (yDiff); 
						*(b + i_in+14)=*(b + i_in+5)-(1.0/6.0)*u_p  + (1.0/2.0) * (yDiff);  
						*(b + i_in+13)=*(b + i_in+4)-(1.0/6.0)*u_p  + (1.0/2.0) * (zDiff);   
						*(b + i_in+11)=*(b + i_in+2)-(1.0/6.0)*u_p  - (1.0/2.0) * (zDiff); 
						*/


				}
			}
		}


		

		/*
		#pragma omp	for private(z, x, i, i_in, i_in2, l, h, m, u_Poi, rho0, ru,yDiff,zDiff, rr, rrr)
		for (y = 51; y <= YLENGTH; y++)
		{
		for (z = 1; z <= ZLENGTH; z++)
		{
		for (x = 1; x <= 1; x++)
		{
		i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

		if (*(bnode + i_in2) == 0)
		{


		for (i = 0; i < 1; i++)
		{
		i_in = Q*(z+ZLENGTH*(y+YLENGTH*x))+i+1;
		l = 2; h = ZLENGTH-2; m = (l+h)/2.0;
		u_Poi =  u_X0/m/m*(z-l)*(z-h);
		rho0 =  (*(b + i_in) + *(b + i_in + 16) + *(b + i_in + 7) + *(b + i_in + 9) + *(b + i_in+ 18) + *(b + i_in + 17) + *(b + i_in + 8) + *(b + i_in + 15) + *(b + i_in + 6) + 2.0 *(*(b + i_in + 3) + *(b + i_in + 1) + *(b + i_in + 5) + *(b + i_in + 4) + *(b + i_in + 2)))/(1.0 + u_Poi);
		ru =  rho0 *  u_Poi;
		*(b + i_in+12)= *(b + i_in+3) - (1.0/3.0)*ru;
		*(b + i_in+10)= *(b + i_in+1) - (1.0/6.0)*ru + (1.0/4.0) * (*(b + i_in+7) - *(b + i_in+16));
		*(b + i_in+14)= *(b + i_in+5) - (1.0/6.0)*ru + (1.0/4.0) * (*(b + i_in+16)  - *(b + i_in+7));
		*(b + i_in+13)= *(b + i_in+4) - (1.0/6.0)*ru + (1.0/4.0) * (*(b + i_in+18) - *(b + i_in+9));
		*(b + i_in+11)= *(b + i_in+2) - (1.0/6.0)*ru + (1.0/4.0) * (*(b + i_in+9)  - *(b + i_in+18));


		i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;



		rrr = (z - 13.0)*(z - 13.0) + (y - 85)*(y - 85);
		rr = pow(rrr, 0.5);

		*(ux + i_in2) = 0.2* uMax*(1.0 - (rr / pipe_r)* (rr / pipe_r));

		rho0 = (*(b + i_in) + *(b + i_in + 16) + *(b + i_in + 7) + *(b + i_in + 9) + *(b + i_in + 18) + *(b + i_in + 17) + *(b + i_in + 8) + *(b + i_in + 15) + *(b + i_in + 6) + 2.0 *(*(b + i_in + 3) + *(b + i_in + 1) + *(b + i_in + 5) + *(b + i_in + 4) + *(b + i_in + 2))) / (1.0 + *(ux + i_in2));

		yDiff = *(b + i_in + 15) + *(b + i_in + 17) + *(b + i_in + 18) - (*(b + i_in + 6) + *(b + i_in + 7) + *(b + i_in + 8));
		zDiff = *(b + i_in + 8) + *(b + i_in + 15) + *(b + i_in + 16) - (*(b + i_in + 6) + *(b + i_in + 9) + *(b + i_in + 17));

		*(b + i_in + 12) = *(b + i_in + 3) + (1.0 / 3.0)**(ux + i_in2)*rho0;
		*(b + i_in + 10) = *(b + i_in + 1) + (1.0 / 6.0)*rho0**(ux + i_in2) - (1.0 / 2.0) * (yDiff);
		*(b + i_in + 14) = *(b + i_in + 5) + (1.0 / 6.0)*rho0**(ux + i_in2) + (1.0 / 2.0) * (yDiff);
		*(b + i_in + 13) = *(b + i_in + 4) + (1.0 / 6.0)*rho0**(ux + i_in2) + (1.0 / 2.0) * (zDiff);
		*(b + i_in + 11) = *(b + i_in + 2) + (1.0 / 6.0)*rho0**(ux + i_in2) - (1.0 / 2.0) * (zDiff);
		}
		}
		}
		}
		}

	
#pragma omp	for private(z, x, i, i_in, i_in2, i_in3, i_in4, l, h, m, u_Poi, rho0, ru,yDiff,zDiff, rr, rrr)          
		for (y = 1; y <= YLENGTH; y++)
		{
			for (z = ZLENGTH; z <= ZLENGTH; z++)
			{
				for (x = 1; x <= XLENGTH; x++)
				{
					i_in2 = z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1)) + 1;

					if (*(bnode + i_in2) == 0)
					{


						for (i = 0; i < 1; i++)
						{


							i_in = Q*(z - 1 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
							i_in3 = Q*(z - 2 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;
							i_in4 = Q*(z - 3 + ZLENGTH*(y - 1 + YLENGTH*(x - 1))) + i + 1;

							*(b + i_in + 2) = 2.0 * *(b + i_in3 + 2) - *(b + i_in4 + 2);
							*(b + i_in + 6) = 2.0 * *(b + i_in3 + 6) - *(b + i_in4 + 6);
							*(b + i_in + 9) = 2.0 * *(b + i_in3 + 9) - *(b + i_in4 + 9);
							*(b + i_in + 13) = 2.0 * *(b + i_in3 + 13) - *(b + i_in4 + 13);
							*(b + i_in + 17) = 2.0 * *(b + i_in3 + 17) - *(b + i_in4 + 17);
						}
					}
				}
			}
		}

		*/

	} //End of omp


	return;

}	/* End of procedure - propagate	*/

void	constructimage(void)
{
	FILE *f3;
	FILE *f4;
	char buffer[50];
	int x1,y1,z1;
	int x, y, z;
	int x2,y2,z2;
	int xp, yp, zp;
	int p1, p2, p3;
	int i_in1, i_in2, i;
	int A1, A2;
	int mod1, mod2;
	
	
	if((f3 = fopen("coronary.plt", "r")) == NULL)
	 {
		 printf("sibal");
		 system("PAUSE");
		 return;
	}

	f4 = fopen("test1.plt", "w");
	fprintf(f4, "Variables = X, Y, Z, boundary\n");
	fprintf(f4, "Zone I= %d,J= %d,K= %d, F=POINT\n", XLENGTH, YLENGTH, ZLENGTH);

	fgets(buffer, sizeof(buffer), f3);
	fgets(buffer, sizeof(buffer), f3);

	  for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH; z++) for(x=0; x<XLENGTH; x++) //outlet_geo
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		*(bnode + i_in1)=1;
		*(b_y + i_in1)=0;
		*(b_z + i_in1)=0;
		*(b_b + i_in1)=0;

	}


  for (i=0; i<=i_main3_2; i++)
  {
	   fscanf( f3, "%d\t%d\t%d\t%d\n", &x1, &y1, &z1, &*(bnode2 + i));
	   {
			

		   x2=XLENGTH-z1-1-10;
		   y2=y1;
		   z2=x1;
		   
		   if(x2>=0) if(y2>=0) if(z2>=0)
		   {
			
		   i_in2 = z2+ZLENGTH*(y2+YLENGTH*x2)+1;
		   *(bnode + i_in2)=1;
		   if(*(bnode2 + i)==1)*(bnode + i_in2)=0;

		   fprintf( f4, "%d\t%d\t%d\t%d\n", x2, y2, z2, *(bnode + i_in2));
		   }
	   }
  }
  
 
  	for(x=0; x<XLENGTH; x++) for(z=0; z<ZLENGTH; z++)   
	{

		mod1=0; mod2=0;
		for(y=YLENGTH-1;y>=0; y--)
		{
			i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
			if(mod1==0) if(*(bnode + i_in1)==0) mod1 = 1;
			if(mod1==1) if(*(bnode + i_in1)==1)
			{
				mod1=2;
				mod2 +=1;
			}
			if(mod1==2)
			{
				if(*(bnode + i_in1)==1) *(b_y + i_in1) = mod2;
				else mod1=1;
			}
		}

		for(y=YLENGTH-1;y>=0; y--)
		{
			i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
			if(*(b_y + i_in1) == mod2) *(b_y + i_in1) =0;
		}


	}


  	for(x=0; x<XLENGTH; x++) for(y=0; y<YLENGTH; y++)   
	{

		mod1=0; mod2=0;
		for(z=0;z<ZLENGTH; z++)
		{
			i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
			if(mod1==0) if(*(bnode + i_in1)==0) mod1 = 1;
			if(mod1==1) if(*(bnode + i_in1)==1)
			{
				mod1=2;
				mod2 +=1;
			}
			if(mod1==2)
			{
				if(*(bnode + i_in1)==1) *(b_z + i_in1) = mod2;
				else mod1=1;
			}
		}

		for(z=0;z<ZLENGTH; z++)
		{
			i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
			if(*(b_z + i_in1) == mod2) *(b_z + i_in1) =0;
		}

	}	

	
	//�޲ٱ�
	 for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH; z++) for(x=0; x<XLENGTH; x++) //outlet_geo
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(b_y + i_in1)==1) if(*(b_z + i_in1)==1) *(bnode + i_in1)=2;
	}

	for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH; z++) for(x=14; x<=15; x++)
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(b_y + i_in1)==1) *(bnode + i_in1)=2;
	}

	for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH; z++) for(x=16; x<=16; x++)
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(b_y + i_in1)==2) if(*(b_z + i_in1)==1) *(bnode + i_in1)=2;
	}

	for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH; z++) for(x=23; x<=23; x++)
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(b_y + i_in1)==1) if(*(b_z + i_in1)==3) *(bnode + i_in1)=2;
	}

	for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH/2; z++) for(x=18; x<=26; x++)
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(b_y + i_in1)==1) if(*(b_z + i_in1)==2) *(bnode + i_in1)=2;
	}

	for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH; z++) for(x=21; x<=31; x++)
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(b_y + i_in1)==2) if(*(b_z + i_in1)==1) *(bnode + i_in1)=2;
	}

	for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH; z++) for(x=40; x<=44; x++)
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(b_y + i_in1)==1) if(*(b_z + i_in1)>1) *(bnode + i_in1)=2;
	}

	for(y=0;y<YLENGTH; y++) for(z=17; z<75; z++) for(x=44; x<=58; x++)
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(b_y + i_in1)==1) *(bnode + i_in1)=2;
	}

	 for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH; z++) for(x=0; x<XLENGTH; x++) //outlet_geo
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		if(*(bnode + i_in1)==2) *(bnode + i_in1)=0;
	}





	 for(y=13;y<=18; y++) for(z=101; z<=106; z++) for(x=50; x<=50; x++) //2
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 *(bnode + i_in1)=1;
	 }


	 for(y=40;y<=46; y++) for(z=105; z<ZLENGTH; z++) for(x=56; x<=56; x++) //3
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 *(bnode + i_in1)=1;
	 }

	 for(y=57;y<=57; y++) for(z=104; z<ZLENGTH; z++) for(x=64; x<=69; x++) //5
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 *(bnode + i_in1)=1;
	 }

	 for(y=32;y<=40; y++) for(z=64; z<=72; z++) for(x=75; x<=75; x++) //6
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 *(bnode + i_in1)=1;
	 }

	 for(y=62;y<=68; y++) for(z=67; z<=67; z++) for(x=85; x<=90; x++) //7
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 *(bnode + i_in1)=1;
	 }

	 for(y=34;y<=42; y++) for(z=88; z<=90; z++) for(x=85; x<=91; x++) //8
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 *(bnode + i_in1)=1;
	 }

	for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH; z++) for(x=0; x<XLENGTH; x++) //outlet_geo
	{
		i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		*(bnode3 + i_in1)=*(bnode + i_in1);
	}

	 //�ܸ� ó��
	 for(y=33;y<=58; y++) for(z=40; z<=65; z++)
	 {
		 mod1=0;
		 for(x=44; x<=58; x++)
		 {
			 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
			 if(*(bnode + i_in1)==0) mod1=1;
			 if(mod1==1) *(bnode + i_in1)=0;
		 }
	 }



	 //Outlet ���� --> inlet으로 바꿔야함
	 //1 : x+ / 2: x- / 3: y- 4: z+
	 for(x=58; x<=58; x++) for(y=33;y<=58; y++) for(z=40; z<=65; z++)
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 if(*(bnode + i_in1)==0) *(b_b + i_in1)=7;//1
	 }

	 //for(x=33; x<=33; x++) for(y=33;y<=40; y++) for(z=99; z<=104; z++) //1 x+

	 /*
	 for(x=34; x<=34; x++) for(y=33;y<=40; y++) for(z=99; z<=104; z++) //1 x+
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 i_in2 = z+ZLENGTH*(y+YLENGTH*(x-1))+1;
		 *(bnode + i_in1) = *(bnode + i_in2);
		 if(*(bnode + i_in1)==0) *(b_b + i_in1)=1;//1
	 }
	 */

    // outlet    
	 for(x=50; x<=50; x++) for(y=13;y<=18; y++) for(z=101; z<=106; z++) //2 x+
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 i_in2 = z+ZLENGTH*(y+YLENGTH*(x-1))+1;
		 *(bnode + i_in1) = *(bnode + i_in2);
		 if(*(bnode + i_in1)==0) *(b_b + i_in1)=1;//1
	 }

    // outlet
	 for(x=56; x<=56; x++) for(y=40;y<=46; y++) for(z=105; z<ZLENGTH; z++) //3 x+
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 i_in2 = z+ZLENGTH*(y+YLENGTH*(x-1))+1;
		 *(bnode + i_in1) = *(bnode + i_in2);
		 if(*(bnode + i_in1)==0) *(b_b + i_in1)=1;//1
	 }

    // outlet
	 for(x=60; x<=60; x++) for(y=0;y<=4; y++) for(z=29; z<=33; z++) //4 x-
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 i_in2 = z+ZLENGTH*(y+YLENGTH*(x+1))+1;
		 *(bnode + i_in1) = *(bnode + i_in2);
		 if(*(bnode + i_in1)==0) *(b_b + i_in1)=2;
	 }

    // outlet
	 for(x=64; x<=69; x++) for(y=57;y<=57; y++) for(z=104; z<ZLENGTH; z++) //5 y-
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 i_in2 = z+ZLENGTH*((y+1)+YLENGTH*x)+1;
		 *(bnode + i_in1) = *(bnode + i_in2);
		 if(*(bnode + i_in1)==0) *(b_b + i_in1)=3;
	 }

    // outlet
	 for(x=75; x<=75; x++) for(y=32;y<=40; y++) for(z=64; z<=72; z++) //6 x-
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 i_in2 = z+ZLENGTH*(y+YLENGTH*(x+1))+1;
		 *(bnode + i_in1) = *(bnode + i_in2);
		 if(*(bnode + i_in1)==0) *(b_b + i_in1)=2;
	 }

    // outlet
	 for(x=84; x<=91; x++) for(y=61;y<=69; y++) for(z=67; z<=67; z++) //7 z+
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 i_in2 = (z-1)+ZLENGTH*(y+YLENGTH*x)+1;
		 *(bnode + i_in1) = *(bnode + i_in2);
		 if(*(bnode + i_in1)==0) *(b_b + i_in1)=4;
	 }

    // outlet
	 for(x=85; x<=93; x++) for(y=22;y<=22; y++) for(z=77; z<=86; z++) //8 y+
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 i_in2 = z+ZLENGTH*((y-1)+YLENGTH*x)+1;
		 *(bnode + i_in1) = *(bnode + i_in2);
		 if(*(bnode + i_in1)==0) *(b_b + i_in1)=5;
	 }

    // outlet
	  for(x=85; x<=91; x++) for(y=34;y<=42; y++) for(z=88; z<=88; z++) //8 z+
	 {
		 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
		 i_in2 = (z-1)+ZLENGTH*(y+YLENGTH*x)+1;
		 *(bnode + i_in1) = *(bnode + i_in2);
		 if(*(bnode + i_in1)==0) *(b_b + i_in1)=4;
	 }

    
	A2=0; 
	for(x=0; x<XLENGTH; x++)
	{
		A1=0;

		for(y=0;y<YLENGTH; y++) for(z=0; z<ZLENGTH; z++)
		{
			 i_in1 = z+ZLENGTH*(y+YLENGTH*x)+1;
			 if(*(bnode + i_in1)==0) A1++; // x수직 단면적
		}
		*(RGradient + x + 1) = A1-A2;
		A2=A1; // 이전 x값에서의 단면적
	}


	  fflush(f3);
	  fclose(f3);
	  fflush(f4);
	  fclose(f4);




}

void	pulsatile	(void)
{
     pulse[0] =0.100904;
     pulse[1] =0.100904;
     pulse[2] =0.096386;
     pulse[3]=0.0512050;
     pulse[4] =0.0165660;
     pulse[5] =0.0662650;
     pulse[6] =0.632530;
     pulse[7] =1.0;
     pulse[8] =0.573795;
     pulse[9] =0.394578;
     pulse[10] =0.144578;
     pulse[11] =-0.03012;
     pulse[12] =0.085843;
     pulse[13] =0.114458;
     pulse[14] =0.138554;
     pulse[15] =0.11747;
     pulse[16] =0.054217;
     pulse[17] =0.069277;
     pulse[18] =0.103916;
     pulse[19] =0.100904;

    return;
}