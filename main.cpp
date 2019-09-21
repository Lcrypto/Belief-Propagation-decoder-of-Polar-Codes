/*
Copyright(c) 2012, Vasiliy Usatyuk
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :
*Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and / or other materials provided with the distribution.
* Neither the name of the <organization> nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <stdio.h>
#include <math.h>
#include <memory.h>
//#include <new.h>
#include "random.h"


//#define _PRINT_TO_FILE  
#define  _PRINT_TO_SCREEN

#define   m				15 // log 2 from N (blocklenght)
#define   N				32768 // codes lenght 

//#define   m				11 // log 2 from N (blocklenght)
//#define   N				2048 // codes lenght 
#define   R_RATE		1./2.  // code rate
#define	  K			    (int)(R_RATE * N) // 




#define MAX(X,Y) (X-Y)>0. ?  X : Y
#define MIN(X,Y) (X-Y)<0. ?  X : Y
#define SIGN(X)  X>0. ?  1. : -1.


//My
#define   EB_NO_0			1.6
#define   delta_EB_NO		0.2
#define   number__EB_NO		30
#define	  iter4astep		100
//int BP_N=N;
const int max_itter=8;
double  P[m+1][N]; // LLR sets
int     B[m+1][N]; // u_i
int    ByPass[m+1][N];


// Array for BP decoder 
//double  P_BP[N][m+1]; 
//double L[max_itter][N][m+1];
//double R[max_itter][N][m+1];
//double scale=1-1>>4; //0.9375

// Array for BP decoder  scratch
//double  P_BP[N][m+1]; 
//double L_0[1][N]; //  initial LLR  to propagation from left to right
//double B_n[N][1];// initial R to propagation from right to left 
//double B_prev[N][N];
//double B_curr[N][N];
//double L_prev[N][N];
//double L_curr[N][N];
//double R[max_itter][N][m+1];
//double scale=1-1>>4; //0.9375

// We have i=N numbers rows and j=log2(N) columns.
// Array for BP decoder third (i think stable)
	int   u[N] ;
	int    out[N]; 
int   u_PB[N] ; //message for BP
double  P_BP[N][m+1]; 
//double R[max_itter][N][m+1];  // version 1 variable
//double scale=1-1>>4; //0.9375
//          lambda   phi          omega
//double L_BP[m+1][pow2(m+1)][pow2(m-pow2(m+1))]; //  initial LLR  to propagation from left to right
//double B_BP[m+1][pow2(m+1)][pow2(m-pow2(m+1))];
//start index from 1, this why +1 to dimmension
// to easy debug use representation with itteration index
//     itteration number|layer| index of nodes in layer
double L_BP[m+1+1][N+1]; //  initial LLR  to propagation from left to right
double B_BP[m+1+1][N+1]; //Insert LLR of message to belief propagation decoder 
//     

double cd_wr_md[N]; // code word after demodulation to check parity-check matrix

	int phi = 0;
	int beta = 0;
	int lambda = 0;
	// Initialization

	int  frozen_bits[N];

//for version 3 
//double R_O_BP[m+1][N][N];// initial R to propagation from right to left odd number phi
//double R_E_BP[m+1][N][N];// initial R to propagation from right to left even number phi


Random rnd;


int pow2(int x){
	return 1<<x;
}

int log2(int val)
{
	int res = 0;
	while(val>0)
	{
		val = val>>1;
		res++;
	}
	return res-1;
}

double jac_xor_bp( double x,  double y)
{
	return	SIGN(x)*SIGN(y)*MIN(fabs(x), fabs(y));
}
double jac_xor_approx_bp( double a,  double b)
{
	double minimum = MIN(fabs(a), fabs(b));
	if((a >= 0.0 && b < 0.0) || (a < 0.0 && b >= 0.0))
	{
		minimum *= -1.0;
	}
	double difference= fabs(a - b);
	double sum= fabs(a + b);
	//double log_probability_ratio = minimum + log(1.0+exp(-sum)) - log(1.0+exp(-difference));
	double log_probability_ratio = minimum + log((1.0+exp(-sum))/(1.0+exp(-difference)));
	//if(log_probability_ratio==0.)
	//{
	// printf("\n stop\n");
	//}
	return log_probability_ratio;	
}

// Blocklenght of Polar codes is N, number of nodes is : (log2(N)+1) columns (i-layers), N rows (j-index of nodes within that level 
// as examples n=1, ll give to us G^2,   
//  R          -->  |+|  <-- R(i+1,2-1)
//  L(i,j)     <--  |+|  --> L   
//                   | 
//  R(i,j+N/2) -->   =   <-- R (i+1,2j)
//  L          <--   =    -->L
//  
  //number of columns log2(N)+1
//                   i   j

//Variant on Arikan's tanner
// to correct work use i,j index of left upper corner of z-subgraph
/*void LLR_upd(long i, long j){   
	L_BP[i][j] = jac_xor_approx_bp(L_BP[i+1][j],L_BP[i+1][j+pow2(m-i)]+B_BP[i][j+pow2(m-i)]); 
	L_BP[i][j+pow2(m-i)] = jac_xor_approx_bp(B_BP[i][j], L_BP[i+1][j]) + L_BP[i+1][j+pow2(m-i)] ; //if row and columns start from 1
	
	printf("\n LLR update L_BP[%d][%d]  = xor(L_BP[%d][%d],L_BP[%d][%d]+B_BP[%d][%d]  = %f    ", i, j, i+1, j,i+1,j+pow2(m-i),i,j+pow2(m-i), L_BP[i][j] );
	printf("\n LLR update L_BP[i][j+N/2]  = %f  L_BP[i+1][2*j] %f  + jac_xor_approx_bp(B_BP[i][j] %f , L_BP[i+1][2*j-1] %f )  ", L_BP[i][j+N/2] , L_BP[i+1][2*j],B_BP[i][j], L_BP[i+1][2*j-1] );
	
	
	getchar();
}
void R_upd(long i, long j){   
B_BP[i+1][j] = jac_xor_approx_bp(B_BP[i][j], L_BP[i+1][j+pow2(m-i)]+B_BP[i][j+pow2(m-i)]) ; 
	B_BP[i+1][j+pow2(m-i)]= jac_xor_approx_bp(B_BP[i][j], L_BP[i+1][j] ) + B_BP[i][j+pow2(m-i)] ;
		printf("\n R update B_BP[%d][%d]  =  jac_xor_approx_bp(B_BP[%d][%d], L_BP[%d][%d]+B_BP[%d][%d] = %f", i+1,j, i, j, i+1, j+pow2(m-i), i, j+pow2(m-i), B_BP[i+1][j] );
		printf("\n R update B_BP[%d][%d]  =   jac_xor_approx_bp(B_BP[%d][%d], L_BP[%d][%d] ) + B_BP[%d][%d] = %f", i+1,j+pow2(m-i), i, j, i+1, j, i, j+pow2(m-i), B_BP[i+1][j] );
		getchar();
}
*/

//M1 representation of Arikan tanner's based on Pamuk's An FPGA Implementation Architecture for Decoding of Polar Codes
// it has  level upper interconnection of z-shape subgraph connected to first N/2 index ON EVERY Layers 
void LLR_upd(long i, long j){   
	L_BP[i][j] = jac_xor_approx_bp(L_BP[i+1][2*j-1],L_BP[i+1][2*j]+B_BP[i][j+N/2]); 
	L_BP[i][j+N/2] = jac_xor_approx_bp(B_BP[i][j], L_BP[i+1][2*j-1]) + L_BP[i+1][2*j] ; //if row and columns start from 1
	
	//printf("\n LLR update L_BP[%d][%d]  = xor(L_BP[%d][%d],L_BP[%d][%d]+B_BP[%d][%d]  = %f ", i, j, i+1,  2*j-1,i+1,2*j,i,j+N/2, L_BP[i][j] );
	//printf("\n LLR update L_BP[%d][%d]  = xor(B_BP[%d]][%d]], L_BP[%d]][%d]]) + L_BP[%d]][%d]]   = %f ", i,j+N/2, i, j ,i+1,2*j-1,i+1,2*j ,L_BP[i][j+N/2]  );
	
	
	//getchar();
}
void R_upd(long i, long j){   
	B_BP[i+1][2*j-1] = jac_xor_approx_bp(B_BP[i][j], L_BP[i+1][2*j]+B_BP[i][j+N/2]) ; 
	B_BP[i+1][2*j] = jac_xor_approx_bp(B_BP[i][j], L_BP[i+1][2*j-1] ) + B_BP[i][j+N/2] ;
	//printf("\n R update B_BP[%d][%d]  =  jac_xor_approx_bp(B_BP[%d][%d], L_BP[%d][%d]+B_BP[%d][%d] = %f", i+1,2*j-1, i, j, i+1, 2*j, i,j+N/2, B_BP[i+1][2*j-1]  );
	//printf("\n R update B_BP[%d][%d]  =   jac_xor_approx_bp(B_BP[%d][%d], L_BP[%d][%d] ) + B_BP[%d][%d] = %f", i+1,2*j, i, j, i+1, 2*j-1, i, j+N/2, B_BP[i+1][2*j] );

	//getchar();
}






//Scrach version 2 

// Kernel of BP
// count going from 0, shift need to decsease by 1 from i+N/(2^j)
// inside after size of connection divide by 2 
/*void belief_propagation_kernel(int curr_itter, long i, long j, double scale){   
	//L[curr_itter][i][j] = scale*jac_xor_bp(L[curr_itter-1][i][j+1],L[curr_itter-1][i+N/pow2(j+1)][j+1]+R[curr_itter][i+N/pow2(j+1)][j]); 
	//L[curr_itter][i+N/(1<<j)][j] = L[curr_itter-1][i+N/(pow2(j+1))][j+1] + scale* jac_xor_bp(L[curr_itter-1][i][j+1], R[curr_itter][i][j]);
	//R[curr_itter][i][j+1] = scale* jac_xor_bp(R[curr_itter][i][j], L[curr_itter-1][i+N/(pow2(j+1))][j+1]+R[curr_itter][i+N/(pow2(j+1))][j]) ;  //if row and columns start from 1
	//R[curr_itter][i+N/(pow2(j+1))][j]=R[curr_itter][i+N/(pow2(j+1))][j] + scale* jac_xor_bp(L[curr_itter-1][i][j+1], R[curr_itter][i][j]);
	printf("\n Current itteration %d, (%d,%d) ", curr_itter, i,j);
	L[curr_itter][i][j] = scale*jac_xor_bp(L[curr_itter-1][i][j+1],L[curr_itter-1][i+N/pow2(j+1)-1][j+1]+R[curr_itter][i+N/pow2(j+1)-1][j]); 
	L[curr_itter][i+N/(pow2(j+1)-1)][j] = L[curr_itter-1][i+N/(pow2(j+1))-1][j+1] + scale* jac_xor_bp(L[curr_itter-1][i][j+1], R[curr_itter][i][j]);  
	R[curr_itter][i][j+1] = scale* jac_xor_bp(R[curr_itter][i][j], L[curr_itter-1][i+N/(pow2(j+1))-1][j+1]+R[curr_itter][i+N/(pow2(j+1))-1][j]) ;
	R[curr_itter][i+N/(pow2(j+1)-1)][j+1]=R[curr_itter][i+N/(pow2(j+1))-1][j] + scale* jac_xor_bp(L[curr_itter-1][i][j+1], R[curr_itter][i][j]);
}
*/


// version 3 space efficient 

/*
void updateLLR(int lambda,int phi){
	if (lambda==0) return;
	 int psi = phi>>1;
	if(phi%2==0) updateLLR(lambda-1, psi);
	for (int w= 0; w<=(pow2(m-lambda)-1);w++){
		
		if(phi%2==0) L_BP[lambda][phi][w]=jac_xor_bp(L_BP[lambda-1][psi][2*w], L_BP[lambda-1][psi][2*w+1]+R_E_BP[lambda][phi+1][w]);
		else L_BP[lambda][phi][w]=L_BP[lambda-1][psi][2*w+1] + jac_xor_bp(L_BP[lambda-1][psi][2*w], R_O_BP[lambda][phi-1][w]);
	}
}
void updateR(int lambda, int phi){
	int psi = phi>>1;
	if(phi %2 != 0){
		for( int w=0; w<=(pow2(m-lambda)-1);w++){
			if (psi%2 ==0){
				R_E_BP[lambda-1][psi][2*w]= jac_xor_bp(R_E_BP[lambda][phi-1][w],R_O_BP[lambda][phi][w]+L_BP[lambda-1][psi][2*w+1]);
				R_E_BP[lambda-1][psi][2*w+1] = R_O_BP[lambda][phi][w]+jac_xor_bp(R_E_BP[lambda][phi-1][w], L_BP[lambda-1][psi][2*w]);
							}
			else{
			R_O_BP[lambda-1][psi][2*w]=jac_xor_bp(R_E_BP[lambda][phi-1][w], R_O_BP[lambda][phi][w]+L_BP[lambda-1][psi][2*w+1]);
			R_O_BP[lambda-1][phi][2*w+1] = R_O_BP[lambda][phi][w] + jac_xor_bp(R_E_BP[lambda][phi-1][w], L_BP[lambda-1][psi][2*w]);
				}
	
		
							}
		if (psi %2 !=0) updateR(lambda-1,psi);
												}

}

*/

//version 2 stable SC sheduler in BP decoder
/*
void updateR(int lambda,int phi){
	//printf("\n R propagate \n");
	int psi = phi>>1;
	if ( phi %2 ==1 ){
		for (int i=0; i<=(pow2(m-lambda)-1);i++){
			
			B_BP[lambda-1][psi][2*i]= jac_xor_approx_bp(B_BP[lambda][phi-1][i],B_BP[lambda][phi][i]+L_BP[lambda-1][psi][2*i+1]);
			B_BP[lambda-1][psi][2*i+1]=B_BP[lambda][phi][i]+  jac_xor_approx_bp(B_BP[lambda][phi-1][i],L_BP[lambda-1][psi][2*i]); // some problem here
			//if(B_BP[lambda-1][psi][2*i]!=1.E+300&&	B_BP[lambda-1][psi][2*i]!=-1.E+300) if ( abs(B_BP[lambda-1][psi][2*i] ) > 30000) B_BP[lambda-1][psi][2*i]= si30000		
			 //printf("\n B_BP[lambda-1][psi][2*omega] is %f.  \n (A) XOR (B+C) \n   %f XOR (%f + %f.) ", B_BP[lambda-1][psi][2*i], B_BP[lambda][phi-1][i], B_BP[lambda][phi][i], L_BP[lambda-1][psi][2*i+1]);
			// printf("\n B_BP[lambda-1][psi][2*omega+1] is %f. \n (A XOR B)+C \n (%f.XOR %f.) + %f. ", B_BP[lambda-1][psi][2*i+1], B_BP[lambda][phi-1][i], L_BP[lambda-1][psi][2*i] , B_BP[lambda][phi][i]);
		//  printf("\n Psi is %d. Lambda is %d. Phi is %d. Omega is %d \n", psi, lambda, phi,i);
	// getchar();
		}
	
	if(psi%2 == 1) updateR(lambda-1,psi);
	}
}

void updateLLR( int lambda,int  phi){
	// printf("\n LLR propagate \n");
	//if (lambda==0)  printf("Return from updateLLR. Lambda is 0. Phi is %d\n", lambda, phi);
	if (lambda==0) return;
	 int psi = phi>>1;
	// printf("Psi is %d. Lambda is %d. Phi is %d.", psi, lambda, phi);
	// getchar();
	 if (phi% 2 ==0) updateLLR(lambda-1,psi );  
	
	for(int i=0; i<=(pow2(m-lambda)-1);i++){
	//printf("\n 2^(log2 N - lambda)-1 = %d \n ", pow2(m-lambda)-1);
		if(phi%2==0) {
			L_BP[lambda][phi][i]= jac_xor_approx_bp(L_BP[lambda-1] [psi][2*i], L_BP[lambda-1][psi][2*i+1]+B_BP[lambda][phi+1][i]);
		   
		
	//	printf("\n L_BP[lambda][phi][omega] is %f. \n (A) XOR (B+C) \n (%f.)XOR(%f. + %f.) \n", L_BP[lambda][phi][i], L_BP[lambda-1] [psi][2*i], L_BP[lambda-1][psi][2*i+1] , B_BP[lambda][phi+1][i]); 
	// printf("Psi is %d. Lambda is %d. Phi is %d. Omega is %d \n", psi, lambda, phi, i);

	// getchar();
		}
		else 
			{
				L_BP[lambda][phi][i] = L_BP[lambda-1] [psi] [2*i+1]+ jac_xor_approx_bp(L_BP[lambda-1][psi][2*i],B_BP[lambda][phi-1][i]); 
	//printf("Result is %f. \n (A XOR B) + (C) : \n (%f XOR %f) + %f. \n", L_BP[lambda][phi][i] , L_BP[lambda-1][psi][2*i],B_BP[lambda][phi-1][i], L_BP[lambda-1] [psi] [2*i+1]);
		
	//	getchar();
			}

	}
	
}


*/
double rand_gauss (double No)
{
	double v1,v2,s; 
	do
	{ 
		v1 = 2. *  rnd.doublerandom() - 1.;  
		v2 = 2. *  rnd.doublerandom() - 1.;  
		s = v1*v1 + v2*v2;  
	} while ( s >= 1. );
	if (s == 0.)     return 0.;  
	else  
		return (sqrt(No/2.)*v1*sqrt(-2. * log(s) / s)); 
} 

double jac_xor_approx(const double x, const double y)
{
	return	SIGN(x)*SIGN(y)*MIN(fabs(x), fabs(y));
}




/***************************************************************************************/

static inline double jac_xor(const double a, const double b)
{
	double minimum = MIN(fabs(a), fabs(b));
	if((a >= 0.0 && b < 0.0) || (a < 0.0 && b >= 0.0))
	{
		minimum *= -1.0;
	}
	double difference= fabs(a - b);
	double sum= fabs(a + b);
	//double log_probability_ratio = minimum + log(1.0+exp(-sum)) - log(1.0+exp(-difference));
	double log_probability_ratio = minimum + log((1.0+exp(-sum))/(1.0+exp(-difference)));
	//if(log_probability_ratio==0.)
	//{
	// printf("\n stop\n");
	//}
	return log_probability_ratio;	
}


//static inline 
//double jac_xor(double x, double y)
//{
//    return (x + y > 0. ? x + y : 0.) + log(1. + exp(-fabs(x + y)))
//         - (x > y ? x : y) - log(1. + exp(-fabs(x - y)));
//}

/***************************************************************************************/

//double tanh(double x)
//{
//	double res;	
//	res = (exp(x)-exp(-x))/(exp(x)+exp(-x));
//	return res;
//} 

/***************************************************************************************/

double atanh(double x)
{
	double res;	
	//res = .5 * (log(1.+x)-log(1.-x)); // |x|< 1.
	res = .5 * log((1.+x)/(1.-x)); // |x|< 1.
	return res;
} 
/***************************************************************************************/

int reversBit(unsigned int val) // Arikan reverse shuffle 
{
	unsigned int rev_val = 0;
	for(int i=0;i<m;i++)
	{
		if(val&(1<<i))
			rev_val |= 1<<((m-1)-i);
	}
	return  rev_val;
}

/*********************************************/
int ptr(int lambda, int phi, int beta)   // return 2^lambda*beta+phi 
{
	int res	= phi + (1<<lambda)*beta;
	return res;
}

/**********************************************/
void recursivelyUpdateB (int lambda, int phi)
{  
	int psi = phi>>1;     // phi/2
	for (int beta=0;beta<=(1<<(m-lambda))-1;beta++) // before beta<=2^(m-lambda)-1
	{
		B[lambda-1][ptr(lambda-1, psi,2*beta)] = B[lambda][ptr(lambda,phi-1,beta)] ^ B[lambda][ptr(lambda,phi,beta)]; // B[lambda-1][2^(lambda-1)*2*beta+phi/2]  = B[lambda] [2^(lambda)*beta+ phi-1]^ B[lambda] [2^lambda*beta+phi] ;
		B[lambda-1][ptr(lambda-1, psi,2*beta+1)] = B[lambda][ptr(lambda,phi,beta)];                                   // B[lambda-1][2^(lambda-1)*2*(beta+1)+phi/2] =  B[lambda] [2^lambda*beta+phi] ;
	}
	if ((psi%2)==1)  
	{
		recursivelyUpdateB (lambda-1, psi);
	}
}


 
// This function decode using 
// 
//
//
/**********************************************/
void recursivelyCalcP (int lambda, int phi)        
{  
	if (lambda==0) return; // Stopping condition 
	int psi = phi>>1;       
	if ((phi%2)==0) recursivelyCalcP (lambda-1, psi); // if phi even resursive calculate P with (lambda-1, phi/2);

	for (int beta=0;beta<=(1<<(m-lambda))-1;beta++)
	{
		if ((phi%2)==0) // phi even 
		{
			//double  x1 =  tanh(.5 * P[lambda-1][ptr(lambda-1, psi,2*beta)]);
			//double  x2 =  tanh(.5 * P[lambda-1][ptr(lambda-1, psi,2*beta+1)]);
			//P[lambda][ptr(lambda, phi,beta)] = 2. * atanh(x1 * x2);
			double  x1 =  P[lambda-1][ptr(lambda-1, psi,2*beta)];      //   x1 =  P[lambda-1][2^(lambda-1)*(2*beta)+phi/2];
			double  x2 =  P[lambda-1][ptr(lambda-1, psi,2*beta+1)];    //   x2 =  P[lambda-1][2^(lambda-1)*(2*beta+1)+phi/2];
			P[lambda][ptr(lambda, phi,beta)] = jac_xor(x1,x2);         //   P[lambda][ 2^lambda*beta+phi] =                      calculate f function 
			// double  x1 =   P[lambda-1][ptr(lambda-1, psi,2*beta)];
			// double  x2 =   P[lambda-1][ptr(lambda-1, psi,2*beta+1)];
			// P[lambda][ptr(lambda, phi,beta)] = (1. + x1 * x2)/(x1+x2);

		}
		else         // phi uneven 
		{
			int u = B[lambda][ptr(lambda, phi-1,beta)];

			double x1 = P[lambda-1][ptr(lambda-1, psi,2*beta)];    // LLR
			double x2 = P[lambda-1][ptr(lambda-1, psi,2*beta+1)];  // LLR
			P[lambda][ptr(lambda, phi,beta)] = x1*pow(-1.,u) +  x2; //g function  

			//double  x1 =  pow(P[lambda-1][ptr(lambda-1, psi,2*beta)],1-2*u);
			//double  x2 =  P[lambda-1][ptr(lambda-1, psi,2*beta+1)];
			//	P[lambda][ptr(lambda, phi,beta)] = x1*x2;

		}
	}
}
/*******************************************************************/


/*******************************************************************/
void encode (int *u,int *x, int n)
{
	int tu[N/2];
	if(log2(n)==1)
	{
		x[0] = u[0]^u[1];
		x[1] = u[1];
		return;
	}
	else
	{
		for(int i=0;i<=n/2-1;i++)
			tu[i] = u[i]^u[i+n/2];

		int x1[N];
		encode(tu, x1,n/2);

		for(int i=0;i<=n/2-1;i++)
			tu[i] = u[i+n/2];

		int x2[N];
		encode(tu, x2,n/2);

		for(int i=0;i<=n/2-1;i++)
		{
			x[i] = x1[i];
			x[i+n/2] = x2[i];
		}
	}
	return;
}

/*******************************************************************/

void main(void)
{


	int	x[N] ;
	double  s[N]; 

	int    ptr_fr_bit;
	double tt; 
	double dB;
	double No;
	double snr,sigma;
	int rev_bit_fr[N];
	char fr_char;
	double	   rate;	
	
	
	//dB = 10*log(Eb/No) where Eb is 1
	//FILE * fdmp = fopen("log.txt","w");
	FILE * fdmp = fopen("Frozen bit.txt","w");
	//FILE * fid = fopen("fr.txt","r"); // used for ARIKAN algorithms
   //FILE * fid1 = fopen("z_2048_05_21.txt","r");
	//FILE * fid1 = fopen("z8_r_0.5.txt","r");   // 8 bit, (4 frozen and source) toy examples for debug
	//FILE * fid1 = fopen("z_256_0.5_3DB.txt","r");   // 256 bit, (42 frozen and source) toy examples for debug
	//FILE * fid1 = fopen("z_256_0.5_4DB.txt","r");
	//FILE * fid1 = fopen("BL512_0.5_2.2 db.txt","r");
	//FILE * fid1 = fopen("8_7_over_8.txt","r");
	//FILE * fid1 = fopen("1024_rate_0_5_2.2.txt","r");



	FILE * fid1 = fopen("32768_1_2.txt","r");
	//FILE * fid1 = fopen("z2048_R0.8_Db2.2.txt","r");

	//FILE * fid1 = fopen("BL_512_rate_0.8_2.2_DB.txt","r");
	//FILE * fid1 = fopen("BL_512_0.5_2.2.txt","r");
	
	
	//FILE * fid1 = fopen("BL_256_0.5_2.2.txt","r");
	//FILE * fid1 = fopen("BL_256_0.8_2.2.txt","r");

	//FILE * fid1 = fopen("BL_128_0.8_2.2.txt","r");
	//FILE * fid1 = fopen("BL_128_0.5.2.2.txt","r");


	//FILE * fid1 = fopen("BL_4_r_1.2_2.2DB.txt","r");



	//FILE * fid1 = fopen("fr_8_3_8.txt","r");
	//FILE * fid1 = fopen("fr_16_4_16.txt","r");
	//FILE * fid1 = fopen("BL_2048_0.8_2.2.txt","r");
	//FILE * fid1 = fopen("fr_8_3_8.txt","r");
	//FILE * fid1 = fopen("BL_2048_rate0_5_2_0DB.txt","r");   // USE INVERSE ORDER FROZEN BIT MASK!!!
	//FILE * fid1 = fopen("BL_4_r_1.4_2.2DB.txt","r");
	
	//FILE * fid1 = fopen("z256_r0.5_2db.txt","r");
	//FILE * fid1 = fopen("4096_0.8_3db.txt","r");
	//FILE * fid1 = fopen("BL_2048_R0.8_3db.txt","r");
	
	//FILE * fid1 = fopen("BL2048_0.8_4.0 db.txt","r")
	//FILE * fid1 = fopen("z2048_R0.8_Db2.2.txt","r");
	//FILE * fid1 = fopen("z16_r_0.75.txt","r");
	//FILE * fid1 = fopen("z_32768_093_7.txt","r");
	//FILE * fid1 = fopen("z256_r0.8_3b.txt","r");
	FILE * cdwr_aft_modul = fopen("code_words_after_encode.txt","w");
	//FILE * sim_res = fopen("result_BL_2048_rate0_5_2_2DB.txt","w");
	FILE * sim_res = fopen("result_BL_8_rate3_8_2_2DB.txt","w");
	FILE * perm = fopen("Permutation matrix.txt","w");
	//FILE * inv_froz = fopen("Frozen bit_inverse.txt","w");
	//FILE * frozen_bit = fopen("Frozen bit.txt","w");

#ifdef _PRINT_TO_SCREEN
{
	for (int i=0;i<=N-1;i++)
	{
				fscanf(fid1,"%d",&fr_char);	
		if(fr_char==1)  
			frozen_bits[i] = 1;        
		else 
			frozen_bits[i] = 0;
			
		//printf("%d", frozen_bits[i]);
	}


	//	int ptr = 0; 
	//	for (int i=0;i<=N-1;i++)
	//		if(info_bits[ptr]==(i+1))
	//		{
	//			frozen_bits[i] = 0;
	//			ptr++;
	//		}
	//		else
	//			frozen_bits[i] = 1;


		for (int i=0;i<=N-1;i++)
			rev_bit_fr[i] =  frozen_bits[reversBit(i)];

		

		for (int i=0;i<=N-1;i++)
			frozen_bits[i] = rev_bit_fr[i];
		
	}
		


	

	int info_bits1=0;
	for (int i=0;i<=N-1;i++)
		if(frozen_bits[i]==1) info_bits1++;
	printf("Frozen_bits=%d\n",info_bits1);
	printf("Blocklenght=%d\n",N);
	printf("BP iteration=%d\n",max_itter);
#else//original Arikan
	int info_bit = 0;
	for (int i=0;i<=N-1;i++)
	{
		fscanf(fid,"%d\n",(frozen_bits+i));
	  if(*(frozen_bits+i)==0)info_bit++;
	
	}
	printf("\n K = %d\n",info_bit);
#endif// _PRINT_TO_SCREEN


	for(int iii = 0; iii < number__EB_NO; iii++) {

			dB = EB_NO_0 + iii * delta_EB_NO;


			//No = 1./pow(10.,(dB-3.)/10.);      //convert Eb/N0 from unit db to normal numbers
		    
			rate = (double)K/(double)N;
			printf("rate=%f\n", rate);
			No = 1./pow(10.,(dB+10.*log10(rate))/10.);      //convert Eb/N0 from unit db to normal numbers
			printf("No=%f \n", No);
			
			int iter = 0;
			int err = 0;
			int fer = 0;
			int err_BP  = 0;
			int fer_BP = 0;
//			printf("\n");

			while(iter < iter4astep)   
			{
				iter++;
//				printf("inter=%d\r", iter);
				//for (int i=0;i<=N-1;i++) u[i] = 0;

				/*for (int i=0;i<=N-1;i++)
					for (int j=0;j<=m;j++)
						P[j][i] = 0;
						*/

		//Let's Rock belief propagation 8)
					
		for(int i=0;i<=m+1;i++)
		for(int j=0;j<=N+1;j++)
		L_BP[i][j]=0., B_BP[i][j]=0. ; //  initial LLR  to propagation from left to right
		
						


				
				for (int i=0;i<=N-1;i++)
				{
					if(frozen_bits[i]==1)
				 {
					 u[i] = 0;
				 }
					else
				 {
					 u[i] = rnd.boolrandom();
					//u[i] = 1;  // for debug
				 }
				}
				

				encode(u, x, N);

			/*	for (int i=0;i<=N-1;i++)  printf("Frozen bit [%d]: %d \n",i, frozen_bits[i]); 
					for (int i=0;i<=N-1;i++)   printf("Message before encode [%d]: %d \n",i, u[i]); 
					getchar();
					for (int i=0;i<=N-1;i++)   printf("Codeword after encode [%d]: %d \n",i, x[i]); 
					getchar();
					*/
				for (int i=0;i<=N-1;i++) cd_wr_md[i]=x[i]; //check x after encode
				
					
					
					//bpsk modulation
				for (int i=0;i<=N-1;i++)	(x[i]==1)? s[i]= -1. : s[i]= 1.; 
			
				
				//add AWGN

				for (int i=0;i<=N-1;i++)	s[i] += rand_gauss(No);


				//demodulation


				for (int i=0;i<=N-1;i++)
				{
					//P[0][i] = exp(2.*s[i]/No);
					P[0][i] = 4.*s[i]/No;                    // log from LR for AWGN 
					//printf("LLR from channel [%d]: %f", i , P[0][i]);
					//getchar();
				}
			 
				//For Belief propagation decoding LLR  
				//in debug ll use last one LLR
					for (int i=0;i<N;i++)
				{
					//P[0][i] = exp(2.*s[i]/No);
					P_BP[i][0] = P[0][i];
			
				} 
				for (int i=0;i<N;i++)
				{
					
					//L[0][i][m+1] =P[0][i];
					L_BP[m+1][i+1] =P_BP[i][0];
					//printf("\n LLR %d is: %f \n",i, L_BP[0][0][i]);
				}
		
				// R 
				for (int i=0;i<N;i++)
				{
					if(frozen_bits[i]==1){ //if  frozen and i odd
					//if(frozen_bits[reversBit(i)]==1){ //if  frozen and i odd   and interconnecting index function without Arikan's shuffle 
				B_BP[1][i+1] = 1.8E+30; 
					 //if(u[i]==0.)  B_BP[m][i][0] = 1.8E+30;    //toy benchmark if we maybe know input 
					  //if (u[i]==1.)  B_BP[m][i][0] = -1.8E+30;
					}
					else
						{
						B_BP[1][i+1] = 0.;
						}
					 //if (u[i]==0.)  R[0][i][0] = 1.8E+305;
					 //if (u[i]==1.0)    R[0][i][0] = -1.8E+305; 
					 
			
					
					//printf("\n R %d is: %f \n",i,  B_BP[m][i][0]);
				}
				/********************************************************/ 
				// DECODING LOOP
				/********************************************************/
		

				for(phi=0;phi<=N-1;phi++)
				{
					recursivelyCalcP(m, phi);

		#ifdef _PRINT_TO_FILE								
					fprintf(fdmp,"\n***********************************\n");
					for(int i=0;i<=phi;i++)
						fprintf(fdmp,"%3d %3d %3d %3.8f \n", i,frozen_bits[reversBit(i)],u[reversBit(i)], P[m][i]);
					fprintf(fdmp,"\n");

					fprintf(fdmp,"\n#########################################################\n");
					for(int i=0;i<=m;i++)
					{				
						for(int j=0;j<=(1<<m)-1;j++)
							fprintf(fdmp,"%3.8f ", P[i][j]);
						fprintf(fdmp,"\n");
					}

					fprintf(fdmp,"\n");	
					for(int i=0;i<=phi-1;i++)
					{
						fprintf(fdmp,"%d",  u[reversBit(i)]);
					}
					fprintf(fdmp,"\n\n");


					for(int i=0;i<=phi-1;i++)
					{
						fprintf(fdmp,"%d", B[m][i]);
					}



		#endif// _PRINT_TO_FILE				


					if(frozen_bits[reversBit(phi)]==1)
					{
						B[m][phi] = 0;
					}
					else
					{
						if(P[m][phi]>0.)   // estimated u_i bits (in B array) using standart rules 
						{
							B[m][phi] = 0;
						}
						else
						{
							B[m][phi] = 1;
						}
					}
					if((phi%2)==1) 
						recursivelyUpdateB(m, phi);
				}
				



				int terr = err;
				int terr_BP = err_BP;
				// number of itterations 
				
				
					
				/*	for(int phi=0; phi<=N-1;phi++){
						updateLLR(m,phi);
						if(phi % 2 !=0) updateR(m,phi); 
					}
					*/

					//BP decoder
				
			
					
				
				
				for(int i=1;i<=max_itter;i++){




						for(phi=m; phi>=1;phi--){

					for(int k=1; k<=N/2;k++){
						LLR_upd(phi,k);
					

					}
					}
				
					
						for(phi=1; phi<=m;phi++){

					for(int k=1; k<=N/2;k++){
						R_upd(phi,k);
						

					}
					}
						
				
					
				
			
				

					
					//printf("Where here");

					
			//		printf("\n Itteration number %d\n ", i);
			//	for(int i=0;i<N;i++) printf("LLR after belief propagate: %8.8f \n", L_BP[m][i][0]); //N
			//		for(int i=0;i<N;i++) printf("R after belief propagate: %8.8f \n", B_BP[m][i][0]); //N
			//	for(int i=0;i<N;i++) printf("R after belief propagate: %8.8f \n", B_BP[0][0][i]); //N
				
				//getchar();
				}
									
				
				float vv[N];
				float debug[N];
					 
						
				
						for(int i=1;i<=N;i++)	 if ((B_BP[m+1][i]+L_BP[m+1][i])>=0.)  u_PB[i-1]=0 ;											
					else 	u_PB[i-1]=1;
					 // Codeword after BP

					 for(int i=0;i<N;i++) debug[i]= u_PB[reversBit(i)];
					 // information symbols
					 for(int i=1;i<=N;i++)	 if ((B_BP[1][i]+L_BP[1][i])>=0.)  vv[i-1]=0 ;											
					else 	vv[i-1]=1;
				
				/*	for(int i=0;i<N;i++)	{
						printf("\n LLR is equal: %f", L_BP[m+1][i+1]);
						printf("\n LLR+R is equal: %f", B_BP[m+1][i+1]+L_BP[m+1][i+1]);
						printf("\n Decoded message %d: %d", i ,  u_PB[i]);
						printf("\n Sent  codeword %d: %d\n",i, u[i]);
						printf("\n Sent  codeword after shuffle %d: %f\n",i,debug[i]);
					}
					*/
				
					//for(int i=0;i<=N-1;i++) printf("\n B_BP[0][0][%d] propagation is equal: %f. \n  ",i, B_BP[0][0][i]); 
					//getchar();

					//for(int i=0;i<N;i++) printf("LLR after belief propagate: %8.8f \n", L_BP[m][i][0]); //N
					//for(int i=0;i<N;i++) printf("R after belief propagate: %8.8f \n", B_BP[m][i][0]); //N
	//for(int i=0;i<N;i++) printf("R after belief propagate: %8.8f \n", B_BP[m][i][0]); //N
	
					//for(int i=0;i<N;i++)  printf(" Input: %d , SC output: %d , BP output: %d \n", u[i], out[i], u_PB[reversBit(i)]); //N
					// BER, FER FOR BP
					
					//errors for codeword
					//for(int i=0;i<=N-1;i++) if(u_PB[i]!=x[i])err_BP++;
					//errors for information symbols
					for(int i=0;i<=N-1;i++) if(vv[i]!=u[i])err_BP++;

			
				if((terr_BP-err_BP)!=0) fer_BP++;
			
		


			//	if((terr_BP-err_BP)==0)	
			//		{
			//			for(int i=0;i<=N-1;i++)  fprintf(cdwr_aft_modul, "%f ",cd_wr_md[i]); 
			//			fprintf(cdwr_aft_modul, "\n ");
			//	}
					
				// BER, FER FOR SC
				for(int i=0;i<=N-1;i++)
				{
					out[i] = B[m][reversBit(i)];
					if(out[i]!=u[i])err++;
					
				}
				

				if((terr-err)!=0)
				{
					fer++;
				}
				//if((terr-err)==0)	
				//	{
						//for(int i=0;i<=N-1;i++)  fprintf(cdwr_aft_modul, "%f ",cd_wr_md[i]); 
						//fprintf(cdwr_aft_modul, "\n ");
			//	}
				if(iter%iter4astep==0)
				//if(iter%100==0)
				{
					

					printf("SC Eb/NO %5.3f  frame  %d ber %e  fer %d  %e \n", dB, iter, (double)err/(double)(iter*K), fer,(double)fer/(double)(iter));
//					fprintf(sim_res, "Eb/NO %5.3f  frame  %d ber %e  fer %d  %e \n", dB, iter, (double)err/(double)(iter*K), fer,(double)fer/(double)(iter));
					fprintf(sim_res, "%5.3f  %e  %e \n", dB, (double)err/(double)(iter*K), (double)fer/(double)(iter));
					fprintf(sim_res, "%BP 5.3f  %e  %e \n", dB, (double)err_BP/(double)(iter*K), (double)fer_BP/(double)(iter));
					printf("BP (%d) Eb/NO %5.3f  frame  %d ber %e  fer %d  %e \n", max_itter,  dB, iter, (double)err_BP/(double)(iter*K), fer_BP,(double)fer_BP/(double)(iter));
					//getchar();
				}
			}
	}
	
		

		


					//for(int i=0;i<N;i++) printf("\n Right propagation is equal: %f. \n Left propagation is equal: %f ",B_BP[m][i][0], L_BP[m][i][0]); 
					
						//for(int i=0;i<N;i++) printf("\n u[%d] from BP: %f \n",i, u_PB[i]); //N
			
					


					//for(int i=0;i<N;i++) printf("u[%d] from BP: %10.8f \n",i, u_PB[i]); //N
	//for(int i=0;i<N;i++) printf("LLR after belief propagate: %10.8f \n", L_BP[0][0][i]); //N
	//for(int i=0;i<N;i++) printf("R after belief propagate: %10.8f \n", B_BP[m][i][0]); //N


		//for(int i=0;i<N;i++) printf("LLR after belief propagate: %10.8f \n", L_BP[m][i][0]); //N
	//for(int i=0;i<N;i++) printf("R after belief propagate: %10.8f \n", B_BP[0][0][i]); //N
					

				
				
					



/* Version 3 decoder
				for(int i=1;i<=max_itter;i++){
					for(phi=0; phi<=N-1;phi++){
						updateLLR(m,phi);
						if (phi % 2 !=0) R_E_BP[m][phi][0]= 1.E+30; 
							else updateR(m,phi); 
					}
					printf("\nItteration number: %d \n",i);
						for(int i=0;i<N;i++) printf("u[%d] from BP: %10.8f \n",i, u_PB[i]); //N
	for(int i=0;i<N;i++) printf("LLR after belief propagate: %10.8f \n", L_BP[0][0][i]); //N
	//for(int i=0;i<N;i++) printf("R after belief propagate: %10.8f \n", R_BP[m][i][0]); //N
		for(int i=0;i<N;i++) printf("R after belief propagate: %10.8f \n", R_E_BP[0][0][i]); //N
					}
				
				for(int i=0;i<N;i++)
					if ((B_BP[m][i][0]+L_BP[m][i][0])>=0.)    u_PB[i]=0 ;
					else u_PB[i]=1;
*/
	
	

	
		
	fclose(sim_res);
	for (int i=0;i<=N-1;i++) 	fprintf(perm,"%d ", reversBit(i));
	fclose(perm);
	fclose(cdwr_aft_modul);
	getchar();
}


