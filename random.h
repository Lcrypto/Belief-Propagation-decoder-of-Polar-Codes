#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#ifndef RANDOM_H
#define RANDOM_H

#define IM1  2147483563
#define IM2  2147483399
#define AM   (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.e-20
#define RNMX  (1.0-EPS)

// a uniform random number generator between zero and 1.
class Random
{
	long idum2;
   long idum;
	long iy;
	long iv[NTAB];

   unsigned memory;
   void init(long seed);
   double ran2();
public:
	Random(long seed);
   Random();
	double doublerandom();
   long longrandom(long range);
   bool boolrandom();
};

#endif  // RANDOM_H
