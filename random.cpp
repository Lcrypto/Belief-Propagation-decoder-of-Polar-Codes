#include "random.h"

/*
Long period (? 2 \Theta 10 18 ) random number generator of L'Ecuyer with Bays­Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1.
--
*/




double Random::ran2()
{
	int j;
	long k;
	double temp;


	k=(idum)/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;  // Compute idum=(IA1*idum) % IM1 without overflows by Schrage's method.
	if (idum < 0)
		idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;	// Compute idum2=(IA2*idum) % IM2 likewise.
	if (idum2 < 0)
   	idum2 += IM2;
   j = iy/NDIV;
   iy=iv[j]-idum2;
//	iy=iv[j]-idum2; 		// Here idum is shuffled, idum and idum2 are combined to generate output.
	iv[j] = idum;
	if (iy < 1)
		iy += IMM1;
	if ((temp=AM*iy) > RNMX)
		return RNMX; 		// Because users don't expect endpoint values.
	else
		return temp;
}


void Random::init(long seed)
{
	//	seed = 1332764272;
	//printf("seed =%d\n",seed);
	idum2=123456789;
   idum=0;
	iy=0;

	if (seed != 0)
		idum = seed;
   else
   	idum = 1;


	for (int j=NTAB+7;j>=0;j--) // Load the shuffle table (after 8 warm­ups).
	{
		long k=(idum)/IQ1;

		idum=IA1*(idum-k*IQ1)-k*IR1;
		if (idum < 0)
			idum += IM1;
		if (j < NTAB)
			iv[j] = idum;
	}
	iy=iv[0];
}

Random::Random(long seed)
{
	init(seed);
}

Random::Random()
{
   time_t t;

   time(&t);

   init((long)t);
}

double Random::doublerandom()
{
	double t = ran2();
   return t;
}

long Random::longrandom(long range)
{
	double t;

   t = doublerandom();
   return((long)(t*(double)range));
}

bool Random::boolrandom()
{
	double t=doublerandom();

   if (t>0.5)
   	return true;
   else
   	return false;
}
