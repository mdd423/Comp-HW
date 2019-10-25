
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
// dl'd from numerical recipe lib
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#define OK		   0		/* no problems */
#define	NOMEMORY	  -1		/* no memory available for buffers */
#define EQALXES		  -2		/* two x-coordinates are equal */

#define true		   1
#define false		   0
#define MAXTXTLEN	  80

// #define finit( fc , len ) { fc.a = malloc( ( len + 1 ) * sizeof( char ) ) ;  \
                            fc.l = len ; }

static int spline( float *x   , float *y   , int  n   , float  yp1 , float  ypn , float *y2 );
static int splint( float *xa , float *ya , float *y2a , int n    , float x   , float *y   );
float simpsons(float *k, float *power, float radius, float (&integrand)(float,float,float), int size);
int getLineNum();
void readIn( float *x_vals, float *pk_vals, int totLines );
float integrand(float k, float power, float r);

int main()
{
  int arraySize = getLineNum();

  float *xs = new float[arraySize];
  float *pks = new float[arraySize];
  readIn(xs, pks, arraySize);

  float *spline_vals =  new float[arraySize+1];
  spline(xs, pks, arraySize, 10e30, 10e30, spline_vals);

  std::ofstream powerFile;
  powerFile.open("powervals.txt");

  int k_divisions{ 100000 };
  ++k_divisions;
  float *k_array = new float[k_divisions];
  float *power_array = new float[k_divisions];

  // bounds on k to be integrated over
  float lowerbound_k = (xs[0]);
  float upperbound_k = (xs[arraySize-1]);
  for (int iii{ 0 }; iii < k_divisions; ++iii)
  {
    k_array[iii] = iii*(upperbound_k - lowerbound_k)/static_cast<float>(k_divisions-1) + lowerbound_k;
    // k_array[iii] = exp(k_array[iii]);
    splint(xs, pks, spline_vals, arraySize, k_array[iii], &power_array[iii]);
    powerFile << k_array[iii] << "\t" << power_array[iii] << std::endl;
    // k_array[iii] = log(k_array[iii]);
  }

  std::ofstream densityFile;
  densityFile.open("densevals.txt");

  int r_divisions{ 10000 };
  ++r_divisions;
  float *radii = new float[r_divisions+1];
  float *epsilon = new float[r_divisions+1];

  float r_min{ 0. };
  float r_max{ 120. };
  for (int iii{ 0 }; iii < r_divisions+1; ++iii)
  {
    radii[iii] = iii*(r_max - r_min)/static_cast<float>(r_divisions-1) + r_min;
    epsilon[iii] = simpsons(k_array, power_array, radii[iii], integrand, k_divisions);
    densityFile << radii[iii] << "\t" << epsilon[iii] << "\t" << radii[iii] * radii[iii] * epsilon[iii] << std::endl;
  }
  densityFile.close();

  delete[] xs;
  delete[] pks;

  delete[] spline_vals;

  delete[] radii;
  delete[] epsilon;

  delete[] k_array;
  delete[] power_array;
	return 0;
}

static int klo = -1 ;
static int khi = -1 ;

/*
    spline constructs a cubic spline given a set of x and y values, through
these values.
*/
static int spline( float *x   , float *y   , int  n   ,
                    float  yp1 , float  ypn , float *y2 )
{
	int  i,k;
  // u is a pointer to float but uninitialized
	float p,qn,sig,un,*u;
  // zeroes this pointer
  u=NULL;
  // allocates memory at pointers location that is the size of n-1 unsigned floats
	u=(float *)malloc((unsigned) (n-1)*sizeof(float));
        if (!u) return( NOMEMORY ) ;
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<=n-2;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free(u);
        return( OK ) ;
}
/*
    splint uses the cubic spline generated with spline to interpolate values
in the XY  table.
*/
static int splint( float *xa , float *ya , float *y2a ,
                    int n    , float x   , float *y   )
{
        int  r   = 0 ;
	int  k;
	float h,b,a;

        if ( klo < 0 ){
  	   klo=0;
 	   khi=n-1;
        } else {
           if ( x < xa[klo] ) klo=0;
           if ( x > xa[khi] ) khi=n-1;
        }
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) {
           //setfblank_c( y ) ;
           r = EQALXES ;
        } else {
	   a=(xa[khi]-x)/h;
	   b=(x-xa[klo])/h;
	   *y=a*ya[klo]+b*ya[khi]+( (a*a*a-a)*y2a[klo]+
                                    (b*b*b-b)*y2a[khi] ) * (h*h) / 6.0;
        }
	return( r ) ;
}
// assume even divisions in k
// radius remains constant
// power array corresponds to values in k array
// thus they must be the same size
// also int size must be even
float simpsons(float *k, float *power, float radius, float (&integrand)(float,float,float), int size)
{
	float integr{ 0. };
  float len = (k[size-1] - k[0])/static_cast<float>(size);
	integr += integrand(k[0],power[0],radius);
	integr += integrand(k[size-1],power[size-1],radius);
	for (int iii{ 1 }; iii < size; iii += 2)
	{
		integr += 4 * integrand(k[iii],power[iii],radius);
	}
	for (int iii{ 2 }; iii < size-1; iii += 2)
	{
		integr += 2 * integrand(k[iii],power[iii],radius);
	}
	return integr * len / 3.0;
}

int getLineNum()
{
  std::string line;
  std::ifstream myfile("lcdm_z0-2.matter_pk");

  int lineNum{ 0 };
  if (myfile.is_open())
  {
    int iii{ 0 };
    while( getline(myfile, line))
    {
      ++lineNum;
    }
    myfile.close();
  }
  else std::cout << "File unable to open." << std::endl;
  return lineNum;
}

void readIn( float *x_vals, float *pk_vals, int totLines )
{
  std::string line;
  std::ifstream myfile("lcdm_z0-2.matter_pk");
  if (myfile.is_open())
  {
    int iii{ 0 };
    while( getline(myfile, line))
    {
      std::stringstream theStream;
      theStream << line;
      theStream >> x_vals[iii];
      theStream >> pk_vals[iii];
      ++iii;
    }
    myfile.close();
  }
  else std::cout << "File unable to open." << std::endl;
}

float integrand(float k, float power, float r)
{
  return (1/(2* M_PI*M_PI))*pow(k,2) * power * sin(k*r)/(k*r);
}

// � g++ -std=c++0x densityfluc.cpp -o pk.exe
