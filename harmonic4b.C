//
// harmonic4b.C Copyright Jun Makino 1998
//
// demonstration program for numerical integration of
// harmonic oscillator
// use C type array

#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <stdiostream.h>

#define NDIM 2
double k;

void calc_df(double x[NDIM], double f[NDIM])
{
    f[0] = x[1];
    f[1] =  -k*x[0];
}
void rk4(double x[NDIM],
	 double h)
{
    double kx1[NDIM],kx2[NDIM],kx3[NDIM],kx4[NDIM], xtmp[NDIM],ftmp[NDIM];
    int i;

    calc_df(x,ftmp);
    for(i=0;i<NDIM;i++)kx1[i] = ftmp[i]*h;

    for(i=0;i<NDIM;i++)xtmp[i] = x[i]+kx1[i]/2;
    calc_df(xtmp,ftmp);
    for(i=0;i<NDIM;i++)kx2[i] = ftmp[i]*h;
    
    for(i=0;i<NDIM;i++)xtmp[i] = x[i]+kx2[i]/2;
    calc_df(xtmp,ftmp);
    for(i=0;i<NDIM;i++)kx3[i] = ftmp[i]*h;
    
    for(i=0;i<NDIM;i++)xtmp[i] = x[i]+kx3[i];
    calc_df(xtmp,ftmp);
    for(i=0;i<NDIM;i++)kx4[i] = ftmp[i]*h;
    
    for(i=0;i<NDIM;i++)x[i] += (kx1[i]+2*kx2[i]+2*kx3[i]+kx4[i])/6;
}


main()
{
    double x[NDIM], t,h ;
    k = 1; x[0] = 1; x[1] = 0;
    int i,nsteps, norbits, print_interval;
    cerr << "Enter steps,  orbits, print_interval:";
    cin >> nsteps >> norbits >> print_interval;
    h = 2*M_PI/nsteps;
    cout.precision(16);
    cerr  << " x= " << x[0] << " v= " << x[1] << " h = " <<  h <<"\n";
    cerr  << " norbits= " << norbits <<"\n";
    for(i=0; i<nsteps*norbits; i++ ){
	double phase = h * (i % nsteps + 1);
        t = h* (i+1);
	rk4(x,h);
	if ( ((i+1) % print_interval) == 0){
	    cout << "t= " << t  << "           err_x= " << x[0]-cos(phase)
		 << " err_v= " << x[1]+sin(phase) << "\n";
	}
    }
}

