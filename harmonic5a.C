#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <stdiostream.h>

const int VLEN = 2;
#include  "vector.h"

double k;
vector dxdt(vector & x)
{
    vector d;
    d[0] = x[1];
    d[1] = -k*x[0];
    return d;
}
void rk4(vector & x,
	 vfunc f,
	 double h)
{
    vector kx1,kx2,kx3,kx4;
    kx1 = f(x)*h;
    kx2 = f(x+kx1/2)*h;
    kx3 = f(x+kx2/2)*h;
    kx4 = f(x+kx3)*h;
    x += (kx1+2*kx2+2*kx3+kx4)/6;
}

void main()
{
    vector x;
    double  t,h ;
    k = 1; x[0] = 1; x[1] = 0;
    int i,nsteps, norbits, print_interval;
    cerr << "Enter steps,  orbits, print_interval:";
    cin >> nsteps >> norbits >> print_interval;
    h = 2*M_PI/nsteps;
    cout.precision(16);
    cerr  << " x= " << x  << " h = " <<  h <<"\n";
    cerr  << " norbits= " << norbits <<"\n";
    double emax = 0;
    for(i=0; i<nsteps*norbits; i++ ){
	vector exact;
	double phase = h * (i % nsteps + 1);
        t = h* (i+1);
	rk4(x,dxdt,h);
        if ( ((i+1) % print_interval) == 0){
            cout << "t= " << t  << "           err_x= " << x[0]-cos(phase)
                 << " err_v= " << x[1]+sin(phase) << "\n";
        }
    }
}
