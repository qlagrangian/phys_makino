//
// harmonic4.C Copyright Jun Makino 1998
//
// demonstration program for numerical integration of
// harmonic oscillator

#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <stdiostream.h>
double k;
double dvdt(double x)
{
    return -k*x;
}
void rk4(double & x,
	 double & v,
	 double h)
{
    double kx1,kx2,kx3,kx4;
    double kv1,kv2,kv3,kv4;
    kx1 = v*h;
    kv1 = dvdt(x)*h;
    kx2 = (v+kv1/2)*h;
    kv2 = dvdt(x+kx1/2)*h;
    kx3 = (v+kv2/2)*h;
    kv3 = dvdt(x+kx2/2)*h;
    kx4 = (v+kv3)*h;
    kv4 = dvdt(x+kx3)*h;
    x += (kx1+2*kx2+2*kx3+kx4)/6;
    v += (kv1+2*kv2+2*kv3+kv4)/6;
}
void leapfrog(double & x,
	      double & v,
	      double h)
{
    double vhalf;
    vhalf = v + dvdt(x)*h/2;
    x += vhalf*h;
    v = vhalf + dvdt(x)*h/2;
}

main()
{
    double x, v, t,h ;
    k = 1; x = 1; v = 0;
    int i,nsteps, norbits, print_interval;
    cerr << "Enter steps,  orbits, print_interval:";
    cin >> nsteps >> norbits >> print_interval;
    h = 2*M_PI/nsteps;
    cout.precision(16);
    cerr  << " x= " << x << " v= " << v << " h = " <<  h <<"\n";
    cerr  << " norbits= " << norbits <<"\n";
    int integrator_type;
    cerr << "Enter integrator type (0: rk4, otherwize: leapfrog):";
    cin >> integrator_type;
    if (integrator_type) {
	cerr << "Leapfrog will be used\n";
    }else{
	cerr << "Rk4 will be used\n";
    }
    for(i=0; i<nsteps*norbits; i++ ){
	double phase = h * (i % nsteps + 1);
        t = h* (i+1);
	if (integrator_type) {
	    leapfrog(x,v,h);
	}else{
	    rk4(x,v,h);
	}
	if ( ((i+1) % print_interval) == 0){
	    cout << "t= " << t  << "           err_x= " << x-cos(phase)
	    << " err_v= " << v+sin(phase) << "\n";
	}
    }
}

