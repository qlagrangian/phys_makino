//
// henon.C Copyright Jun Makino 1998
//
// demonstration program for numerical integration of
// Henon-Heiles Hamiltonian
// H = \frac{1}{2}(p_1^2 + p_2^2) +  \frac{1}{2}(q_1^2 + q_2^2) +
// q_1^2q_2 - \frac{1}{3}q_2^3
//
// dp_1/dt= -q_1 - 2q_1q_2
// dp_2/dt= -q_2 - q_1^2 + q_2^2
// uses a simple VECTOR class
// x_0 = q_1
// x_1 = p_1
// x_2 = q_2
// x_3 = p_2
//
// Similar to HENON2 but integrate two orbits with small initial distance
//

#include "cpgplot.h"

#include  <stdlib.h>
#include  <math.h>
#include  <stdiostream.h>
const int VLEN = 4;
#include  "vector.h"

#include "integrators.C"
#include "henonlib.C"

main()
{
    vector x, x2;
    double  t,h ;
    int i,nsteps, ntime, print_interval;
    static vintegfunc integtab[4]={rk4,yoshida4,yoshida6b,gauss4};
    static char * integname[4] = {"rk4","yoshida4","yoshida6b","gauss4"};
    cerr << "Enter steps,  time to stop, print_interval:";
    cin >> nsteps >> ntime >> print_interval;
    cout.precision(6);
    h = 1.0/nsteps;
    cerr  << " x= " << x  << " h = " <<  h <<"\n";
    cerr  << " tend= " << ntime <<"\n";
    int integrator_type;
    cerr << "Enter integrator type (0: rk4, 1: yoshida4, 2: yoshida6b,3:gauss4):";
    cin >> integrator_type;
    cerr << integname[integrator_type] << " will be used\n";
    double h0;
    cerr << "Enter h0 (=p1) ";
    cin >> h0;
    double theta, dtheta;
    cerr << "Enter theta, dtheta(note: dtheta is radian):";
    cin >> theta >> dtheta;
    theta *= M_PI/180;
    x =x2 = 0; 
    x[1] = cos(theta)*sqrt(2*h0);
    x[3] = sin(theta)*sqrt(2*h0);
    x2[1] = cos(theta+dtheta)*sqrt(2*h0);
    x2[3] = sin(theta+dtheta)*sqrt(2*h0);
    double e0 = energy(x);
    int first_call = 1;
    int first_call2 = 1;
    for(i=0; i<nsteps*ntime; i++ ){
	t = h* (i+1);
	integtab[integrator_type](x,dxdt,h, first_call);
	integtab[integrator_type](x2,dxdt,h, first_call2);
	first_call = first_call2 = 1;
	if ((i+1) % print_interval == 0){
	    cout << t << " " << x  << " " << energy(x) - e0 ;
	    cout << " "<< sqrt((x-x2)*(x-x2)) <<"\n";
	}
    } 
}

