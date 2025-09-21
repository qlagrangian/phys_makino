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

#include "cpgplot.h"

#include  <stdlib.h>
#include  <math.h>
#include  <stdiostream.h>
const int VLEN = 4;
#include  "vector.h"

typedef void (*vintegfunc)(vector &, vfunc_ptr, double, int&);

#include "integrators.C"
#include "henonlib.C"

main()
{
    vector x;
    double  t,h ;
    x[0] = 0; x[1] = 0; x[2] = 0; x[3] = 0;
    int i,nsteps, ntime, print_interval;
    static vintegfunc integtab[4]={rk4,yoshida4,yoshida6b,gauss4};
    static char * integname[4] = {"rk4","yoshida4","yoshida6b","gauss4"};
    cerr << "Enter steps,  time to stop, print_interval:";
    cin >> nsteps >> ntime >> print_interval;
    cout.precision(16);
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
#ifndef NOGRAPHICS    
    initgraph(sqrt(h0*6),h0);
#endif    
    double theta;
    cerr << "Enter new theta (<0) to finish):";
    cin >> theta;
    while (theta >= 0){
	theta *= M_PI/180;
	x[0] = 0; x[1] = 0; x[2] = 0; x[3] = 0;
	x[1] = cos(theta)*sqrt(2*h0);
	x[3] = sin(theta)*sqrt(2*h0);
	double e0 = energy(x);
	int first_call = 1;
	for(i=0; i<nsteps*ntime; i++ ){
	    t = h* (i+1);
	    vector xold = x;
	    integtab[integrator_type](x,dxdt,h, first_call);
	    if ((i+1) % print_interval == 0){
		cerr << t << " " << x  << " " << energy(x) <<"\n";
	    }
	    //	    cerr <<"Before Plot "; PRC(x);PRL(energy(x));
#ifndef NOGRAPHICS    
	    first_call |=plot(x,xold,dxdt,h);
#endif	    
	    //	    PRL(first_call);
	    
	    //	    cerr <<"After  Plot ";PRC(x);PRL(energy(x));
	    //	    plot(x,xold);
	} 
	cout << e0 <<" " << theta << " " << nsteps << " " << integrator_type << " ";
	cout << t << " " << x  << " " << energy(x) - e0 <<"\n";
	cerr << "Enter new theta (<0) to finish):";
	cin >> theta;
    }
#ifndef NOGRAPHICS    
    cpgend();
#endif    
}

