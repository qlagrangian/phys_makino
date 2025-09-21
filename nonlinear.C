//
// nonlinear.C Copyright Jun Makino 1998
//
// demonstration program for numerical integration of
// harmonic oscillator
// uses a simple VECTOR class

#include  <stdlib.h>
#include  <math.h>
#include  <stdiostream.h>

#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

const int VLEN = 2;
#include  "vector.h"
double k1,k2;
const vector dxdt( vector & x)
{
    vector d;
    d[0] = x[1];
    d[1] = -(k2*x[0]*x[0]+k1)*x[0];
    return d;
}

void rk4(vector & x,
	 vfunc_ptr f,
	 double h)
{
static    vector kx1,kx2,kx3,kx4;
    kx1 = f(x)*h;
    kx2 = f(x+kx1*0.5)*h;
    kx3 = f(x+kx2*0.5)*h;
    kx4 = f(x+kx3)*h;
    x += (kx1+2*(kx2+kx3)+kx4)*0.166666666666666666666666667;
}

double energy( vector &x)
{
    double x2 = x[0]*x[0];
    return 0.5*x[1]*x[1] + 0.5*k1*x2 + 0.25*k2*x2*x2;
}
void leapfrog(vector & x,
	      vfunc_ptr f,
	      double h,
	      int & first_call)
{
    static vector fprev;
    if (first_call != 0){
	fprev = f(x);
	first_call = 0;
    }
    x[1] += h*fprev[1]*0.5;
    x[0] += h*x[1];
    fprev = f(x);
    x[1] += h*fprev[1]*0.5;
}

main()
{
    vector x;
    double  t,h ;
    x[0] = 1; x[1] = 0;
    int i,nsteps, ntime, print_interval;
    
    cerr << "Enter k1, k2:";
    cin >> k1 >> k2;
    
    cerr << "Enter steps,  time to stop, print_interval:";
    cin >> nsteps >> ntime >> print_interval;
    cout.precision(16);
    h = 1.0/nsteps;
    cerr  << " x= " << x  << " h = " <<  h <<"\n";
    cerr  << " tend= " << ntime <<"\n";
    cerr  << " k1, k2 =  " << k1  << " " << k2 << "\n";
        int integrator_type;
    cerr << "Enter integrator type (0: rk4, otherwize: leapfrog):";
    cin >> integrator_type;
    if (integrator_type) {
	cerr << "Leapfrog will be used\n";
    }else{
	cerr << "Rk4 will be used\n";
    }
    double e0 = energy(x);
    for(i=0; i<nsteps*ntime; i++ ){
	static int first_call = 1;
        t = h* (i+1);
	if (integrator_type) {
	    leapfrog(x,dxdt,h, first_call);
	}else{
	    rk4(x,dxdt,h);
	}
	if ((i+1) % print_interval == 0){
	    cerr << t << " " << x  << " " << energy(x) <<"\n";
	}
    } 
    cout << k1 << " " << k2 << " "   << nsteps << " ";
    cout << t << " " << x  << " " << energy(x) - e0 <<"\n";
}

