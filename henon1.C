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

#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

typedef struct plot_control{
    int xvarid;
    int yvarid;
    int secvarid;
    int mode;
}PLOT_CONTROL;

struct plot_control pctr;

void initgraph(const double xmax, const double h0)
{
    cerr << "Enter var to plot (x, y, and sos):";
    cin >> pctr.xvarid >>  pctr.yvarid >> pctr.secvarid ;
    if ( pctr.secvarid < 0){
	pctr.mode = 0;
    }else{
	pctr.mode = 1;
    }
    char  pgdevname[255];
    cerr << "Enter devname";
    cin >> pgdevname ;
    
    if(cpgopen(pgdevname) != 1) exit(EXIT_FAILURE);
    cpgask(0);
    cpgenv(-xmax, xmax, -xmax, xmax, 1, 0);
    const char * name[4]={"q1","p1","q2","p2"};
    char label[255];
    sprintf(label, "H0 = %g", h0);
    cpglab(name[pctr.xvarid],name[pctr.yvarid], label);
}

#ifndef WRONG_PLOT
void plot( vector x,  vector xold)
{

    if (pctr.mode == 0 ){
	cpgpt1(x[pctr.xvarid],x[pctr.yvarid],-1);
    }else    if (x[pctr.secvarid]*xold[pctr.secvarid] <= 0) {
	vector xinterp = (-x*xold[pctr.secvarid]+xold*x[pctr.secvarid])
	    /(x[pctr.secvarid]-xold[pctr.secvarid]);
	cpgpt1(xinterp[pctr.xvarid],xinterp[pctr.yvarid],-1);
    }	


}
#else
void plot( vector x,  vector xold)
{

    if (pctr.mode == 0 ){
	cpgpt1(x[pctr.xvarid],x[pctr.yvarid],-1);
    }else    if (x[pctr.secvarid]*xold[pctr.secvarid] <= 0) {
	cpgpt1(x[pctr.xvarid],x[pctr.yvarid],-1);
    }	


}
#endif

const vector dxdt( vector & x)
{
    vector d;
    d[0] = x[1];
    d[1] = -x[0] - 2* x[0]*x[2];
    d[2] = x[3];
    d[3] = -x[2] - x[0]*x[0] + x[2]*x[2];
    return d;
}


double energy( vector &x)
{
    return 0.5*(x[1]*x[1]+x[3]*x[3]) + 0.5*(x[0]*x[0]+x[2]*x[2])
	+ x[0]*x[0]*x[2]-x[2]*x[2]*x[2]/3;
}

void rk4(vector & x,
	 vfunc_ptr f,
	 double h,
	 int & first_call)
{
static    vector kx1,kx2,kx3,kx4;
    kx1 = f(x)*h;
    kx2 = f(x+kx1*0.5)*h;
    kx3 = f(x+kx2*0.5)*h;
    kx4 = f(x+kx3)*h;
    x += (kx1+2*(kx2+kx3)+kx4)*0.166666666666666666666666667;
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
    for(int i = 0; i<4;i+=2){
	x[1+i] += h*fprev[1+i]*0.5;
	x[i] += h*x[1+i];
    }
    fprev = f(x);
    for(int i = 0; i<4;i+=2){
	x[1+i] += h*fprev[1+i]*0.5;
    }
}

const int niter = 6;
void trapezoidal(vector & x,
	      vfunc_ptr f,
	      double h,
	      int & first_call)
{
    static vector fprev;
    if (first_call != 0){
	fprev = f(x);
	first_call = 0;
    }
    vector xnew = x + fprev*h;
    vector fnew;
    for(int i = 0; i<niter; i++){
	fnew = f(xnew);
	xnew = x + (fprev+fnew)*(h*0.5);
	//	PRC(i);PRC(xnew);PRL(fnew);
    }
    fprev = fnew;
    x =xnew;
}

void gauss4(vector & x,
	      vfunc_ptr f,
	      double h,
	    int & first_call)
{
    vector f1,f2;
    f1 = f2 = f(x);
    double a11, a12, a21, a22;
    a11 = 0.25;
    a12 = 0.25 - sqrt(3)/6;
    a21 = 0.25 + sqrt(3)/6;
    a22 = 0.25;
    for(int i = 0; i<niter; i++){
	vector xg1 = x+(a11*f1+a12*f2)*h;
	vector xg2 = x+(a21*f1+a22*f2)*h;
	f1 = f(xg1);
	f2 = f(xg2);
    }
    x = x+(f1+f2)*0.5*h;
}


void yoshida4(vector & x,
	      vfunc_ptr f,
	      double h,
	      int & first_call)
{
    static double d1, d2;
    if (first_call != 0){
	d1 = 1.0 / (2-pow(2,1.0/3.0));
	d2 = 1 - 2*d1;
	PRC(d1); PRL(d2);
    }
    leapfrog(x,dxdt,h*d1, first_call);
    leapfrog(x,dxdt,h*d2, first_call);
    leapfrog(x,dxdt,h*d1, first_call);
}
void yoshida6b(vector & x,
	      vfunc_ptr f,
	      double h,
	      int & first_call)
{
    const double d[4] = {0.784513610477560e0, 0.235573213359357e0,
			  -1.17767998417887e0, 1.31518632068391e0};
    for(int i=0;i<3;i++)leapfrog(x,dxdt,h*d[i], first_call);
    leapfrog(x,dxdt,h*d[3], first_call);
    for(int i=0;i<3;i++)leapfrog(x,dxdt,h*d[2-i], first_call);
}

typedef void (*vintegfunc)(vector &, vfunc_ptr, double, int&); 

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
    initgraph(sqrt(h0*6),h0);
    
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
	    plot(x,xold);
	} 
	cout << e0 <<" "  << nsteps << " ";
	cout << t << " " << x  << " " << energy(x) - e0 <<"\n";
	cerr << "Enter new theta (<0) to finish):";
	cin >> theta;
    }
    cpgend();
}

