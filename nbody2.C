//
// nbody2.C Copyright Jun Makino 1998
//
// phase space coordinates are expressed in a single big vector
//
//

#include "cpgplot.h"

#include  <stdlib.h>
#include  <math.h>
#include  <stdiostream.h>

#ifndef N_BODY
const int NBODY = 2;
#else
const int NBODY = N_BODY;
#endif


const int NDIM = 3;
const int VLEN = NBODY*NDIM*2;
#include  "vector.h"

typedef void (*vintegfunc)(vector &, vfunc_ptr, double, int&);
const double EPS = 1e-10;

#include "integrators.C"




typedef void (*vintegfunc)(vector &, vfunc_ptr, double, int&);


typedef struct plot_control{
    int xvarid;
    int yvarid;
    int mode;
    double xmax;
}PLOT_CONTROL;

struct plot_control pctr;

void initgraph(const double xmax)
{
    cerr << "Enter var to plot (x, y):";
    cin >> pctr.xvarid >>  pctr.yvarid ;
    cerr << "Enter mode (0: snapplot, 1: trakplot):";
    cin >> pctr.mode ;
    char  pgdevname[255];
    cerr << "Enter devname";
    cin >> pgdevname ;
    if(cpgopen(pgdevname) != 1) exit(EXIT_FAILURE);
    cpgask(0);
    cpgenv(-xmax, xmax, -xmax, xmax, 1, 0);
    pctr.xmax = xmax;
    const char * name[6]={"x","vx","y","vy","z","vz"};
    char label[255];
    cpglab(name[pctr.xvarid],name[pctr.yvarid], label);
    if(pctr.mode == 0){
	cpgslw(10);
    }
}

inline real & pos(vector &x, int i, int k)       {return  x[6*i+2*k];}
inline real & vel(vector &x, int i, int k)       {return  x[6*i+2*k+1];}
inline real & acc(vector &x, int i, int k)       {return  x[6*i+2*k+1];}

double eps2;

calc_force(vector & x, vector & d, int i, int j)
{
    int k;
    double dx[NDIM];
    double r2 = eps2;
    //    PRC(i);PRL(j);
    for(k=0;k<NDIM;k++) {
	dx[k] = pos(x,i,k)-pos(x,j,k);
	r2 += dx[k]*dx[k];
	//	PRC(k),PRC(pos(x,i,k)); PRC(pos(x,j,k)); PRL(dx[k]);
    }
    double r3inv = 1/(sqrt(r2)*r2);
    for(k=0;k<NDIM;k++) {
	acc(d,i,k) -= dx[k]*r3inv;
	acc(d,j,k) += dx[k]*r3inv;
    }
}

double pair_energy(vector & x,  int i, int j)
{
    int k;
    double dx[NDIM];
    double r2 = eps2;
    for(k=0;k<NDIM;k++) {
	dx[k] = pos(x,i,k)-pos(x,j,k);
	r2 += dx[k]*dx[k];
    }
    return - 1/sqrt(r2);
}




const vector dxdt( vector & x)
{
    vector d;
    int i,j;
    for(i = 0; i<VLEN; i+= 2){
	d[i] = x[i+1];
	d[i+1] = 0;
    }
    for(i=0;i<NBODY-1; i++){
	for(j=i+1;j<NBODY; j++){
	    calc_force(x,d,i,j);
	}
    }
    return d;
}


double energy( vector &x)
{
    vector d;
    int i,j;
    double ke = 0;
    double  pe = 0;
    for(i = 0; i<VLEN; i+= 2){
	ke += x[i+1]* x[i+1];
    }
    ke *= 0.5;
    for(i=0;i<NBODY-1; i++){
	for(j=i+1;j<NBODY; j++){
	    pe += pair_energy(x,i,j);
	}
    }
    return ke+pe;
}


void plot(vector &x)
{
    int i;
    if(pctr.mode == 0){
	double xmax = pctr.xmax;
	cpgeras();  
    } 
    
    for(i=0; i < NBODY; i++){
	cpgsci( i % 15 + 1);
	cpgpt1(x[i*6+pctr.xvarid],x[i*6+pctr.yvarid],-1);
    }
}

void getbody (vector &x)
{
    int n;
    cin >> x;
}


void dump(vector & x)
{
    cerr << "N = " << NBODY << " eps2 = " << eps2 << "\n";
    cerr << x << "\n"; 
}

#ifndef TESTMAIN
#define STANDARD_MAIN
#endif

#ifdef STANDARD_MAIN
main()
{
    vector x;
    double  t, toutstep,  tend, h ;
    cerr << "Enter tend,toustep, h, eps2:";
    cin >> tend >> toutstep >> h >> eps2;

    PRC(tend); PRC(h); PRL(eps2);
    cerr << "Enter xmax:";
    double xmax;
    cin >> xmax;
    initgraph(xmax);
    getbody(x);
    plot(x);
    double e0;
    e0 = energy(x);
    cout << "Energy " << e0  <<"\n";
    dump(x);
    plot(x);
    int first_call = 1;
    for(t=0; t<tend; t+=h ){
	yoshida6b(x,dxdt,h,first_call); 
	if (fmod(t,toutstep) == 0.0|| t == tend){
	    double e;
	    e = energy(x);
	    cout << "TIME = "<<t << " " <<  e << " " << e-e0 <<"\n";
	    dump(x);
	}  
	plot(x);
    }
    int dum;
    cin >> dum;
    cpgend();

}

#endif

#ifdef TESTMAIN
//
// test main routine for error plot
//

main()
{
    vector x;
    double  t, toutstep,  tend, h ;
    cerr << "Enter tend,toustep, h, eps2:";
    cin >> tend >> toutstep >> h >> eps2;

    PRC(tend); PRC(h); PRL(eps2);
    cerr << "Enter xmax:";
    double xmax;
    cin >> xmax;
    initgraph(xmax);
    getbody(x);
    plot(x);
    double e0;
    e0 = energy(x);
    cout << "Energy " << e0  <<"\n";
    dump(x);
    plot(x);
    int first_call = 1;
    for(t=0; t<tend; t+=h ){
	yoshida6b(x,dxdt,h,first_call); 
	if (fmod(t,toutstep) == 0.0|| t == tend){
	    double e;
	    e = energy(x);
	    cout << "Time = " <<t << " " <<  e << " DE=" << e-e0 <<"\n";
	    dump(x);
	}  
	plot(x);
    }
    int dum;
    cin >> dum;
    cpgend();

}

#endif

