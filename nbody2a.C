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
#ifndef NBODY
const int NBODY = 2;
#endif
const int NDIM = 3;
const int VLEN = NBODY*NDIM*2;
#include  "vector.h"
typedef void (*vintegfunc)(vector &, vfunc_ptr, double, int&);

#include "integrators.C"


const double EPS = 1e-10;



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
#ifndef TESTMAIN2
#ifndef TESTMAIN3
#define STANDARD_MAIN
#endif
#endif
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
	    cout << t << " " <<  e << " " << e-e0 <<"\n";
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

void setup_2body(vector & x, double ecc)
{
    if (NBODY != 2){
	cerr << "setup_2body NBODY = " << NBODY << "unspoorted" <<endl;
	exit(1);
    }
    for(int i = 0; i<VLEN; i++){
	x[i] = 0;
    }
    x[0] = 1; x[6] = -x[0];
    double vx = ecc * 0.5;
    double vy = sqrt(1-ecc*ecc)*0.5;
    x[1] = vx; x[7] = -x[1];
    x[3] = vy; x[9] = -x[3];
}

main()
{
    vector x;
    double  t, h;
    int   norbits ;
    cerr << "Enter norbits, eps2:";
    cin >> norbits >>   eps2;
    static vintegfunc integtab[]={rk4,leapfrog,yoshida4,yoshida6b,gauss4};
    static char * integname[] = {"rk4","leapfrog","yoshida4","yoshida6b","gauss4"};
    cerr << "Enter steps,   print_interval:";
    int nsteps, print_interval;
    cin >> nsteps >>  print_interval;
    double ecc;
    cerr << "enter ecc"; cin >> ecc;
    setup_2body(x,ecc);
    cout.precision(6);
    h = 4*M_PI/nsteps;
    cerr  << " x= " << x  << " h = " <<  h <<"\n";
    cerr  << " norbits= " << norbits <<"\n";
    int integrator_type;
    cerr << "Enter integrator type (0: rk4, 1: leapfrog, 2: yoshida4, 3: yoshida6b,4:gauss4):";
    cin >> integrator_type;
    cerr << integname[integrator_type] << " will be used\n";
    cerr << "Enter xmax:";
    double xmax;
    cin >> xmax;
    double e0;
    e0 = energy(x);
    vector x0 = x;
    cerr << "Energy " << e0  <<"\n";
    int first_call = 1;
    int ntotalsteps = nsteps*norbits;
    for(int istep = 0; istep<ntotalsteps; istep++){
	t = istep*h;
	integtab[integrator_type](x,dxdt,h,first_call);
	if ( (istep % print_interval) == 0 || istep == ntotalsteps-1){
	    double e;
	    e = energy(x);
	    cerr << t << " " <<  e << " " << e-e0 <<"\n";
	    PRL(x);
	}  
    }
    double e;
    e = energy(x);
    cout << norbits << " " << nsteps << " " << integname[integrator_type] << " " << ecc  << " " <<e << " " << e-e0 << " ";
    cout << sqrt((x-x0)*(x-x0))<< endl;
    

}

#endif

#ifdef TESTMAIN2
//
// test main routine for error plot
//

void setup_2body(vector & x, double ecc)
{
    if (NBODY != 2){
	cerr << "setup_2body NBODY = " << NBODY << "unspoorted" <<endl;
	exit(1);
    }
    for(int i = 0; i<VLEN; i++){
	x[i] = 0;
    }
    
    x[0] = 1; x[6] = -x[0];
    double vx = ecc * 0.5;
    double vy = sqrt(1-ecc*ecc)*0.5;
    x[1] = vx; x[7] = -x[1];
    x[3] = vy; x[9] = -x[3];

    double xfact = pow(4*M_PI,-2.0/3);
    double vfact = pow(4*M_PI,1.0/3);
    for(int i = 0; i<VLEN; i+= 2) x[i] *= xfact;
    for(int i = 1; i<VLEN; i+= 2) x[i] *= vfact;
}

main()
{
    vector x;
    double  t, h, hmax;
    int   norbits ;
    cerr << "Enter norbits, eps2:";
    cin >> norbits >>   eps2;
    static vintegfunc integtab[]={rk4,leapfrog,yoshida4,yoshida6b,gauss4};
    static char * integname[] = {"rk4","leapfrog","yoshida4","yoshida6b","gauss4"};
    cerr << "Enter steps,  print_interval, tol, tol_fact,  :";
    int nsteps, print_interval;
    double tol, tol_fact;
    cin >> nsteps >>  print_interval >> tol >> tol_fact;
    double ecc;
    cerr << "enter ecc"; cin >> ecc;
    setup_2body(x,ecc);
    cout.precision(6);
    h = hmax = 1.0/nsteps;
    cerr  << " x= " << x  << " h = " <<  h <<"\n";
    cerr  << " norbits= " << norbits <<"\n";
    int integrator_type;
    cerr << "Enter integrator type (0: rk4, 1: leapfrog, 2: yoshida4, 3: yoshida6b,4:gauss4):";
    cin >> integrator_type;
    cerr << integname[integrator_type] << " will be used\n";
    double e0;
    e0 = energy(x);
    vector x0 = x;
    cerr << "Energy " << e0  <<"\n";
    int first_call = 1;
    int stepcount = 0;
    int tout = print_interval;
    for(t = 0; t < norbits; ){
	t += h;
	variable_step_integrator(x,integtab[integrator_type],
				 dxdt,h,hmax, tol, tol_fact);
	while (fmod(t,h) != 0.0) {
	    h /= 2;
	}
	if ( t>= tout ){
	    double e;
	    e = energy(x);
	    PRC(t); PRC(e); PRL( e-e0);
	    PRL(x);
	    tout += print_interval;
	}
	stepcount ++;
    }
    double e;
    e = energy(x);
    cout << t << " " << tol << " " << integname[integrator_type] << " " << ecc  << " " << stepcount << " " <<e << " " << e-e0 << " ";
    cout << sqrt((x-x0)*(x-x0))<< endl;
    

}

#endif



#ifdef TESTMAIN3
//
// DEMO program
//

main()
{
    vector x;
    double  t, h, hmax;
    int   tmax ;
    cerr << "Enter tmax";
    cin >> tmax;
    cerr  << " tmax= " << tmax <<"\n";
    eps2 = 0;
    static vintegfunc integtab[]={rk4,leapfrog,yoshida4,yoshida6b,gauss4};
    static char * integname[] = {"rk4","leapfrog","yoshida4","yoshida6b","gauss4"};
    cerr << "Enter steps,  print_interval, tol, tol_fact,  :";
    int nsteps, print_interval;
    double tol, tol_fact;
    cin >> nsteps >>  print_interval >> tol >> tol_fact;
    int integrator_type;
    cerr << "Enter integrator type (0: rk4, 1: leapfrog, 2: yoshida4, 3: yoshida6b,4:gauss4):";
    cin >> integrator_type;
    cerr << integname[integrator_type] << " will be used\n";
    initgraph(10);
    getbody(x);
    cout.precision(6);
    h = hmax = 1.0/nsteps;
    cerr  << " x= " << x  << " h = " <<  h <<"\n";
    double e0;
    e0 = energy(x);
    vector x0 = x;
    cerr << "Energy " << e0  <<"\n";
    int first_call = 1;
    int stepcount = 0;
    int tout = print_interval;
    for(t = 0; t < tmax; ){
	t += h;
	variable_step_integrator(x,integtab[integrator_type],
				 dxdt,h,hmax, tol, tol_fact);
	while (fmod(t,h) != 0.0) {
	    h /= 2;
	}
	plot(x);
	if ( t>= tout ){
	    double e;
	    e = energy(x);
	    PRC(t); PRC(e); PRL( e-e0);
	    PRL(x);
	    tout += print_interval;
	}
	stepcount ++;
    }
    double e;
    e = energy(x);
    cout << t << " " << tol << " " << integname[integrator_type] << stepcount << " " <<e << " " << e-e0 << " ";
    

}

#endif



