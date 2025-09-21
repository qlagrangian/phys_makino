//
// nbody.C Copyright Jun Makino 1998
//
// demonstration program for numerical integration of
//
// particles are treated as objects
//

#include "cpgplot.h"

#include  <stdlib.h>
#include  <math.h>
#include  <stdiostream.h>
const int NDIM =3;
const int VLEN = NDIM;
#include  "nbody.h"

const int niter = 5;

typedef void (*vintegfunc)(nbody *,  double, int&);
const double EPS = 1e-10;


#ifndef NOGRAPHICS
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
    const char * name[6]={"x","y","z","vx","vy","vz"};
    char label[255];
    cpglab(name[pctr.xvarid],name[pctr.yvarid], label);
    if(pctr.mode == 0){
	cpgslw(10);
    }
}
#else
void initgraph(const double xmax){}

#endif
#if 0
void calc_and_acc_force_between_two_particles(body * p1,
					      body * p2,
					      double eps2)
{
    vector dx = p1->pos-p2->pos;
    double r2inv = 1/(dx*dx+eps2);
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;
    p1->phi -= p2->mass*rinv;
    p2->phi -= p1->mass*rinv;
    p1->acc -= p2->mass*r3inv*dx;
    p2->acc += p1->mass*r3inv*dx;
}
#endif
#if 0
void calc_and_acc_force_between_two_particles(body * p1,
					      body * p2,
					      double eps2)
{
    double dx[3];
#if 0
    double r2 = eps2;
    for(int i = 0; i<3;i++){
	dx[i] = p1->pos[i]-p2->pos[i];
	r2 += dx[i]*dx[i];
    }
#else
    dx[0] = p1->pos[0]-p2->pos[0];
    dx[1] = p1->pos[1]-p2->pos[1];
    dx[2] = p1->pos[2]-p2->pos[2];
    double r2 = eps2 + dx[0]*dx[0]+ dx[1]*dx[1]+ dx[2]*dx[2];
#endif    
    double r2inv = 1/r2;
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;
    p1->phi -= p2->mass*rinv;
    p2->phi -= p1->mass*rinv;
    double ifact = p2->mass*r3inv;
    double jfact = p1->mass*r3inv;
#if 0    
    for(int i = 0; i<3;i++){
	p1->acc[i] -= ifact*dx[i];
	p2->acc[i] += jfact*dx[i];
    }
#else
	p1->acc[0] -= ifact*dx[0];
	p1->acc[1] -= ifact*dx[1];
	p1->acc[2] -= ifact*dx[2];
	p2->acc[0] += jfact*dx[0];
	p2->acc[1] += jfact*dx[1];
	p2->acc[2] += jfact*dx[2];
#endif
}

#endif
#if 1
extern "C"{
    
void force_(double * xi,
	    double * xj,
	    double * mi,
	    double * mj,
	    double * ai,
	    double * aj,
	    double * phii,
	    double * phij,
	    double * eps2);
}

inline void calc_and_acc_force_between_two_particles(body * p1,
					      body * p2,
					      double eps2)
{
    
    force_((double*)&(p1->pos),(double*)&(p2->pos),
	   &(p1->mass),	   &(p2->mass),
	   (double*)&(p1->acc),(double*)&(p2->acc),
	   &(p1->phi),	   &(p2->phi),
	   &eps2);
}


#endif

calc_force_for_system(nbody* nbp)
{
    int i, j;
    body * bpi,* bpj;
    for(i=0, bpi= nbp->bp; i < nbp->n; i++,bpi++)bpi->clear_acc_and_phi();
    for(i=0, bpi= nbp->bp; i < nbp->n - 1; i++,bpi++){
	for(j=i+1, bpj = bpi+1; j< nbp->n; j++, bpj++){
	    calc_and_acc_force_between_two_particles
		(bpi, bpj,nbp->eps2);
	}
    }
}

push_pos(nbody * nbp, double h)
{
    int i;
    body * bpi;
    for(i=0, bpi= nbp->bp; i < nbp->n; i++,bpi++)bpi->push_pos(h);
}
push_vel(nbody * nbp, double h)
{
    int i;
    body * bpi;
    for(i=0, bpi= nbp->bp; i < nbp->n; i++,bpi++)bpi->push_vel(h);
}
    

void leapfrog(nbody * nbp,
	      double h,
	      int & first_call)
{
    if (first_call != 0){
	calc_force_for_system(nbp);
	first_call = 0;
    }
    push_vel(nbp, h/2);
    push_pos(nbp,h);
    calc_force_for_system(nbp);
    push_vel(nbp, h/2);
}

void yoshida4(nbody * nbp,
	      double h,
	      int & first_call)
{
    static double d1, d2;
    if (first_call != 0){
	d1 = 1.0 / (2-pow(2,1.0/3.0));
	d2 = 1 - 2*d1;
    }
    leapfrog(nbp,h*d1, first_call);
    leapfrog(nbp,h*d2, first_call);
    leapfrog(nbp,h*d1, first_call);
}
void yoshida6b(nbody * nbp,
	      double h,
	      int & first_call)
{
    const double d[4] = {0.784513610477560e0, 0.235573213359357e0,
			  -1.17767998417887e0, 1.31518632068391e0};
    for(int i=0;i<3;i++)leapfrog(nbp,h*d[i], first_call);
    leapfrog(nbp,h*d[3], first_call);
    for(int i=0;i<3;i++)leapfrog(nbp,h*d[2-i], first_call);
}


void body::trapezoidal_iterate(body & prev_val, double h)
{
    pos = prev_val.pos + (h*0.5)*(vel + prev_val.vel);
    vel = prev_val.vel + (h*0.5)*(acc + prev_val.acc);
}



void trapezoidal(nbody * nbp,
	      double h,
	      int & first_call)
{
    static nbody nbprev;
    nbody * nbprevp = & nbprev;
    if (first_call != 0){
	calc_force_for_system(nbp);
	nbprevp->n = nbp->n;
	nbprevp->eps2 = nbp->eps2;
	nbprevp->bp  = new body[nbp->n];
	first_call = 0;
    }
    int i;
    body * bpi;
    body * bpj;
    for(i=0, bpi= nbp->bp, bpj = nbprevp->bp; i < nbp->n; i++,bpi++,bpj++){
	*bpj = *bpi;
    }
    for(int k=0;k < niter; k++){
	for(i=0, bpi= nbp->bp, bpj = nbprevp->bp; i < nbp->n; i++,bpi++,bpj++){
	    bpi->trapezoidal_iterate(*bpj, h);
	}
	calc_force_for_system(nbp);
    }

}


void calc_energy(nbody * nbp, double & ke, double & pe, double& etot)
{
    int i;
    body * bpi;
    ke= 0;
    pe= 0;
    for(i=0, bpi= nbp->bp; i < nbp->n; i++,bpi++){
	ke += bpi->kinetic_energy();
	pe += bpi->potential_energy();
    }
    pe *= 0.5;
    etot = pe+ke;
}

#ifndef NOGRAPHICS
void plot(nbody * nbp)
{
    int i;
    body * bpi;
    if(pctr.mode == 0){
	double xmax = pctr.xmax;
	cpgeras();  
    } 
    
    for(i=0, bpi= nbp->bp; i < nbp->n; i++,bpi++){
	cpgsci( i % 15 + 1);
	cpgpt1(bpi->phasevar(pctr.xvarid),
	       bpi->phasevar(pctr.yvarid),-1);
    }
}
#else
void plot(nbody * nbp)
{
}
#endif

void getbody (nbody * nbp)
{
    int n;
    cin >> n;
    nbp->n = n;
    nbp->bp = new body[n];
    body * bpi = nbp->bp;
    for(int i = 0; i<n;i++){
	cin >> *bpi;
	bpi++;
    }
}

const vector uniform_random_position_in_sphere(double power_index)
{
    vector x;
    do{
	for(int i = 0; i<VLEN;i++) x[i] = drand48()*2-1;
    }while (x*x >= 1);
    x *=  pow(x*x, 3.0/(power_index+3)-1);
    return x;
}


void create_uniform_sphere(nbody * nbp)
{
    cerr << "Enter n, power_index:";
    int n;
    double power_index;
    cin >> n >> power_index;
    PRC(n); PRL(power_index);
    nbp->n = n;
    nbp->bp = new body[n];
    body * bpi = nbp->bp;
    for(int i = 0; i<n;i++){
	bpi->set_pos(uniform_random_position_in_sphere(power_index));
	bpi->set_vel(0);
	bpi->set_mass(1.0/n);
	bpi++;
    }
}

void dump(nbody * nbp)
{
    cerr << "N = " << nbp->n << " eps2 = " << nbp->eps2 << "\n";
    int i;
    body * bpi;
    for(i=0, bpi= nbp->bp; i < nbp->n; i++,bpi++)bpi->dump();
}
    
main()
{
    nbody nb;
    nbody * nbp = &nb;
    double  t, toutstep,  tend, h, eps ;
    cerr << "Enter tend,toustep, h, eps:";
    cin >> tend >> toutstep >> h >> eps;
    nbp->eps2 = eps*eps;

    PRC(tend); PRC(h); PRL(nbp->eps2);
    static vintegfunc integtab[]={leapfrog,yoshida4,yoshida6b};
    static char * integname[] = {"leapfrog","yoshida4","yoshida6b"};
    static int  integsteps[] = {1,3,7};
    cerr << "Enter integrator type (0: leapfrog, 1: yoshida4, 2: yoshida6b:";
    int integrator_type;
    cin >> integrator_type;
    cerr << "Enter xmax:";
    double xmax;
    cin >> xmax;
    initgraph(xmax);
    cerr << "before getbody" ;  PRL(nbp->eps2);
    create_uniform_sphere(nbp);
    cerr << "after getbody" ;  PRL(nbp->eps2);
    plot(nbp);
    double e0, ke, pe;
    calc_force_for_system(nbp);
    calc_energy(nbp, ke, pe, e0);
    cerr << "Energies " << e0 << " " << ke << " " << pe <<"\n";
    plot(nbp);
    int first_call = 1;
    double e;
    for(t=0; t<tend; t+=h ){
	integtab[integrator_type](nbp,h,first_call);
	calc_energy(nbp,ke,pe,e);
	if (fmod(t,toutstep) == 0.0){
	    cerr << "T,e,de="<< t << " " <<  e << " " << e-e0 <<"\n";
	}
	plot(nbp);
    }
    int dum;
    cin >> dum;
#ifndef NOGRAPHICS    
    cpgend();
#endif
    calc_energy(nbp,ke,pe,e);
    cout << integname[integrator_type] << " " <<eps << " " << 1/h << " " << eps/h<< " " << eps/h*integsteps[integrator_type] << " "<<  e-e0 <<endl; 
}





