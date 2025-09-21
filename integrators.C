//
// integrators.C
//
// basic integrators for the vector class.
//

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
    for(int i = 0; i<VLEN;i+=2){
	x[1+i] += h*fprev[1+i]*0.5;
	x[i] += h*x[1+i];
    }
    fprev = f(x);
    for(int i = 0; i<VLEN;i+=2){
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
    }
    leapfrog(x,f,h*d1, first_call);
    leapfrog(x,f,h*d2, first_call);
    leapfrog(x,f,h*d1, first_call);
}
void yoshida6b(vector & x,
	      vfunc_ptr f,
	      double h,
	      int & first_call)
{
    const double d[4] = {0.784513610477560e0, 0.235573213359357e0,
			  -1.17767998417887e0, 1.31518632068391e0};
    for(int i=0;i<3;i++)leapfrog(x,f,h*d[i], first_call);
    leapfrog(x,f,h*d[3], first_call);
    for(int i=0;i<3;i++)leapfrog(x,f,h*d[2-i], first_call);
}


void variable_step_integrator(vector & x,
			      vintegfunc integ,
			      vfunc_ptr f,
			      double &h,
			      double hmax,
			      double tol,
			      double tolmin_factor )
{
    vector x0 = x;
    int first_call = 1;
    integ(x0,f,h,first_call);
    first_call = 1;
    integ(x,f,h/2,first_call);
    integ(x,f,h/2,first_call);

    double err = sqrt((x-x0)*(x-x0));

    if((err < tol*tolmin_factor) && (h < hmax)) h*= 2;
    if(err > tol) h /= 2;
}

